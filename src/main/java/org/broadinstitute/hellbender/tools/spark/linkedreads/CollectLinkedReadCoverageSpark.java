package org.broadinstitute.hellbender.tools.spark.linkedreads;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import org.apache.spark.api.java.JavaDoubleRDD;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.ReadMetadata;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.SVReadFilter;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalTree;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVKmerShort;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVKmerizer;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.seqdoop.hadoop_bam.util.NIOFileUtil;
import scala.Tuple2;

import java.io.*;
import java.nio.file.FileVisitResult;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.SimpleFileVisitor;
import java.nio.file.attribute.BasicFileAttributes;
import java.util.*;
import java.util.stream.Collectors;

@CommandLineProgramProperties(
        summary = "Computes the coverage by long molecules from linked-read data",
        oneLineSummary = "CollectLinkedReadCoverage on Spark",
        programGroup = SparkProgramGroup.class
)
public class CollectLinkedReadCoverageSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    @Argument(doc = "uri for the output file",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            optional = true)
    public String out;

    @Argument(doc = "cluster size",
            shortName = "clusterSize", fullName = "clusterSize",
            optional = true)
    public int clusterSize = 5000;

    @Argument(doc = "molSizeHistogramFile",
            shortName = "molSizeHistogramFile", fullName = "molSizeHistogramFile",
            optional = true)
    public String molSizeHistogramFile;

    @Argument(doc = "gapHistogramFile",
            shortName = "gapHistogramFile", fullName = "gapHistogramFile",
            optional = true)
    public String gapHistogramFile;

    @Argument(doc = "metadataFile",
            shortName = "metadataFile", fullName = "metadataFile",
            optional = true)
    public String metadataFile;

    @Argument(doc = "molSizeReadDensityFile",
            shortName = "molSizeReadDensityFile", fullName = "molSizeReadDensityFile",
            optional = true)
    public String molSizeReadDensityFile;

    @Argument(doc = "barcodeFragmentCountsFile",
            shortName = "barcodeFragmentCountsFile", fullName = "barcodeFragmentCountsFile",
            optional = true)
    public String barcodeFragmentCountsFile;

    @Argument(doc = "startMotifOutputCountFile",
            shortName = "startMotifOutputCountFile", fullName = "startMotifOutputCountFile",
            optional = true)
    public String startMotifOutputCountFile;

    @Argument(doc = "phaseSetIntervalsFile",
            shortName = "phaseSetIntervalsFile", fullName = "phaseSetIntervalsFile",
            optional = true)
    public String phaseSetIntervalsFile;

    @Argument(fullName = "minEntropy", shortName = "minEntropy", doc="Minimum trigram entropy of reads for filter", optional=true)
    public double minEntropy = 0;

    @Argument(fullName = "minReadCountPerMol", shortName = "minReadCountPerMol", doc="Minimum number of reads to call a molecule", optional=true)
    public int minReadCountPerMol = 2;

    @Argument(fullName = "edgeReadMapqThreshold", shortName = "edgeReadMapqThreshold", doc="Mapq below which to trim reads from starts and ends of molecules", optional=true)
    public int edgeReadMapqThreshold = 30;

    @Argument(fullName = "minMaxMapq", shortName = "minMaxMapq", doc="Minimum highest mapq read to create a fragment", optional=true)
    public int minMaxMapq = 30;

    @Argument(fullName = "molSizeDeviationFile", shortName = "molSizeDeviationFile", doc="file containing estiamtes of molecule size deviations across the genome", optional=true)
    public String molSizeDeviationFile;

    @Argument(fullName = "molSizeDeviationTestBinSize", shortName = "molSizeDeviationTestBinSize", doc="width of windows within which to sample mol sizes", optional=true)
    public int molSizeDeviationTestBinSize = 10000;

    private static final int REF_RECORD_LEN = 10000;
    // assuming we have ~1Gb/core, we can process ~1M kmers per partition
    private static final int REF_RECORDS_PER_PARTITION = 1024*1024 / REF_RECORD_LEN;


    @Override
    public boolean requiresReads() { return true; }

    @Override
    public boolean requiresReference() { return true; }

    @Override
    public ReadFilter makeReadFilter() {
        return super.makeReadFilter().
                and(new LinkedReadAnalysisFilter(minEntropy));
    }

    @Override
    protected void runTool(final JavaSparkContext ctx) {
        final JavaRDD<GATKRead> reads = getReads();

        final int finalClusterSize = clusterSize;

        final JavaRDD<GATKRead> mappedReads =
                reads.filter(read ->
                        !read.isDuplicate() && !read.failsVendorQualityCheck() && !read.isUnmapped() && ! read.isSecondaryAlignment() && ! read.isSupplementaryAlignment());
        final ReadMetadata readMetadata = new ReadMetadata(Collections.emptySet(), getHeaderForReads(), 2000,
                mappedReads, new SVReadFilter(new StructuralVariationDiscoveryArgumentCollection.FindBreakpointEvidenceSparkArgumentCollection()));

        if (metadataFile != null) {
            ReadMetadata.writeMetadata(readMetadata, metadataFile);
        }

        final Broadcast<ReadMetadata> broadcastMetadata = ctx.broadcast(readMetadata);
        final ReferenceMultiSource reference = getReference();
        final Broadcast<ReferenceMultiSource> broadcastReference = ctx.broadcast(reference);

        final JavaPairRDD<String, SVIntervalTree<List<ReadInfo>>> barcodeIntervals;
        if (out != null || barcodeFragmentCountsFile != null || molSizeHistogramFile != null || gapHistogramFile != null) {
            barcodeIntervals = getBarcodeIntervals(finalClusterSize, mappedReads, broadcastMetadata, minReadCountPerMol, minMaxMapq, edgeReadMapqThreshold);
        } else {
            barcodeIntervals = null;
        }

        if (barcodeFragmentCountsFile != null) {
            computeFragmentCounts(barcodeIntervals, barcodeFragmentCountsFile);
        }

        final Tuple2<double[], long[]> molSizeHistogram = molSizeHistogramFile != null ? computeMolSizeHistogram(barcodeIntervals, molSizeHistogramFile) : null;


        final Map<String, SVIntervalTree<Integer>> localBarcodeIntervals = barcodeIntervals.mapValues(oldTree -> {
            final SVIntervalTree<Integer> newTree = new SVIntervalTree<>();
            oldTree.forEach(e -> newTree.put(e.getInterval(), null));
            return newTree;
        }).collectAsMap();

        Broadcast<Map<String, SVIntervalTree<Integer>>> broadcastAllIntervals = ctx.broadcast(localBarcodeIntervals);

        JavaPairRDD<Integer, SVIntervalTree<String>> barcodeIntervalsByContig = barcodeIntervals.flatMapToPair(entry -> {
            final String barcode = entry._1();
            final SVIntervalTree<List<ReadInfo>> tree = entry._2();
            final Map<Integer, SVIntervalTree<String>> intervalsByChrom = new HashMap<>();
            tree.iterator().forEachRemaining(e -> {
                final int contig = e.getInterval().getContig();
                final SVInterval interval = e.getInterval();
                if (!intervalsByChrom.containsKey(contig)) intervalsByChrom.put(contig, new SVIntervalTree<>());
                intervalsByChrom.get(contig).put(interval, barcode);
            });
            return Utils.stream(intervalsByChrom.entrySet()).map(e -> new Tuple2<>(e.getKey(), e.getValue())).iterator();
        }).reduceByKey((tree1, tree2) -> {
            tree2.iterator().forEachRemaining(e -> tree1.put(e.getInterval(), e.getValue()));
            return tree1;
        }) .repartition(readMetadata.getContigNameMap().size());


        final JavaPairRDD<SVInterval, Set<String>> barcodeSetsByIntervals = barcodeIntervalsByContig.flatMapToPair(kv -> {
            final int contig = kv._1();
            final SVIntervalTree<String> intervals = kv._2();
            final List<Tuple2<Integer, String>> starts = new ArrayList<>(intervals.size());
            intervals.iterator().forEachRemaining(e -> starts.add(new Tuple2<>(e.getInterval().getStart(), e.getValue())));
            final List<Tuple2<Integer, String>> stops = new ArrayList<>(intervals.size());
            intervals.iterator().forEachRemaining(e -> stops.add(new Tuple2<>(e.getInterval().getEnd(), e.getValue())));
            stops.sort(Comparator.comparingInt(Tuple2::_1));

            int startIdx = 0;
            int stopIdx = 0;
            int prevBoundary = 0;

            final Set<String> currentBarcodes = new HashSet<>();

            final List<Tuple2<SVInterval, Set<String>>> results = new ArrayList<>(starts.size() + stops.size());
            while (stopIdx < stops.size()) {
                final int nextStart = startIdx < starts.size() ? starts.get(startIdx)._1 : -1;
                final int nextStop = stops.get(stopIdx)._1;
                if (nextStart != -1 && nextStart <= nextStop) {
                    // process a start
                    //System.out.println("Process a start at " + nextStart + ": " + starts.get(startIdx)._2);
                    if (!currentBarcodes.isEmpty()) {
                        final SVInterval newInterval = new SVInterval(contig, prevBoundary, nextStart);
                        final Set<String> barcodes = new HashSet<>(currentBarcodes);
                        results.add(new Tuple2<>(newInterval, barcodes));
                    }

                    currentBarcodes.add(starts.get(startIdx)._2);
                    prevBoundary = starts.get(startIdx)._1;
                    startIdx++;
                } else {
                    // process a stop
                    //System.out.println("Process a stop at " + nextStop + ": " + stops.get(stopIdx)._2);
                    final SVInterval newInterval = new SVInterval(contig, prevBoundary, nextStop);
                    final Set<String> barcodes = new HashSet<>(currentBarcodes);
                    results.add(new Tuple2<>(newInterval, barcodes));

                    currentBarcodes.remove(stops.get(stopIdx)._2);
                    prevBoundary = stops.get(stopIdx)._1;
                    stopIdx++;
                }

            }
            return results.iterator();
        });

        barcodeSetsByIntervals.saveAsTextFile("foo");

        JavaPairRDD<Tuple2<SVInterval, SVInterval>, Integer> allIntervalPairOverlaps = barcodeSetsByIntervals.flatMapToPair(kv -> {
            final SVInterval interval = kv._1();
            final Set<String> barcodes = kv._2();

            final SVIntervalTree<Integer> sharedBarcodeIntervals = new SVIntervalTree<>();
            barcodes.iterator().forEachRemaining(barcode -> {
                broadcastAllIntervals.getValue().get(barcode).iterator().forEachRemaining(barcodeIntervalEntry -> {
                    if (!barcodeIntervalEntry.getInterval().overlaps(interval)) {
                        sharedBarcodeIntervals.put(barcodeIntervalEntry.getInterval(), barcodeIntervalEntry.getValue());
                    }
                });
            });

            return new IntervalDepthIterator(interval, sharedBarcodeIntervals);

        });

        JavaPairRDD<Tuple2<SVInterval, SVInterval>, Integer> filteredPairedOverlaps = allIntervalPairOverlaps.filter(kv -> kv._2() > 1);
        filteredPairedOverlaps.saveAsTextFile("foobar");

        final SVIntervalTree<String> allIntervals = new SVIntervalTree<>();
        localBarcodeIntervals.forEach((barcode, value) -> {
            // todo: could be a bug if exact overlapping intervals
            value.forEach(intervalEntry -> allIntervals.put(intervalEntry.getInterval(), barcode));
        });




        //final Set<String> currentBarcodes = new HashSet<>();



//        if (molSizeReadDensityFile != null) {
//            barcodeIntervals.flatMapToPair(kv -> {
//                final SVIntervalTree<List<ReadInfo>> SVIntervalTree<List<ReadInfo>> = kv._2();
//                List<Double> results = new ArrayList<>(SVIntervalTree<List<ReadInfo>>.myTree.size());
//                final Iterator<SVIntervalTree.Entry<List<ReadInfo>>> iterator = SVIntervalTree<List<ReadInfo>>.myTree.iterator();
//                Utils.stream(iterator).map(e -> e.getInterval().getLength()).forEach(l -> results.add(new Tuple2<>(l % 1000)));
//            }
//        }
//        final List<Double> molSizeSample = moleculeSizes.takeSample(false, 100000);
//        final EmpiricalDistribution molSizeDistribution = new EmpiricalDistribution();
//        molSizeDistribution.load(molSizeSample.);
//        molSizeDistribution.

        if (gapHistogramFile != null) {
            computeGapSizeHistogram(barcodeIntervals, gapHistogramFile);
        }

        if (out != null) {
            writeIntervalsAsBed12(broadcastMetadata, barcodeIntervals, shardedOutput, out);
        }

        if (startMotifOutputCountFile != null) {
            computeStartMotifs(mappedReads, broadcastReference, startMotifOutputCountFile);

        }

        if (phaseSetIntervalsFile != null) {
            final JavaRDD<String> phaseSetIntervals = mappedReads.mapPartitions(iter -> {
                final ReadMetadata metadata = broadcastMetadata.getValue();
                List<Tuple2<Integer, SVInterval>> partitionPhaseSets = new ArrayList<>();
                int currentPS = -1;
                int currentContig = -1;
                int currentStart = -1;
                int currentEnd = -1;
                while (iter.hasNext()) {
                    GATKRead read = iter.next();
                    if (read.hasAttribute("PS")) {
                        int ps = read.getAttributeAsInteger("PS");
                        if (ps != currentPS) {
                            if (currentPS != -1) {
                                partitionPhaseSets.add(new Tuple2<>(currentPS, new SVInterval(currentContig, currentStart, currentEnd)));
                            }
                            currentPS = ps;
                            currentContig = metadata.getContigID(read.getContig());
                            currentStart = read.getStart();
                            currentEnd = read.getEnd();
                        } else {
                            currentEnd = read.getEnd();
                        }
                    }
                }
                return partitionPhaseSets.iterator();
            }).mapToPair(v -> v).reduceByKey((svInterval1, svInterval2) -> new SVInterval(svInterval1.getContig(),
                    Math.min(svInterval1.getStart(), svInterval2.getStart()),
                    Math.max(svInterval1.getEnd(), svInterval2.getEnd())))
                    .mapToPair(kv -> new Tuple2<>(kv._2(), kv._1())).sortByKey()
                    .map(kv -> {
                        final SVInterval svInterval = kv._1();
                        final Integer ps = kv._2();
                        final String contigName = lookupContigName(svInterval.getContig(), broadcastMetadata.getValue());
                        return contigName + "\t" + svInterval.getStart() + "\t" + svInterval.getEnd() + "\t" + ps;
                    });


            if (shardedOutput) {
                phaseSetIntervals.saveAsTextFile(phaseSetIntervalsFile);
            } else {
                final String shardedOutputDirectory = phaseSetIntervalsFile + ".parts";
                phaseSetIntervals.saveAsTextFile(shardedOutputDirectory);
                unshardOutput(phaseSetIntervalsFile, shardedOutputDirectory, phaseSetIntervals.getNumPartitions());
            }
        }

        final int molSizeDeviationTestBinSizeFinal = molSizeDeviationTestBinSize;

        if (molSizeHistogram != null && molSizeDeviationFile != null) {
            final JavaPairRDD<SVInterval, Integer> molSizesBySamplePoints = barcodeIntervals.flatMapToPair(t -> {
                        final List<Tuple2<SVInterval, Integer>> results = new ArrayList<>();
                        final Iterator<SVIntervalTree.Entry<List<ReadInfo>>> iterator = t._2().iterator();
                        while (iterator.hasNext()) {
                            final SVIntervalTree.Entry<List<ReadInfo>> next = iterator.next();
                            final List<ReadInfo> readInfos = next.getValue();
                            readInfos.sort((o1, o2) -> new Integer(o1.getStart()).compareTo(o2.getStart()));
                            final int contig = readInfos.get(0).getContig();
                            final int minStart = readInfos.get(0).getStart();
                            final int maxEnd = readInfos.get(0).getStart();
                            final int length = maxEnd - minStart;

                            int sample = minStart + minStart % molSizeDeviationTestBinSizeFinal;
                            while (sample < maxEnd) {
                                results.add(new Tuple2<>(new SVInterval(contig, sample, sample), length));
                                sample += molSizeDeviationTestBinSizeFinal;
                            }
                        }
                        return results.iterator();
                    }
            );
            //molSizesBySamplePoints.
        }
    }

    private static void computeStartMotifs(final JavaRDD<GATKRead> mappedReads, final Broadcast<ReferenceMultiSource> broadcastReference, final String startMotifOutputCountFile) {

        final Map<String, Long> countByMotif = mappedReads.filter(GATKRead::isFirstOfPair).flatMap(r -> {

            final ReferenceMultiSource ref = broadcastReference.getValue();
            final SAMSequenceDictionary referenceSequenceDictionary = ref.getReferenceSequenceDictionary(null);
            final SAMSequenceRecord sequence = referenceSequenceDictionary.getSequence(r.getContig());
            if (sequence == null) {
                return Collections.emptyListIterator();
            }

            if (r.isReverseStrand() && (r.getUnclippedEnd() < 1 || r.getUnclippedEnd() > sequence.getSequenceLength() - 7)) {
                return Collections.emptyListIterator();
            }

            if (!r.isReverseStrand() && (r.getUnclippedStart() < 8 || r.getUnclippedStart() > sequence.getSequenceLength() - 7)) {
                return Collections.emptyListIterator();
            }

            final SimpleInterval motifInterval = new SimpleInterval(r.getContig(), r.isReverseStrand() ? r.getUnclippedEnd() + 1 : r.getUnclippedStart() - 7,
                    r.isReverseStrand() ? r.getUnclippedEnd() + 7 : r.getUnclippedStart() - 1 );
            if (motifInterval.getStart() < 1) return Collections.emptyListIterator();

            if (motifInterval.getEnd() >= sequence.getSequenceLength()) {
                return Collections.emptyListIterator();
            }

            final ReferenceBases referenceBases = ref.getReferenceBases(null, motifInterval);
            final String baseString = new String(referenceBases.getBases());
            if (baseString.contains("N")) return Collections.emptyListIterator();

            final String kmer = SVKmerizer.toKmer(baseString, new SVKmerShort(7)).canonical(7).toString(7);

            return Collections.singletonList(kmer).iterator();

        }).countByValue();

        try (final Writer writer =
                     new BufferedWriter(new OutputStreamWriter(BucketUtils.createFile(startMotifOutputCountFile)))) {
            writer.write("# motif counts\n");
            for (String motif : countByMotif.keySet()) {
                writer.write(motif + "\t" + countByMotif.get(motif) + "\n");
            }
        } catch (final IOException ioe) {
            throw new GATKException("Can't start motif file.", ioe);
        }
    }

    private static void writeIntervalsAsBed12(final Broadcast<ReadMetadata> broadcastMetadata, final JavaPairRDD<String, SVIntervalTree<List<ReadInfo>>> barcodeIntervals, final boolean shardedOutput, final String out) {
        final JavaPairRDD<SVInterval, String> bedRecordsByBarcode;
        bedRecordsByBarcode = barcodeIntervals.flatMapToPair(x -> {
            final String barcode = x._1;
            final SVIntervalTree<List<ReadInfo>> svIntervalTree = x._2;

            final List<Tuple2<SVInterval, String>> results = new ArrayList<>();
            for (final SVIntervalTree.Entry<List<ReadInfo>> next : svIntervalTree) {
                results.add(new Tuple2<>(next.getInterval(), intervalTreeToBedRecord(barcode, next, broadcastMetadata.getValue())));
            }

            return results.iterator();
        });

        if (shardedOutput) {
            bedRecordsByBarcode.values().saveAsTextFile(out);
        } else {
            final String shardedOutputDirectory = out + ".parts";
            final int numParts = bedRecordsByBarcode.getNumPartitions();
            bedRecordsByBarcode.sortByKey().values().saveAsTextFile(shardedOutputDirectory);
            //final BlockCompressedOutputStream outputStream = new BlockCompressedOutputStream(out);
            unshardOutput(out, shardedOutputDirectory, numParts);
        }
    }

    private static void unshardOutput(final String out, final String shardedOutputDirectory, final int numParts) {
        final OutputStream outputStream;

        outputStream = new BlockCompressedOutputStream(new BufferedOutputStream(BucketUtils.createFile(out)), null);
        for (int i = 0; i < numParts; i++) {
            String fileName = String.format("part-%1$05d", i);
            try {
                final BufferedInputStream bufferedInputStream = new BufferedInputStream(Files.newInputStream(NIOFileUtil.asPath(shardedOutputDirectory + System.getProperty("file.separator") + fileName)));
                int bite;
                while ((bite = bufferedInputStream.read()) != -1) {
                    outputStream.write(bite);
                }
                bufferedInputStream.close();
            } catch (IOException e) {
                throw new GATKException(e.getMessage());
            }
        }
        try {
            outputStream.close();
        } catch (IOException e) {
            throw new GATKException(e.getMessage());
        }
        try {
            deleteRecursive(NIOFileUtil.asPath(shardedOutputDirectory));
        } catch (IOException e) {
            throw new GATKException(e.getMessage());
        }
    }

    private static void computeGapSizeHistogram(final JavaPairRDD<String, SVIntervalTree<List<ReadInfo>>> barcodeIntervals, final String gapHistogramFile) {
        final Tuple2<double[], long[]> readGapHistogram = barcodeIntervals.flatMapToDouble(kv -> {
            final SVIntervalTree<List<ReadInfo>> intervalTree = kv._2();
            List<Double> results = new ArrayList<>();
            for (final SVIntervalTree.Entry<List<ReadInfo>> next : intervalTree) {
                List<ReadInfo> readInfoList = next.getValue();
                readInfoList.sort(Comparator.comparing(ReadInfo::getStart));
                for (int i = 1; i < readInfoList.size(); i++) {
                    results.add((double) (readInfoList.get(i).start - readInfoList.get(i - 1).start));
                }
            }
            return results.iterator();
        }).histogram(1000);


        try (final Writer writer =
                     new BufferedWriter(new OutputStreamWriter(BucketUtils.createFile(gapHistogramFile)))) {
            writer.write("# Read gap histogram\n");
            for (int i = 1; i < readGapHistogram._1().length; i++) {
                writer.write(readGapHistogram._1()[i - 1] + "-" + readGapHistogram._1()[i] + "\t" + readGapHistogram._2()[i - 1] + "\n");
            }
        } catch (final IOException ioe) {
            throw new GATKException("Can't write gap histogram file.", ioe);
        }
    }

    private static Tuple2<double[], long[]> computeMolSizeHistogram(final JavaPairRDD<String, SVIntervalTree<List<ReadInfo>>> barcodeIntervals, final String molSizeHistogramFile) {
        final JavaDoubleRDD moleculeSizes = barcodeIntervals.flatMapToDouble(kv -> {
            final SVIntervalTree<List<ReadInfo>> intervalTree = kv._2();
            List<Double> results = new ArrayList<>(intervalTree.size());
            final Iterator<SVIntervalTree.Entry<List<ReadInfo>>> iterator = intervalTree.iterator();
            Utils.stream(iterator).map(e -> e.getInterval().getLength()).forEach(l -> results.add(new Double(l)));
            return results.iterator();
        });

        final Tuple2<double[], long[]> moleculeLengthHistogram = moleculeSizes.histogram(1000);

        try (final Writer writer =
                     new BufferedWriter(new OutputStreamWriter(BucketUtils.createFile(molSizeHistogramFile)))) {
            writer.write("# Molecule length histogram\n");
            for (int i = 1; i < moleculeLengthHistogram._1().length; i++) {
                writer.write(moleculeLengthHistogram._1()[i - 1] + "-" + moleculeLengthHistogram._1()[i] + "\t" + moleculeLengthHistogram._2()[i - 1] + "\n");
            }
        } catch (final IOException ioe) {
            throw new GATKException("Can't write read histogram file.", ioe);
        }

        return moleculeLengthHistogram;
    }

    private static void computeFragmentCounts(final JavaPairRDD<String, SVIntervalTree<List<ReadInfo>>> barcodeIntervals, final String barcodeFragmentCountsFile) {
        final JavaPairRDD<Integer, String> barcodeReadCounts = barcodeIntervals.mapToPair(kv ->
                new Tuple2<>(kv._1(), Utils.stream(kv._2().iterator()).mapToInt(e -> e.getValue().size()).sum())).mapToPair(Tuple2::swap).sortByKey(true, 1);
        barcodeReadCounts.saveAsTextFile(barcodeFragmentCountsFile);
    }

    private static JavaPairRDD<String, SVIntervalTree<List<ReadInfo>>> getBarcodeIntervals(final int finalClusterSize, final JavaRDD<GATKRead> mappedReads, final Broadcast<ReadMetadata> broadcastMetadata, final int minReadCountPerMol, final int minMaxMapq, final int edgeReadMapqThreshold) {
        final JavaPairRDD<String, SVIntervalTree<List<ReadInfo>>> readsClusteredByBC = getClusteredReadIntervalsByTag(finalClusterSize, mappedReads, broadcastMetadata, "BX");
        final JavaPairRDD<String, SVIntervalTree<List<ReadInfo>>> edgefilteredClusteredReads = edgeFilterFragments(readsClusteredByBC, edgeReadMapqThreshold);
        final JavaPairRDD<String, SVIntervalTree<List<ReadInfo>>> minMoleculeCountFragments = minMoleculeCountFilterFragments(edgefilteredClusteredReads, minReadCountPerMol);
        final JavaPairRDD<String, SVIntervalTree<List<ReadInfo>>> maxMapqFilteredFragments = minMaxMapqFilterFragments(minMoleculeCountFragments, minMaxMapq);
        return maxMapqFilteredFragments
                .cache();
    }

    private static JavaPairRDD<String, SVIntervalTree<List<ReadInfo>>> minMoleculeCountFilterFragments(final JavaPairRDD<String, SVIntervalTree<List<ReadInfo>>> clusteredReads, final int minReadCountPerMol) {
        if (minReadCountPerMol > 0) {
            final JavaPairRDD<String, SVIntervalTree<List<ReadInfo>>> filteredIntervalsByKey = clusteredReads
                    .mapValues(intervalTree -> {
                        final Iterator<SVIntervalTree.Entry<List<ReadInfo>>> iterator = intervalTree.iterator();
                        while (iterator.hasNext()) {
                            final SVIntervalTree.Entry<List<ReadInfo>> next = iterator.next();
                            if (next.getValue().size() < minReadCountPerMol) {
                                iterator.remove();
                            }
                        }
                        return intervalTree;
                    })
                    .filter(kv -> {
                        final SVIntervalTree<List<ReadInfo>> intervalTree = kv._2();
                        return intervalTree.size() > 0;
                    });
            return filteredIntervalsByKey;
        } else {
            return clusteredReads;
        }
    }

    private static JavaPairRDD<String, SVIntervalTree<List<ReadInfo>>> minMaxMapqFilterFragments(final JavaPairRDD<String, SVIntervalTree<List<ReadInfo>>> readsClusteredByBC, final int minMaxMapq) {
        if (minMaxMapq > 0) {
            return readsClusteredByBC.flatMapToPair(kv -> {
                final SVIntervalTree<List<ReadInfo>> intervalTree = kv._2();
                final Iterator<SVIntervalTree.Entry<List<ReadInfo>>> treeIterator = intervalTree.iterator();
                while (treeIterator.hasNext()) {
                    SVIntervalTree.Entry<List<ReadInfo>> entry = treeIterator.next();
                    final List<ReadInfo> value = entry.getValue();
                    if (!value.stream().anyMatch(readInfo -> readInfo.getMapq() >= minMaxMapq)) {
                        treeIterator.remove();
                    }
                }
                if (intervalTree.size() == 0) {
                    return Collections.emptyIterator();
                } else {
                    return Collections.singletonList(new Tuple2<>(kv._1(), intervalTree)).iterator();
                }
            });
        } else {
            return readsClusteredByBC;
        }

    }

    private static JavaPairRDD<String, SVIntervalTree<List<ReadInfo>>> edgeFilterFragments(final JavaPairRDD<String, SVIntervalTree<List<ReadInfo>>> readsClusteredByBC, final int edgeReadMapqThreshold) {
        if (edgeReadMapqThreshold > 0) {
            return readsClusteredByBC.flatMapToPair(kv -> {
                final SVIntervalTree<List<ReadInfo>> myTree = kv._2();
                final Iterator<SVIntervalTree.Entry<List<ReadInfo>>> treeIterator = myTree.iterator();
                while (treeIterator.hasNext()) {
                    SVIntervalTree.Entry<List<ReadInfo>> entry =  treeIterator.next();
                    final List<ReadInfo> value = entry.getValue();
                    final Iterator<ReadInfo> iterator = value.iterator();
                    while (iterator.hasNext()) {
                        final ReadInfo readInfo = iterator.next();
                        if (readInfo.getMapq() > edgeReadMapqThreshold) {
                            break;
                        } else {
                            iterator.remove();
                        }
                    }
                    final ListIterator<ReadInfo> revIterator = value.listIterator(value.size());
                    while (revIterator.hasPrevious()) {
                        final ReadInfo readInfo = revIterator.previous();
                        if (readInfo.getMapq() > edgeReadMapqThreshold) {
                            break;
                        } else {
                            revIterator.remove();
                        }
                    }
                    if (entry.getValue().isEmpty()) {
                        treeIterator.remove();
                    }
                }
                if (myTree.size() == 0) {
                    return Collections.emptyIterator();
                } else {
                    return Collections.singletonList(new Tuple2<>(kv._1(), myTree)).iterator();
                }
            });
        } else {
            return readsClusteredByBC;
        }

    }

    private static JavaPairRDD<String, SVIntervalTree<List<ReadInfo>>> getPhaseSetIntervals(final int finalClusterSize, final JavaRDD<GATKRead> mappedReads, final Broadcast<ReadMetadata> broadcastMetadata) {
        return getClusteredReadIntervalsByTag(1000000, mappedReads, broadcastMetadata, "PS");
    }

    private static JavaPairRDD<String, SVIntervalTree<List<ReadInfo>>> getClusteredReadIntervalsByTag(final int finalClusterSize, final JavaRDD<GATKRead> mappedReads, final Broadcast<ReadMetadata> broadcastMetadata, final String tag) {

        final JavaPairRDD<String, SVIntervalTree<List<ReadInfo>>> intervalsByKey = mappedReads.filter(GATKRead::isFirstOfPair)
                .mapToPair(read -> new Tuple2<>(read.getAttributeAsString(tag), new ReadInfo(broadcastMetadata.getValue(), read)))
                .combineByKey(
                        readInfo -> {
                            SVIntervalTree<List<ReadInfo>> intervalTree = new SVIntervalTree<>();
                            return addReadToIntervals(intervalTree, readInfo, finalClusterSize);
                        }
                        ,
                        (aggregator, read) -> addReadToIntervals(aggregator, read, finalClusterSize),
                        (intervalTree1, intervalTree2) -> combineIntervalLists(intervalTree1, intervalTree2, finalClusterSize)
                );

        return intervalsByKey;
    }

    /**
     * Delete the given directory and all of its contents if non-empty.
     * @param directory the directory to delete
     */
    private static void deleteRecursive(Path directory) throws IOException {
        Files.walkFileTree(directory, new SimpleFileVisitor<Path>() {
            @Override
            public FileVisitResult visitFile(Path file, BasicFileAttributes attrs) throws IOException {
                Files.delete(file);
                return FileVisitResult.CONTINUE;
            }
            @Override
            public FileVisitResult postVisitDirectory(Path dir, IOException exc) throws IOException {
                Files.delete(dir);
                return FileVisitResult.CONTINUE;
            }
        });
    }


    @VisibleForTesting
    static SVIntervalTree<List<ReadInfo>> addReadToIntervals(final SVIntervalTree<List<ReadInfo>> inputIntervalTree, final ReadInfo read, final int clusterSize) {
        final SVInterval sloppedReadInterval = new SVInterval(read.getContig(), read.getStart() - clusterSize, read.getEnd()+clusterSize);
        SVIntervalTree<List<ReadInfo>> myTree = inputIntervalTree;
        if (inputIntervalTree == null) {
            myTree = new SVIntervalTree<>();
        }
        final Iterator<SVIntervalTree.Entry<List<ReadInfo>>> iterator = myTree.overlappers(
                sloppedReadInterval);
        int start = read.getStart();
        int end = read.getEnd();
        final List<ReadInfo> newReadList = new ArrayList<>();
        newReadList.add(read);
        if (iterator.hasNext()) {
            final SVIntervalTree.Entry<List<ReadInfo>> existingNode = iterator.next();
            final int currentStart = existingNode.getInterval().getStart();
            final int currentEnd = existingNode.getInterval().getEnd();
            final List<ReadInfo> currentValue = existingNode.getValue();
            start = Math.min(currentStart, read.getStart());
            end = Math.max(currentEnd, read.getEnd());
            newReadList.addAll(currentValue);
            iterator.remove();
        }
        while (iterator.hasNext()) {
            final SVIntervalTree.Entry<List<ReadInfo>> next = iterator.next();
            final int currentEnd = next.getInterval().getStart();
            final List<ReadInfo> currentValue = next.getValue();
            end = Math.max(end, currentEnd);
            newReadList.addAll(currentValue);
            iterator.remove();
        }
        myTree.put(new SVInterval(read.getContig(), start, end), newReadList);
        return myTree;
    }

    private static SVIntervalTree<List<ReadInfo>> combineIntervalLists(final SVIntervalTree<List<ReadInfo>> intervalTree1,
                                                                                      final SVIntervalTree<List<ReadInfo>> intervalTree2,
                                                                                      final int clusterSize) {
        return mergeIntervalTrees(intervalTree1, intervalTree2, clusterSize);

    }

    @VisibleForTesting
    static SVIntervalTree<List<ReadInfo>> mergeIntervalTrees(final SVIntervalTree<List<ReadInfo>> tree1, final SVIntervalTree<List<ReadInfo>> tree2, final int clusterSize) {

        if (tree1 == null || tree1.size() == 0) return tree2;
        if (tree2 == null || tree2.size() == 0) return tree1;

        final SVIntervalTree<List<ReadInfo>> mergedTree = new SVIntervalTree<>();

        PriorityQueue<SVIntervalTree.Entry<List<ReadInfo>>> nodes = new PriorityQueue<>(tree1.size() + tree2.size(),
                (o1, o2) -> {
                    if (o1.getInterval().getContig() == o2.getInterval().getContig()) {
                        return new Integer(o1.getInterval().getStart()).compareTo(o2.getInterval().getStart());
                    } else {
                        return new Integer(o1.getInterval().getContig()).compareTo(o2.getInterval().getContig());
                    }
                });

        tree1.iterator().forEachRemaining(nodes::add);
        tree2.iterator().forEachRemaining(nodes::add);

        int currentContig = -1;
        int currentStart = -1;
        int currentEnd = -1;
        List<ReadInfo> values = new ArrayList<>();
        while (nodes.size() > 0) {
            final SVIntervalTree.Entry<List<ReadInfo>> next = nodes.poll();
            final SVInterval newInterval = next.getInterval();
            final int newContig = newInterval.getContig();
            if (currentContig != newContig) {
                if (currentContig != -1) {
                    mergedTree.put(new SVInterval(currentContig, currentStart, currentEnd), values);
                }
                currentContig = newContig;
                currentStart = newInterval.getStart();
                currentEnd = newInterval.getEnd();
                values = new ArrayList<>(next.getValue());
            } else {
                if (overlaps(newInterval, currentEnd, clusterSize)) {
                    currentEnd = Math.max(currentEnd, newInterval.getEnd());
                    values.addAll(next.getValue());
                } else {
                    // no overlap, so put the previous node in and set up the next set of values
                    mergedTree.put(new SVInterval(currentContig, currentStart, currentEnd), values);
                    currentStart = newInterval.getStart();
                    currentEnd = newInterval.getEnd();
                    values = new ArrayList<>(next.getValue());
                }
            }

        }
        if (currentStart != -1) {
            mergedTree.put(new SVInterval(currentContig, currentStart, currentEnd), values);
        }

        return mergedTree;
    }

    private static boolean overlaps(final SVInterval newInterval, final int currentEnd, final int clusterSize) {
        return currentEnd + clusterSize > newInterval.getStart() - clusterSize;
    }

    static String intervalTreeToBedRecord(final String barcode, final SVIntervalTree.Entry<List<ReadInfo>> node, final ReadMetadata readMetadata) {
        final StringBuilder builder = new StringBuilder();
        builder.append(lookupContigName(node.getInterval().getContig(), readMetadata));
        builder.append("\t");
        builder.append(node.getInterval().getStart());
        builder.append("\t");
        builder.append(node.getInterval().getEnd());
        builder.append("\t");
        builder.append(barcode);
        builder.append("\t");
        builder.append(node.getValue().size());
        builder.append("\t");
        builder.append("+");
        builder.append("\t");
        builder.append(node.getInterval().getStart());
        builder.append("\t");
        builder.append(node.getInterval().getEnd());
        builder.append("\t");
        builder.append("0,0,255");
        builder.append("\t");
        builder.append(node.getValue().size());
        final List<ReadInfo> reads = node.getValue();
        reads.sort((o1, o2) -> new Integer(o1.getStart()).compareTo(o2.getStart()));
        builder.append("\t");
        builder.append(reads.stream().map(r -> String.valueOf(r.getEnd() - r.getStart() + 1)).collect(Collectors.joining(",")));
        builder.append("\t");
        builder.append(reads.stream().map(r -> String.valueOf(r.getStart() - node.getInterval().getStart())).collect(Collectors.joining(",")));
        builder.append("\t");
        builder.append(reads.stream().mapToInt(ReadInfo::getMapq).max().orElse(-1));
        return builder.toString();
    }

    private static String lookupContigName(final int contig, final ReadMetadata readMetadata) {
        for (Map.Entry<String, Integer> entry : readMetadata.getContigNameMap().entrySet()) {
            if (entry.getValue().equals(contig)) return entry.getKey();
        }
        throw new GATKException("Invalid contig ID: "+ contig);
    }

    @DefaultSerializer(ReadInfo.Serializer.class)
    static class ReadInfo {
        ReadInfo(final ReadMetadata readMetadata, final GATKRead gatkRead) {
            this.contig = readMetadata.getContigID(gatkRead.getContig());
            this.start = gatkRead.getStart();
            this.end = gatkRead.getEnd();
            this.forward = !gatkRead.isReverseStrand();
            this.mapq = gatkRead.getMappingQuality();
        }

        ReadInfo(final int contig, final int start, final int end, final boolean forward, final int mapq) {
            this.contig = contig;
            this.start = start;
            this.end = end;
            this.forward = forward;
            this.mapq = mapq;
        }

        int contig;
        int start;
        int end;
        boolean forward;
        int mapq;

        public ReadInfo(final Kryo kryo, final Input input) {
            contig = input.readInt();
            start = input.readInt();
            end = input.readInt();
            forward = input.readBoolean();
            mapq = input.readInt();
        }

        public int getContig() {
            return contig;
        }

        public int getStart() {
            return start;
        }

        public int getEnd() {
            return end;
        }

        public boolean isForward() {
            return forward;
        }

        public int getMapq() {
            return mapq;
        }

        public static final class Serializer extends com.esotericsoftware.kryo.Serializer<ReadInfo> {
            @Override
            public void write( final Kryo kryo, final Output output, final ReadInfo interval ) {
                interval.serialize(kryo, output);
            }

            @Override
            public ReadInfo read(final Kryo kryo, final Input input, final Class<ReadInfo> klass ) {
                return new ReadInfo(kryo, input);
            }
        }

        private void serialize(final Kryo kryo, final Output output) {
            output.writeInt(contig);
            output.writeInt(start);
            output.writeInt(end);
            output.writeBoolean(forward);
            output.writeInt(mapq);
        }

    }

    private class IntervalDepthIterator implements Iterator<Tuple2<Tuple2<SVInterval, SVInterval>, Integer>> {
        private final SVInterval interval;
        private final Iterator<Integer> contigIterator;
        private int currentDepth;
        private int currentContig;

        private int startIdx;
        private int stopIdx;
        private int prevBoundary;

        private final Map<Integer, List<Integer>> stops;
        private final Map<Integer, List<Integer>> starts;

        public IntervalDepthIterator(final SVInterval interval, final SVIntervalTree<Integer> sharedBarcodeIntervals) {
            //System.out.println("DepthIterator for " + interval);
            this.interval = interval;

            starts = new HashMap<>();
            stops = new HashMap<>();

            if (! (sharedBarcodeIntervals.size() == 0)) {
                sharedBarcodeIntervals.iterator().forEachRemaining(e -> {
                    final int contig = e.getInterval().getContig();
                    if (! starts.containsKey(contig)) starts.put(contig, new ArrayList<>());
                    starts.get(contig).add(e.getInterval().getStart());
                    if (! stops.containsKey(contig)) stops.put(contig, new ArrayList<>());
                    stops.get(contig).add(e.getInterval().getEnd());
                });

                for (Integer contig : stops.keySet()) {
                    Collections.sort(stops.get(contig));
                }

                contigIterator = starts.keySet().iterator();
                currentContig = contigIterator.next();

                startIdx = 1;
                stopIdx = 0;
                prevBoundary = starts.get(currentContig).get(0);
                currentDepth = 1;

            } else {
                contigIterator = null;
            }

        }

        @Override
        public boolean hasNext() {
            return contigIterator != null && (stopIdx < stops.get(currentContig).size() || contigIterator.hasNext());
        }

        @Override
        public Tuple2<Tuple2<SVInterval, SVInterval>, Integer> next() {
            if (stopIdx == stops.get(currentContig).size()) {
                currentContig = contigIterator.next();

                startIdx = 1;
                stopIdx = 0;
                prevBoundary = starts.get(currentContig).get(0);
                currentDepth = 1;

            }

            final int nextStart = startIdx < starts.get(currentContig).size() ? starts.get(currentContig).get(startIdx) : -1;
            final int nextStop = stops.get(currentContig).get(stopIdx);

            Tuple2<Tuple2<SVInterval, SVInterval>, Integer> result;
            if (nextStart != -1 && nextStart <= nextStop) {
                // process a start
                //System.out.println("Process a start at " + nextStart + ": " + starts.get(startIdx));

                final SVInterval newInterval = new SVInterval(currentContig, prevBoundary, nextStart);
                result = new Tuple2<>(new Tuple2<>(interval, newInterval), currentDepth);

                currentDepth++;
                prevBoundary = starts.get(currentContig).get(startIdx);
                startIdx++;

            } else {
                // process a stop
                //System.out.println("Process a stop at " + nextStop + ": " + stops.get(stopIdx));

                currentDepth--;
                final SVInterval newInterval = new SVInterval(currentContig, prevBoundary, nextStart);
                result = new Tuple2<>(new Tuple2<>(interval, newInterval), currentDepth);


                currentDepth--;
                prevBoundary = stops.get(currentContig).get(stopIdx);
                stopIdx++;
            }

            return result;
        }
    }
}

