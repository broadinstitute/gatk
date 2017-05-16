package org.broadinstitute.hellbender.tools.spark.linkedreads;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.*;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.Locatable;
import org.apache.commons.io.output.ByteArrayOutputStream;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.ReadMetadata;
import org.broadinstitute.hellbender.tools.spark.sv.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.SVIntervalTree;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
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

    @Argument(doc = "writeSAM",
            shortName = "writeSAM", fullName = "writeSAM",
            optional = true)
    public boolean writeSAM = true;

    @Argument(doc = "metadataFile",
            shortName = "metadataFile", fullName = "metadataFile",
            optional = true)
    public File metadataFile;

    @Argument(fullName = "minEntropy", shortName = "minEntropy", doc="Minimum trigram entropy of reads for filter", optional=true)
    public double minEntropy = 4.5;

    @Override
    public boolean requiresReads() { return true; }

    @Override
    public ReadFilter makeReadFilter() {

        return super.makeReadFilter().
                and(new LinkedReadAnalysisFilter(minEntropy));
    }

    @Override
    protected void runTool(final JavaSparkContext ctx) {
        final JavaRDD<GATKRead> reads = getReads();
        final SAMFileHeader headerForReads = getHeaderForReads();

        final int finalClusterSize = clusterSize;

        final JavaRDD<GATKRead> mappedReads =
                reads.filter(read ->
                        !read.isDuplicate() && !read.failsVendorQualityCheck() && !read.isUnmapped());
        final ReadMetadata readMetadata = new ReadMetadata(Collections.emptySet(), getHeaderForReads(), mappedReads);
        ReadMetadata.writeMetadata(readMetadata, metadataFile.getAbsolutePath());

        final Broadcast<ReadMetadata> broadcastMetadata = ctx.broadcast(readMetadata);

        if (writeSAM) {
            final JavaPairRDD<String, SVIntervalTree<List<GATKRead>>> barcodeIntervals =
                    reads.mapToPair(read -> new Tuple2<>(read.getAttributeAsString("BX"), read))
                            .aggregateByKey(
                                    new SVIntervalTree<List<GATKRead>>(),
                                    (aggregator, read) -> addReadToIntervals(aggregator, read, finalClusterSize, broadcastMetadata.getValue()),
                                    (intervalTree1, intervalTree2) -> combineIntervalLists(intervalTree1, intervalTree2, finalClusterSize)
                            );

            final JavaRDD<GATKRead> intervalsByBarcode =
                    barcodeIntervals.flatMap(x -> {
                        final String barcode = x._1;
                        final SVIntervalTree<List<GATKRead>> intervalTree = x._2;
                        final List<GATKRead> results = new ArrayList<>();
                        for (final SVIntervalTree.Entry<List<GATKRead>> next : intervalTree) {
                            results.add(intervalTreeToGATKRead(barcode, headerForReads, next, broadcastMetadata.getValue()));
                        }
                        return results.iterator();
                    });
            writeReads(ctx, out, intervalsByBarcode);

        } else {
            final JavaPairRDD<String, SVIntervalTree<List<ReadInfo>>> barcodeIntervals =
                    reads.filter(GATKRead::isFirstOfPair)
                            .mapToPair(read -> new Tuple2<>(read.getAttributeAsString("BX"), new ReadInfo(read.getContig(), read.getStart(), read.getEnd())))
                            .aggregateByKey(
                                    new SVIntervalTree<>(),
                                    (aggregator, read) -> addReadToIntervals(aggregator, read, finalClusterSize, broadcastMetadata.getValue()),
                                    (intervalTree1, intervalTree2) -> combineIntervalLists(intervalTree1, intervalTree2, finalClusterSize)
                            );

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
                final String shardedOutputDirectory = this.out + ".parts";
                final int numParts = bedRecordsByBarcode.getNumPartitions();
                bedRecordsByBarcode.sortByKey().values().saveAsTextFile(shardedOutputDirectory);
                //final BlockCompressedOutputStream outputStream = new BlockCompressedOutputStream(out);
                final OutputStream outputStream;

                outputStream = new BlockCompressedOutputStream(out);
                for (int i = 0; i < numParts; i++) {
                    String fileName = String.format("part-%1$05d", i);
                    try {
                        final BufferedInputStream bufferedInputStream = new BufferedInputStream(new FileInputStream(shardedOutputDirectory + System.getProperty("file.separator") + fileName));
                        int bite;
                        while ((bite= bufferedInputStream.read()) != -1) {
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
        }

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
    static <T extends Locatable> SVIntervalTree<List<T>> addReadToIntervals(final SVIntervalTree<List<T>> intervalList, final T read, final int clusterSize, final ReadMetadata readMetadata) {
        final SVInterval sloppedReadInterval = new SVInterval(readMetadata.getContigID(read.getContig()), read.getStart() - clusterSize, read.getEnd()+clusterSize);
        final Iterator<SVIntervalTree.Entry<List<T>>> iterator = intervalList.overlappers(
                sloppedReadInterval);
        int start = read.getStart();
        int end = read.getEnd();
        final List<T> newReadList = new ArrayList<>();
        newReadList.add(read);
        if (iterator.hasNext()) {
            final SVIntervalTree.Entry<List<T>> existingNode = iterator.next();
            final int currentStart = existingNode.getInterval().getStart();
            final int currentEnd = existingNode.getInterval().getEnd();
            final List<T> currentValue = existingNode.getValue();
            start = Math.min(currentStart, read.getStart());
            end = Math.max(currentEnd, read.getEnd());
            newReadList.addAll(currentValue);
            iterator.remove();
        }
        while (iterator.hasNext()) {
            final SVIntervalTree.Entry<List<T>> next = iterator.next();
            final int currentEnd = next.getInterval().getStart();
            final List<T> currentValue = next.getValue();
            end = Math.max(end, currentEnd);
            newReadList.addAll(currentValue);
            iterator.remove();
        }
        intervalList.put(new SVInterval(readMetadata.getContigID(read.getContig()), start, end), newReadList);
        return intervalList;
    }

    private static <T extends Locatable> SVIntervalTree<List<T>> combineIntervalLists(final SVIntervalTree<List<T>> intervalTree1,
                                                                                      final SVIntervalTree<List<T>> intervalTree2,
                                                                                      final int clusterSize) {
        return mergeIntervalTrees(intervalTree1, intervalTree2, clusterSize);

    }

    @VisibleForTesting
    static <T extends Locatable> SVIntervalTree<List<T>> mergeIntervalTrees(final SVIntervalTree<List<T>> tree1, final SVIntervalTree<List<T>> tree2, final int clusterSize) {

        if (tree1 == null || tree1.size() == 0) return tree2;
        if (tree2 == null || tree2.size() == 0) return tree1;

        final SVIntervalTree<List<T>> mergedTree = new SVIntervalTree<>();

        PriorityQueue<SVIntervalTree.Entry<List<T>>> nodes = new PriorityQueue<>(tree1.size() + tree2.size(),
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
        List<T> values = new ArrayList<T>();
        while (nodes.size() > 0) {
            final SVIntervalTree.Entry<List<T>> next = nodes.poll();
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

    static <T extends Locatable>  String intervalTreeToBedRecord(final String barcode, final SVIntervalTree.Entry<List<T>> node, final ReadMetadata readMetadata) {
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
        final List<T> reads = node.getValue();
        reads.sort((o1, o2) -> new Integer(o1.getStart()).compareTo(o2.getStart()));
        builder.append("\t");
        builder.append(reads.stream().map(r -> String.valueOf(r.getEnd() - r.getStart() + 1)).collect(Collectors.joining(",")));
        builder.append("\t");
        builder.append(reads.stream().map(r -> String.valueOf(r.getStart() - node.getInterval().getStart())).collect(Collectors.joining(",")));
        return builder.toString();
    }

    static GATKRead intervalTreeToGATKRead(final String barcode, final SAMFileHeader samHeader, final SVIntervalTree.Entry<List<GATKRead>> node, final ReadMetadata readMetadata) {

        final List<GATKRead> reads = node.getValue();
        reads.sort((o1, o2) -> new Integer(o1.getUnclippedStart()).compareTo(o2.getUnclippedStart()));

        final SVInterval interval = node.getInterval();
        int minUnclippedStart = interval.getStart();
        int currentEnd = 0;
        final ByteArrayOutputStream seqOutputStream = new ByteArrayOutputStream();
        final ByteArrayOutputStream qualOutputStream = new ByteArrayOutputStream();
        final Cigar uberCigar = new Cigar();

        int haplotype = -1;
        final Set<String> tenxMolecularIDs = new HashSet<>();
        final Set<String> tenxPhaseSets = new HashSet<>();

        for (final GATKRead read : reads) {
            if (read.hasAttribute("HP")) {
                final Integer readHP = read.getAttributeAsInteger("HP");
                if (haplotype == -1) {
                    haplotype = readHP;
                } else {
                   if (readHP != haplotype) {
                       haplotype = 0;
                   }
                }
            }

            if (read.hasAttribute("MI")) {
                tenxMolecularIDs.add(read.getAttributeAsString("MI"));
            }

            if (read.hasAttribute("PS")) {
                tenxPhaseSets.add(read.getAttributeAsString("PS"));
            }

            final int readStart = read.getUnclippedStart();
            if (readStart < minUnclippedStart) {
                minUnclippedStart = readStart;
            }
            final int readEnd = read.getUnclippedEnd();

            final int chopAmount = Math.max(0, currentEnd - readStart + 1);

            if (chopAmount >= read.getLength()) {
                continue;
            }

            if (chopAmount == 0) {
                if (currentEnd != 0) {
                    final int gap = readStart - currentEnd - 1;
                    if (gap > 0) {
                        uberCigar.add(new CigarElement(gap, CigarOperator.N));
                    }
                }

                seqOutputStream.write(read.getBases(), 0, read.getBases().length);
                qualOutputStream.write(read.getBaseQualities(), 0, read.getBaseQualities().length);

                for (final CigarElement cigarElement : read.getCigarElements()) {
                    if (cigarElement.getOperator() == CigarOperator.H) {
                        continue;
                    }
                    uberCigar.add(translateSoftClip(cigarElement));
                }
            } else {
                int refBasesConsumed = 0;
                int readBasesConsumed = 0;
                for (final CigarElement cigarElement : read.getCigarElements()) {
                    final CigarElement translatedCigarElement = translateSoftClip(cigarElement);
                    if (translatedCigarElement.getOperator() == CigarOperator.H) {
                        continue;
                    }
                    if (refBasesConsumed >= chopAmount) {
                        // no need to chop further, just add cigar operators and bases
                        uberCigar.add(translatedCigarElement);
                        if (translatedCigarElement.getOperator().consumesReadBases()) {
                            seqOutputStream.write(read.getBases(), readBasesConsumed, translatedCigarElement.getLength());
                            qualOutputStream.write(read.getBaseQualities(), readBasesConsumed, translatedCigarElement.getLength());
                            readBasesConsumed = readBasesConsumed + translatedCigarElement.getLength();
                        }
                    } else if (translatedCigarElement.getOperator().consumesReferenceBases() && refBasesConsumed + translatedCigarElement.getLength() > chopAmount) {
                        // need to chop cigar element
                        final int newCigarLength = translatedCigarElement.getLength() - (chopAmount - refBasesConsumed);
                        uberCigar.add(new CigarElement(newCigarLength, translatedCigarElement.getOperator()));
                        if (translatedCigarElement.getOperator().consumesReadBases()) {
                            seqOutputStream.write(read.getBases(), readBasesConsumed + translatedCigarElement.getLength() - newCigarLength, newCigarLength);
                            qualOutputStream.write(read.getBaseQualities(), readBasesConsumed + translatedCigarElement.getLength() - newCigarLength, newCigarLength);
                            readBasesConsumed = readBasesConsumed + translatedCigarElement.getLength();
                        }
                    } else {
                        // skipping cigar element, just add to read bases consumed total
                        if (cigarElement.getOperator().consumesReadBases()) {
                            readBasesConsumed = readBasesConsumed + translatedCigarElement.getLength();
                        }
                    }

                    if (translatedCigarElement.getOperator().consumesReferenceBases()) {
                        refBasesConsumed = refBasesConsumed + translatedCigarElement.getLength();
                    }
                }
            }

            if (currentEnd < readEnd) {
                currentEnd = readEnd;
            }

        }

        final SAMRecord samRecord = new SAMRecord(samHeader);
        samRecord.setReadName(barcode);
        samRecord.setFlags(0x1);
        samRecord.setReferenceName(lookupContigName(interval.getContig(), readMetadata));
        samRecord.setAlignmentStart(minUnclippedStart);
        samRecord.setMappingQuality(60);
        samRecord.setCigar(uberCigar);
        samRecord.setReadBases(seqOutputStream.toByteArray());
        samRecord.setBaseQualities(qualOutputStream.toByteArray());

        samRecord.setAttribute("BX", barcode);

        if (haplotype >= 0) {
            samRecord.setAttribute("HP", haplotype);
        }

        if (! tenxMolecularIDs.isEmpty()) {
            samRecord.setAttribute("MI",  tenxMolecularIDs.stream().collect(Collectors.joining(",")));
        }

        if (! tenxPhaseSets.isEmpty()) {
            samRecord.setAttribute("PS",  tenxPhaseSets.stream().collect(Collectors.joining(",")));
        }

        return new SAMRecordToGATKReadAdapter(samRecord);

    }

    private static String lookupContigName(final int contig, final ReadMetadata readMetadata) {
        for (Map.Entry<String, Integer> entry : readMetadata.getContigNameMap().entrySet()) {
            if (entry.getValue().equals(contig)) return entry.getKey();
        }
        throw new GATKException("Invalid contig ID: "+ contig);
    }

    private static CigarElement translateSoftClip(final CigarElement cigarElement) {
        if (cigarElement.getOperator() == CigarOperator.S) {
            return new CigarElement(cigarElement.getLength(), CigarOperator.M);
        } else {
            return cigarElement;
        }
    }

    static class ReadInfo implements Locatable {
        ReadInfo(final String contig, final int start, final int end) {
            this.contig = contig;
            this.start = start;
            this.end = end;
        }

        String contig;
        int start;
        int end;

        @Override
        public String getContig() {
            return contig;
        }

        @Override
        public int getStart() {
            return start;
        }

        @Override
        public int getEnd() {
            return end;
        }
    }

}

