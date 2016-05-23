package org.broadinstitute.hellbender.tools.spark.linkedreads;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.*;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTree;
import org.apache.commons.io.output.ByteArrayOutputStream;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import scala.Tuple2;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

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

    @Argument(doc = "write SAM",
            shortName = "writeSAM", fullName = "writeSAM",
            optional = true)
    public boolean writeSAM = true;

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

        final JavaPairRDD<String, Map<String, IntervalTree<List<GATKRead>>>> barcodeIntervals =
                reads.mapToPair(read -> new Tuple2<>(read.getAttributeAsString("BX"), read))
                .aggregateByKey(
                        new HashMap<>(),
                        (aggregator, read) -> addReadToIntervals(aggregator, read, clusterSize),
                        (intervalTree1, intervalTree2) -> combineIntervalLists(intervalTree1, intervalTree2, clusterSize)
                );
        if (writeSAM) {
            final JavaRDD<GATKRead> intervalsByBarcode;
            intervalsByBarcode = barcodeIntervals.flatMap(x -> {
                final String barcode = x._1;
                final Map<String, IntervalTree<List<GATKRead>>> contigIntervalTreeMap = x._2;
                final List<GATKRead> results = new ArrayList<>();
                for (final String contig : contigIntervalTreeMap.keySet()) {
                    for (final IntervalTree.Node<List<GATKRead>> next : contigIntervalTreeMap.get(contig)) {
                        results.add(intervalTreeToGATKRead(barcode, contig, headerForReads, next));
                    }
                }
                return results;
            });
            writeReads(ctx, out, intervalsByBarcode);

        } else {
            final JavaRDD<String> intervalsByBarcode;
            intervalsByBarcode = barcodeIntervals.flatMap(x -> {
                final String barcode = x._1;
                final Map<String, IntervalTree<List<GATKRead>>> contigIntervalTreeMap = x._2;
                final List<String> results = new ArrayList<>();
                for (final String contig : contigIntervalTreeMap.keySet()) {
                    for (final IntervalTree.Node<List<GATKRead>> next : contigIntervalTreeMap.get(contig)) {
                        results.add(intervalTreeToBedRecord(barcode, contig, next));
                    }
                }
                return results;
            });
            intervalsByBarcode.saveAsTextFile(out);
        }

    }


    @VisibleForTesting
    static Map<String, IntervalTree<List<GATKRead>>> addReadToIntervals(final Map<String, IntervalTree<List<GATKRead>>> intervalList, final GATKRead read, final int clusterSize) {
        final Interval sloppedReadInterval = new Interval(read.getContig(), read.getStart() - clusterSize, read.getStart()+clusterSize);
        if (! intervalList.containsKey(read.getContig())) {
            intervalList.put(read.getContig(), new IntervalTree<>());
        }
        final IntervalTree<List<GATKRead>> contigIntervals = intervalList.get(sloppedReadInterval.getContig());

        final Iterator<IntervalTree.Node<List<GATKRead>>> iterator = contigIntervals.overlappers(sloppedReadInterval.getStart(), sloppedReadInterval.getEnd());
        int start = read.getStart();
        int end = read.getEnd();
        final List<GATKRead> value = new ArrayList<>();
        value.add(read);
        if (iterator.hasNext()) {
            final IntervalTree.Node<List<GATKRead>> next = iterator.next();
            final int currentStart = next.getStart();
            final int currentEnd = next.getStart();
            final List<GATKRead> currentValue = next.getValue();
            start = Math.min(currentStart, read.getStart());
            end = Math.max(currentEnd, read.getEnd());
            value.addAll(currentValue);
            iterator.remove();
        }
        while (iterator.hasNext()) {
            final IntervalTree.Node<List<GATKRead>> next = iterator.next();
            final int currentEnd = next.getStart();
            final List<GATKRead> currentValue = next.getValue();
            end = currentEnd;
            value.addAll(currentValue);
            iterator.remove();
        }
        contigIntervals.put(start, end, value);
        return intervalList;
    }

    static Map<String, IntervalTree<List<GATKRead>>> combineIntervalLists(final Map<String, IntervalTree<List<GATKRead>>> intervalList1,
                                                                   final Map<String, IntervalTree<List<GATKRead>>> intervalList2,
                                                                   final int clusterSize) {
        final Map<String, IntervalTree<List<GATKRead>>> combinedList = new HashMap<>();
        Stream.concat(intervalList1.keySet().stream(), intervalList2.keySet().stream()).
                forEach(contigName -> combinedList.put(contigName, mergeIntervalTrees(intervalList1.get(contigName), intervalList2.get(contigName), clusterSize)));
        return combinedList;

    }

    @VisibleForTesting
    static IntervalTree<List<GATKRead>> mergeIntervalTrees(final IntervalTree<List<GATKRead>> tree1, final IntervalTree<List<GATKRead>> tree2, final int clusterSize) {

        if (tree1 == null || tree1.size() == 0) return tree2;
        if (tree2 == null || tree2.size() == 0) return tree1;

        tree1.iterator().forEachRemaining(node -> tree2.put(node.getStart(), node.getEnd(), node.getValue()));

        final IntervalTree<List<GATKRead>> mergedTree = new IntervalTree<>();

        while (tree2.size() > 0) {
            final IntervalTree.Node<List<GATKRead>> current = tree2.min();
            final int currentStart = current.getStart();
            int currentEnd = current.getEnd();
            final List<GATKRead> currentValue = current.getValue();

            final Iterator<IntervalTree.Node<List<GATKRead>>> overlappers = tree2.overlappers(current.getStart() - clusterSize, current.getEnd() + clusterSize);
            while (overlappers.hasNext()) {
                final IntervalTree.Node<List<GATKRead>> overlapper = overlappers.next();
                if (overlapper == current) {
                    continue;
                }
                currentEnd = overlapper.getEnd();
                currentValue.addAll(overlapper.getValue());
                overlappers.remove();
            }
            mergedTree.put(currentStart, currentEnd, currentValue);
            tree2.remove(current.getStart(), current.getEnd());
        }
        return mergedTree;
    }

    public static String intervalTreeToBedRecord(final String barcode, final String contig, final IntervalTree.Node<List<GATKRead>> node) {
        final StringBuilder builder = new StringBuilder();
        builder.append(contig);
        builder.append("\t");
        builder.append(node.getStart());
        builder.append("\t");
        builder.append(node.getEnd());
        builder.append("\t");
        builder.append(barcode);
        builder.append("\t");
        builder.append(node.getValue().size());
        builder.append("\t");
        builder.append("+");
        builder.append("\t");
        builder.append(node.getStart());
        builder.append("\t");
        builder.append(node.getEnd());
        builder.append("\t");
        builder.append("0,0,255");
        builder.append("\t");
        builder.append(node.getValue().size());
        final List<GATKRead> reads = node.getValue();
        reads.sort((o1, o2) -> new Integer(o1.getStart()).compareTo(o2.getStart()));
        builder.append("\t");
        builder.append(reads.stream().map(r -> String.valueOf(r.getEnd() - r.getStart() + 1)).collect(Collectors.joining(",")));
        builder.append("\t");
        builder.append(reads.stream().map(r -> String.valueOf(r.getStart() - node.getStart())).collect(Collectors.joining(",")));
        return builder.toString();
    }

    public static GATKRead intervalTreeToGATKRead(final String barcode, final String contig, final SAMFileHeader samHeader, final IntervalTree.Node<List<GATKRead>> node) {

        final List<GATKRead> reads = node.getValue();
        reads.sort((o1, o2) -> new Integer(o1.getUnclippedStart()).compareTo(o2.getUnclippedStart()));

        int minUnclippedStart = node.getStart();
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
        samRecord.setReferenceName(contig);
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

    private static CigarElement translateSoftClip(final CigarElement cigarElement) {
        if (cigarElement.getOperator() == CigarOperator.S) {
            return new CigarElement(cigarElement.getLength(), CigarOperator.M);
        } else {
            return cigarElement;
        }
    }

}

