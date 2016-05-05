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
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.read.GATKRead;
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

    @Override
    public boolean requiresReads() { return true; }

    @Override
    public ReadFilter makeReadFilter() {
        return super.makeReadFilter().and(new LinkedReadAnalysisFilter());
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
        final JavaRDD<String> intervalsByBarcode;
        if (writeSAM) {
            intervalsByBarcode = barcodeIntervals.flatMap(x -> barcodeLocationsToSam(x, headerForReads));
        } else {
            intervalsByBarcode = barcodeIntervals.flatMap(x -> barcodeLocationsToBed(x));
        }

        intervalsByBarcode.saveAsTextFile(out);
    }

    static List<String> barcodeLocationsToBed(final Tuple2<String, Map<String, IntervalTree<List<GATKRead>>>> barcodeLocations) {
        return barcodeLocationsToString(barcodeLocations, new BedRecordStringifier());
    }

    public static List<String> barcodeLocationsToSam(final Tuple2<String, Map<String, IntervalTree<List<GATKRead>>>> barcodeLocations, final SAMFileHeader samHeader) {
        return barcodeLocationsToString(barcodeLocations, new SamRecordStringifier(samHeader));
    }

    static List<String> barcodeLocationsToString(final Tuple2<String, Map<String, IntervalTree<List<GATKRead>>>> barcodeLocations, final RecordStringifier intervalTreeStringifier) {
        final List<String> outputLines = new ArrayList<>();
        final String barcode = barcodeLocations._1;
        final Map<String, IntervalTree<List<GATKRead>>> stringIntervalTreeMap = barcodeLocations._2;
        for (final String contig : stringIntervalTreeMap.keySet()) {
            final IntervalTree<List<GATKRead>> contigIntervalTree = stringIntervalTreeMap.get(contig);
            for (final IntervalTree.Node<List<GATKRead>> node : contigIntervalTree) {
                final String bedRecord = intervalTreeStringifier.stringify(barcode, contig, node);
                outputLines.add(bedRecord);
            }
        }

        return outputLines;
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

    public interface RecordStringifier {
        String stringify(final String barcode, final String contig, final IntervalTree.Node<List<GATKRead>> node);
    }

    public static class BedRecordStringifier implements RecordStringifier {
        @Override
        public String stringify(final String barcode, final String contig, final IntervalTree.Node<List<GATKRead>> node) {
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
    }

    public static class SamRecordStringifier implements RecordStringifier {
        private final SAMFileHeader samHeader;

        public SamRecordStringifier(final SAMFileHeader samHeader) {
            this.samHeader = samHeader;
        }

        @Override
        public String stringify(final String barcode, final String contig, final IntervalTree.Node<List<GATKRead>> node) {

            final List<GATKRead> reads = node.getValue();
            reads.sort((o1, o2) -> new Integer(o1.getStart()).compareTo(o2.getStart()));

            int minUnclippedStart = node.getStart();
            int currentEnd = 0;
            final ByteArrayOutputStream seqOutputStream = new ByteArrayOutputStream();
            final ByteArrayOutputStream qualOutputStream = new ByteArrayOutputStream();
            final Cigar uberCigar = new Cigar();

            for (GATKRead read : reads) {
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

                    for (CigarElement cigarElement : read.getCigarElements()) {
                        if (cigarElement.getOperator() == CigarOperator.H) {
                            continue;
                        }
                        uberCigar.add(translateSoftClip(cigarElement));
                    }
                } else {
                    int refBasesConsumed = 0;
                    int readBasesConsumed = 0;
                    for (CigarElement cigarElement : read.getCigarElements()) {
                        final CigarElement translatedCigarElement = translateSoftClip(cigarElement);
                        if (translatedCigarElement.getOperator() == CigarOperator.H) {
                            continue;
                        }
                        if (refBasesConsumed >= chopAmount) {
                            // no need to chop further, just add cigar operators and bases
                            uberCigar.add(translatedCigarElement);
                            if (translatedCigarElement.getOperator().consumesReadBases()) {
                                try {
                                    seqOutputStream.write(read.getBases(), readBasesConsumed, translatedCigarElement.getLength());
                                    qualOutputStream.write(read.getBaseQualities(), readBasesConsumed, translatedCigarElement.getLength());
                                } catch (IndexOutOfBoundsException e) {
                                    throw new GATKException("read = " + read + "\n" +
                                            "cigar = " + read.getCigar() + "\n" +
                                            "cigarElement  = " + cigarElement + "\n" +
                                            "translatedCigarElement = " + translatedCigarElement + "\n" +
                                            "uberCigar = " + uberCigar + "\n" +
                                            "readBasesConsumed = " + readBasesConsumed + "\n" +
                                            "refBasesConsumed = " + refBasesConsumed + "\n" +
                                            "chopAmount = " + chopAmount + "\n", e);
                                }
                                readBasesConsumed = readBasesConsumed + translatedCigarElement.getLength();
                            }
                        } else if (translatedCigarElement.getOperator().consumesReferenceBases() && refBasesConsumed + translatedCigarElement.getLength() > chopAmount) {
                            // need to chop cigar element
                            final int newCigarLength = translatedCigarElement.getLength() - (chopAmount - refBasesConsumed);
                            uberCigar.add(new CigarElement(newCigarLength, translatedCigarElement.getOperator()));
                            if (translatedCigarElement.getOperator().consumesReadBases()) {
                                seqOutputStream.write(read.getBases(), readBasesConsumed + translatedCigarElement.getLength() - newCigarLength, newCigarLength);
                                qualOutputStream.write(read.getBaseQualities(), readBasesConsumed + translatedCigarElement.getLength() - newCigarLength, newCigarLength);
                                readBasesConsumed = readBasesConsumed + translatedCigarElement.getLength() - chopAmount;
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

            return samRecord.getSAMString();

        }

        private CigarElement translateSoftClip(final CigarElement cigarElement) {
            if (cigarElement.getOperator() == CigarOperator.S) {
                return new CigarElement(cigarElement.getLength(), CigarOperator.M);
            } else {
                return cigarElement;
            }
        }

    }
}
