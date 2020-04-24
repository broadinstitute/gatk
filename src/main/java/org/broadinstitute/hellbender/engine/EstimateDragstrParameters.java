package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.IntervalUtil;
import htsjdk.samtools.util.Locatable;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntList;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.ExampleProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.IntervalPileup;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.ZipUtils;
import org.broadinstitute.hellbender.utils.pairhmm.DragstrCasesSamplerArgumentCollection;
import org.broadinstitute.hellbender.utils.pairhmm.DragstrLocus;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import picard.util.IntervalListTools;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

@CommandLineProgramProperties(
        summary = "summary",
        oneLineSummary = "summary",
        programGroup = ExampleProgramGroup.class
)
public class EstimateDragstrParameters extends GATKTool {

    @Argument(optional = true)
    private int padding = 0;

    @ArgumentCollection
    private DragstrCasesSamplerArgumentCollection config = new DragstrCasesSamplerArgumentCollection();

    @Argument(shortName="sites", fullName="sampling-loci-path", doc="location of the zip that contains the sampling sites for the reference")
    private String lociPath;


    private SAMSequenceDictionary dictionary;
    private long noReads;
    private long lowMinMQ;
    private long lowDepth;

    private DragstrLocus.BinaryTableIndex lociIndex;
    private File lociFile;
    private File lociDir;

    private final Comparator<Locatable> SUBJECT_COMPARATOR = (s1, s2) -> {
        int cmp;
        if ((cmp = Integer.compare(s1.getEnd(), s2.getEnd())) != 0) {
            return cmp;
        } else {
            return Integer.compare(s1.getStart(), s2.getStart());
        }
    };

    private static final Comparator<GATKRead> READ_COMPARATOR = (r1, r2) -> {
        int cmp;
        if (!r1.isUnmapped() && !r2.isUnmapped()) {
            if ((cmp = Integer.compare(r1.getEnd(), r2.getEnd())) != 0) {
                return cmp;
            } else if ((cmp = Integer.compare(r1.getStart(), r2.getStart())) != 0) {
                return cmp;
            } else if ((cmp = Integer.compare(r1.getFlags(), r2.getFlags())) != 0) {
                return cmp;
            } else {
                return r1.getName().compareTo(r2.getName());
            }
        } else {
            if ((cmp = Integer.compare(r1.getAssignedStart(), r2.getAssignedStart())) != 0) {
                return cmp;
            } else if ((cmp = Integer.compare(r1.getFlags(), r2.getFlags())) != 0) {
                return cmp;
            } else if ((cmp = Integer.compare(r1.getFlags(), r2.getFlags())) != 0) {
                return cmp;
            } else {
                return r1.getName().compareTo(r2.getName());
            }
        }
    };

    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    public boolean requiresReads() { return true; }

    @Override
    public void traverse() {
        dictionary = getBestAvailableSequenceDictionary();
        stageLociData();

        final StratifiedResults results = new StratifiedResults(config.maxPeriod, config.maxRepeats);

        if (!hasUserSuppliedIntervals()) {
            dictionary.getSequences().stream()
                    .map(s -> new SimpleInterval(s.getSequenceName(), s.getStart(), s.getEnd()))
                    .map(this::traverse)
                    .forEach(results::addAll);
        } else {
            final List<SimpleInterval> sorted = userIntervals.stream()
                    .sorted(Comparator.comparing(SimpleInterval::getContig).thenComparing(SimpleInterval::getStart).thenComparing(SimpleInterval::getEnd))
                    .collect(Collectors.toList());
            if (sorted.isEmpty()) {
                logger.warn("No intervals to analyze");
                return;
            } else {
                final Deque<SimpleInterval> combined = new ArrayDeque<>(sorted.size());
                combined.add(sorted.get(0));
                for (final SimpleInterval si : sorted) {
                    if (si.overlapsWithMargin(combined.getLast(), 0)) {
                        combined.add(si.mergeWithContiguous(combined.removeLast()));
                    } else {
                        combined.add(si);
                    }
                }
                combined.stream()
                        .map(s -> new SimpleInterval(s.getContig(), s.getStart(), s.getEnd()))
                        .map(this::traverse)
                        .forEach(results::addAll);
            }

        }
        if (logger.isDebugEnabled()) {
            logger.debug("Number of usable cases per period and repeat length:\n");
            final int[] columnWidths = IntStream.range(1, config.maxPeriod + 1).map(period -> {
                final int max = IntStream.range(1, config.maxRepeats).map(repeat -> results.get(period,repeat).size())
                        .max().orElse(0);
                return (int) Math.ceil(Math.log10(max)) + 1; }).toArray();
            logger.debug("      " + IntStream.range(0, config.maxPeriod).mapToObj(i -> String.format("-%" + columnWidths[i] + "s", (i + 1))).collect(Collectors.joining()));
            for (int i = 1; i <= config.maxRepeats; i++) {
                final int repeat = i;
                logger.debug(String.format("%-4s", repeat) + "  " + IntStream.range(1, config.maxPeriod + 1)
                        .mapToObj(period -> String.format("%-" + columnWidths[period - 1] + "s",
                                results.get(period, repeat).size())).collect(Collectors.joining("")));
            }
        }
        unstageLociData();
    }

    private void unstageLociData() {
        Utils.deleteFileTree(lociDir);
    }

    private void stageLociData() {
        try {
            lociDir = File.createTempFile("loci", ".dir");
        } catch (final IOException e) {
            throw new UserException.CouldNotReadInputFile(lociPath, "could not unzip it", e);
        }
        lociDir.delete();
        lociDir.mkdir();
        ZipUtils.unzip(lociPath, lociDir);
        lociFile = new File(lociDir, "all.bin");
        lociFile.deleteOnExit();
        try {
            lociIndex = DragstrLocus.BinaryTableIndex.load(new File(lociDir, "all.idx").toString());
        } catch (final IOException ex) {
            throw new UserException.CouldNotReadInputFile(lociPath, "could not open the index", ex);
        }
    }

    private static class TraversalResult {
        public List<DragstrLocus> loci;
        public IntList depth;
        public IntList nonref;

        public TraversalResult() {
            loci = new ArrayList<>();
            depth = new IntArrayList();
            nonref = new IntArrayList();
        }

        public void add(final DragstrLocus locus, final int total, final int nonref) {
            loci.add(locus);
            depth.add(total);
            this.nonref.add(nonref);
        }

        public int size() {
            return loci.size();
        }
    }

    private static class StratifiedResults {
        public TraversalResult[][] perPeriodAndRepeat;

        public StratifiedResults(final int maxPeriod, final int maxRepeats) {
            perPeriodAndRepeat = new TraversalResult[maxPeriod][maxRepeats];
            for (int i = 0; i < maxPeriod; i++) {
                for (int j = 0; j < maxRepeats; j++) {
                    perPeriodAndRepeat[i][j] = new TraversalResult();
                }
            }
        }

        public TraversalResult get(final int period, final int repeat) {
            return perPeriodAndRepeat[period - 1][repeat - 1];
        }

        final void addAll(final TraversalResult traversalResult) {
            final int size = traversalResult.loci.size();
            for (int i = 0; i < size; i++) {
                final DragstrLocus locus = traversalResult.loci.get(i);
                final int depth = traversalResult.depth.getInt(i);
                final int nonref = traversalResult.nonref.getInt(i);
                final TraversalResult[] periodResults = perPeriodAndRepeat[Math.min(perPeriodAndRepeat.length - 1, locus.getPeriod() - 1)];
                final TraversalResult results = periodResults[Math.min(periodResults.length - 1, locus.getRepeats() - 1)];
                results.add(locus, depth, nonref);
            }
        }
    }

    private static enum LocusResultCode {
        OK, LOW_DEPTH, LOW_MQ, LOW_DEPTH2
    }

    private static class LocusResult {
        public final DragstrLocus locus;
        public final LocusResultCode code;
        public final boolean qualifies;
        public final int depth;
        public final int nonrefDepth;

        public LocusResult(DragstrLocus locus, LocusResultCode code, int total, int nonRef) {
            this.locus = locus;
            this.code = code;
            this.depth = total;
            this.nonrefDepth = nonRef;
            this.qualifies = this.code == LocusResultCode.OK;
        }

        public static LocusResult qualifying(final DragstrLocus locus, final int total, final int nonRef) {
            return new LocusResult(locus, LocusResultCode.OK, total, nonRef);
        }

        public static LocusResult noReads(final DragstrLocus locus) {
            return new LocusResult(locus, LocusResultCode.LOW_DEPTH, 0, 0);
        }

        public static LocusResult lowDepth(final DragstrLocus locus, final int total) {
            return new LocusResult(locus, LocusResultCode.LOW_DEPTH, total, -1);
        }

        public static LocusResult lowDepth2(final DragstrLocus locus, final int total, final int nonref) {
            return new LocusResult(locus, LocusResultCode.LOW_DEPTH2, total, nonref);
        }

        public static LocusResult lowMq(final DragstrLocus locus, final int total) {
            return new LocusResult(locus, LocusResultCode.LOW_MQ, total, -1);
        }
    }

    private TraversalResult traverse(final SimpleInterval interval) {
        final Stream<GATKRead> readStream = getTransformedReadStream(interval, makeReadFilter());
        final Stream<DragstrLocus> lociStream = getSubjectStream(interval);
        final Iterator<GATKRead> readIterator = readStream.iterator();
        final Iterator<DragstrLocus> subjectIterator = lociStream.iterator();
        final IntervalBuffer<GATKRead> readBuffer = new IntervalBuffer<>();
        final TraversalResult result = new TraversalResult();
        GATKRead lastRead = readIterator.hasNext() ? readIterator.next() : null;
        DragstrLocus lastSubject = subjectIterator.hasNext() ? subjectIterator.next() : null;
        long[] caseCounts = new long[LocusResultCode.values().length];
        long okCaseWithNonref = 0;
        long refReads = 0;
        long nonRefReads = 0;
        while (lastSubject != null) {
            readBuffer.removeUpstreamFrom((int) lastSubject.getStart());
            while (lastRead != null && lastRead.getAssignedStart() <= lastSubject.getEnd()) {
                readBuffer.add(lastRead.getAssignedStart(), readEnd(lastRead), lastRead);
                lastRead = readIterator.hasNext() ? readIterator.next() : null;
            }
            final List<GATKRead> reads = readBuffer.overlapping((int) lastSubject.getStart(), (int) lastSubject.getEnd());
            final LocusResult locusResult = evaluateLocus(lastSubject, reads);
            progressMeter.update(new SimpleInterval(lastSubject.getStartInterval(dictionary, 0)));
            if (locusResult.qualifies) {
                if (locusResult.nonrefDepth > 0) okCaseWithNonref++;
                refReads += locusResult.depth - locusResult.nonrefDepth;
                nonRefReads += locusResult.nonrefDepth;
                result.add(lastSubject, locusResult.depth, locusResult.nonrefDepth);
            }
            caseCounts[locusResult.code.ordinal()]++;
            lastSubject = subjectIterator.hasNext() ? subjectIterator.next() : null;
        }
        if (logger.isDebugEnabled()) {
            logger.debug("Summary data collection " + interval + ":\n"
                    + "    " + Arrays.stream(LocusResultCode.values()).map(c -> "" + c + " = " + caseCounts[c.ordinal()]).collect(Collectors.joining("\n    "))
                    + "\n"
                    + "    ref/non-ref reads = " + refReads + "/" + nonRefReads + "\n"
                    + "    ref-only/non-ref-containing reads cases = " + (caseCounts[0] - okCaseWithNonref) + "/" + okCaseWithNonref + "\n");
        }
        return result;
    }

    private int readEnd(final GATKRead read) {
        return read.isUnmapped() ? read.getAssignedStart() : read.getEnd();
    }
    private Stream<DragstrLocus> getSubjectStream(SimpleInterval interval) {
        try {
            return DragstrLocus.binaryReader(new File(lociDir, "all.bin").toString(), lociIndex,
                    dictionary.getSequenceIndex(interval.getContig()), interval.getStart(), interval.getEnd()).stream();
        } catch (final IOException ex) {
            throw new UserException.CouldNotReadInputFile(lociFile.toString(), ex);
        }
    }

    private LocusResult evaluateLocus(DragstrLocus locus, List<GATKRead> rawReads) {
        final List<GATKRead> reads = rawReads.stream()
                .filter(m -> (m.getFlags() & 0x0f04) == 0)
                .filter(m -> m.getMappingQuality() != 0)
                .filter(m -> !m.isUnmapped())
                .collect(Collectors.toList());
        if (reads.isEmpty()) {
            noReads++;
            lowDepth++;
            return LocusResult.noReads(locus);
        } else if (reads.stream().mapToInt(GATKRead::getMappingQuality).min().getAsInt() < config.samplingMinMQ) {
            lowMinMQ++;
            return LocusResult.lowMq(locus, reads.size());
        } else if (reads.size() < config.minDepth) {
            return LocusResult.lowDepth(locus, reads.size());
        } else {
            final int period = locus.getPeriod();
            final SimpleInterval repeatInterval = locus.getRepeatInterval(dictionary, config.pileupPadding, config.pileupPadding + period);
            byte[] referenceSequence = directlyAccessEngineReferenceDataSource().queryAndPrefetch(repeatInterval).getBases();
            final ReferenceBases referenceBases = new ReferenceBases(referenceSequence, repeatInterval);
            final IntervalPileup pileup = IntervalPileup.of(reads, referenceBases);
            final int startColumn = (int) locus.getStart() - repeatInterval.getStart();
            int j, k, l;
            for (j = startColumn, k = j + locus.getRepeats() * period, l = 0; l < period && k < referenceSequence.length; l++, j++, k++) {
                if (referenceSequence[j] != referenceSequence[k]) {
                    break;
                }
            }
            final int endOfRepeatWithExtraBases = l + startColumn + locus.getRepeats() * period - 1;
            final int firstRelevantColumn = 0;
            final int lastRelevantColumn = Math.min(endOfRepeatWithExtraBases + config.pileupPadding, pileup.width() - 1);
            int total = 0;
            int nonRefTotal = 0;
            int readsWithNoBases = 0;
            int readsWithTooManyBadBQs = 0;
            int readsWithNoQuals = 0;
            row_for:
            for (int row = 0; row < pileup.height(); row++) {
                if (!pileup.reads().get(row).hasBaseQualities()) {
                    readsWithNoBases++;
                    continue;
                }
                int disqualifyingBaseCalls = 0;
                int indelLength = 0;
                final IntervalPileup.Element element = pileup.element(row);
                for (int column = firstRelevantColumn; column <= lastRelevantColumn; column++) {
                    final byte qual = pileup.qualAt(row, column);
                    final byte base = pileup.baseAt(row, column);
                    if (base == IntervalPileup.NO_BASE) {
                        //if (++disqualifyingBaseCalls > config.baseQualExceptionsAllowed) {
                        //readsWithTooManyBadBQs++;
                        readsWithNoBases++;
                        continue row_for;
                        //}
                    } else if (base == IntervalPileup.GAP) {
                        if (column >= startColumn && column <= endOfRepeatWithExtraBases) {
                            indelLength--;
                        }
                    } else if (qual != IntervalPileup.NO_BQ) {
                        if (qual < config.baseQualThreshold
                                && ++disqualifyingBaseCalls > config.baseQualExceptionsAllowed) {
                            readsWithTooManyBadBQs++;
                            continue row_for;
                        }
                    } else { //if (qual == IntervalPileup.NO_BQ) {
                        readsWithNoQuals++;
                        continue row_for;
                    }
                }
                for (final IntervalPileup.Insert insert : element.inserts(startColumn - 1, endOfRepeatWithExtraBases)) {
                    indelLength += insert.length();
                }
                total++;
                if (indelLength != 0) {
                    nonRefTotal++;
                }
            }
            if (total < config.minDepth) {
                lowDepth++;
                return LocusResult.lowDepth2(locus, total, nonRefTotal);
            } else {
                return LocusResult.qualifying(locus, total, nonRefTotal);
            }
        }
    }


    private static class IntervalBuffer<E> extends IntervalTree<List<E>> {

        private static <E> List<E> mergeLists(final List<E> l1, final List<E> l2) {
            if (l1 instanceof ArrayList) {
                l1.addAll(l2);
                return l1;
            } else if (l2 instanceof ArrayList) {
                l2.addAll(l1);
                return l2;
            } else {
                final List<E> l3 = new ArrayList<>(l1.size() + l2.size());
                l3.addAll(l1);
                l3.addAll(l2);
                return l3;
            }
        }

        public void add(final int start, final int end, final E elem) {
            merge(start, end, Collections.singletonList(elem), IntervalBuffer::mergeLists);
        }

        public List<E> removeUpstreamFrom(final int start) {
            final Iterator<Node<List<E>>> it = iterator();
            List<E> result = null;
            while (it.hasNext()) {
                final Node<List<E>> node = it.next();
                if (node.getStart() >= start) {
                    break;
                } else if (node.getEnd() < start) {
                    final List<E> values = node.getValue();
                    if (result == null) {
                        result = new ArrayList<>(values.size() + 10);
                    }
                    result.addAll(values);
                    it.remove();
                }
            }
            return result == null ? Collections.emptyList() : result;
        }

        public List<E> overlapping(final int start, final int end) {
            Iterator<Node<List<E>>> it = this.overlappers(start, end);
            if (!it.hasNext()) {
                return Collections.emptyList();
            } else {
                final List<E> result = new ArrayList<>(100);
                do {
                    final Node<List<E>> node = it.next();
                    result.addAll(node.getValue());
                } while (it.hasNext());
                return result;
            }
        }
    }


}
