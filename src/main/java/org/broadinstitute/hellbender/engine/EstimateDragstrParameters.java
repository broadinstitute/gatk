package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.IntervalTree;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntList;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.ExampleProgramGroup;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.transformers.ReadTransformer;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.pairhmm.*;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.atomic.AtomicLong;
import java.util.function.BiConsumer;
import java.util.function.BinaryOperator;
import java.util.function.Function;
import java.util.function.Supplier;
import java.util.stream.Collector;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

@CommandLineProgramProperties(
        summary = "summary",
        oneLineSummary = "summary",
        programGroup = ExampleProgramGroup.class
)
public class EstimateDragstrParameters extends GATKTool {

    @ArgumentCollection
    private DragstrCasesSamplerArgumentCollection config = new DragstrCasesSamplerArgumentCollection();

    @Argument(shortName="sites", fullName="sampling-loci-path", doc="location of the zip that contains the sampling sites for the reference")
    private String lociPath = null;

    @Argument(fullName="parallel", doc="run estimation in parallel")
    private boolean runInParallel;

    @Argument(fullName="threads", shortName="N", minValue = 0.0, doc="suggested number of parallel threads to perform the estimation, "
            + "the default 0 leave it up to the VM to decide. When set to more than 1, this will activate parallel in the absence of --parallel")
    private int threads = 0;

    @Argument(fullName="shard-size", doc="when running in parallel this is the suggested shard size in base pairs. " +
            "The actual shard-size may vary to adapt to small contigs and the requested number of threads")
    private int shardSize = 1_000_000;

    @Argument(fullName="downsample", doc="Targeted maximum number of cases per combination period repeat count, " +
            "the larger the more precise but also the slower estimation.")
    private int downsample = 4096;

    @Argument(fullName="output", shortName = "O", doc = "where to write the parameter output file.")
    private String output;

    private SAMSequenceDictionary dictionary;
    private AbsoluteCoordinates absoluteCoordinates;
    private long noReads;
    private long lowMinMQ;
    private long lowDepth;


    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    public boolean requiresReads() { return true; }

    @Override
    public void traverse() {
        dictionary = getBestAvailableSequenceDictionary();
        absoluteCoordinates = AbsoluteCoordinates.of(dictionary);
        try (final AutoCloseableReference<File> stagedLociData = stageLociData()) {
            final File lociDir = stagedLociData.get();

            final StratifiedDragstrLocusCases allSites;
            final List<SimpleInterval> intervals = !hasUserSuppliedIntervals() ?
                    dictionary.getSequences().stream()
                            .map(s -> new SimpleInterval(s.getSequenceName(), s.getStart(), s.getEnd())).collect(Collectors.toList())
                    : sortAndMergeOverlappingIntervals(userIntervals);

            if (intervals.isEmpty()) {
                logger.warn("No intervals to analyze");
                return;
            }

            // If we aren't running in parallel there is no need to sharding:

            if (runInParallel || threads > 1) {
                allSites = traverseParallel(intervals, stagedLociData.get(), shardSize);
            } else {
                allSites = traverseSerial(intervals, stagedLociData.get());
            }
            if (isThereEnoughCases(allSites)) {
                logCounts(allSites);
                final StratifiedDragstrLocusCases finalSites = downSample(allSites, lociDir);
                logCounts(finalSites);

                logger.info("Estimating parameters used sampled down cases");
                final DragstrParams estimate = estimateParams(finalSites);
                logger.info("Done with estimation, printing output");
                estimate.print(output);
            } else {
                DragstrParams.DEFAULT.print(output);
            }
        }
    }

    /**
     * Check that a minimum number of cases are available in key bins (combo period, repeat).
     */
    private boolean isThereEnoughCases(final StratifiedDragstrLocusCases allSites) {
        for (int i = 1; i <= 9; i++) {
            if (allSites.get(1, i).size() < 200) {
                return false;
            }
        }
        for (int i = 2; i <= 5; i++) {
            if (allSites.get(2, i).size() < 200) {
                return false;
            }
        }
        for (int i = 2; i <= 4; i++) {
            if (allSites.get(3, i).size() < 200) {
                return false;
            }
        }
        if (allSites.get(4, 3).size() < 200) {
            return false;
        }
        for (int i = 4; i <= 8; i++) {
            if (allSites.get(i, 2).size() < 200) {
                return false;
            }
        }
        return true;
    }

    private DragstrParams estimateParams(StratifiedDragstrLocusCases finalSites) {
        final DragstrParametersEstimator estimator = new DragstrParametersEstimator(config);
        if (runInParallel || threads > 1) {
            return estimator.estimateParallel(finalSites, threads);
        } else {
            return estimator.estimateSequencial(finalSites);
        }
    }

    private StratifiedDragstrLocusCases downSample(final StratifiedDragstrLocusCases allSites, final File lociFile) {
        final STRDecimationTable decimationTable = new STRDecimationTable(new File(lociFile, "decimation.txt").getAbsolutePath());
        final List<PeriodRepeatCombo> prCombos = new ArrayList<>(config.maxPeriod * config.maxRepeats);
        for (int i = 1; i <= config.maxPeriod; i++) {
            for (int j = 1; j <= config.maxRepeats; j++) {
                prCombos.add(PeriodRepeatCombo.of(i, j));
            }
        }

        final StratifiedDragstrLocusCases finalSites;
        if (runInParallel) {
            finalSites = Utils.runInParallel(threads, () -> prCombos.stream().parallel()
                .map(combo -> {
                    final DragstrLocusCases all = allSites.perPeriodAndRepeat[combo.period - 1][combo.repeat - 1];
                    final int decimationBit = decimationTable.decimationBit(combo.period, combo.repeat);
                    return downsample(all, decimationBit, downsample);
                }).collect(Stratificator.make(config.maxPeriod, config.maxRepeats)));

        } else {
            finalSites = prCombos.stream().sequential()
                    .map(combo -> {
                        final DragstrLocusCases all = allSites.perPeriodAndRepeat[combo.period - 1][combo.repeat - 1];
                        final int decimationBit = decimationTable.decimationBit(combo.period, combo.repeat);
                        return downsample(all, decimationBit, downsample);
                    }).collect(Stratificator.make(config.maxPeriod, config.maxRepeats));
        }
        return finalSites;
    }

    private static final long[] DECIMATION_MASKS_BY_BIT = new long[Long.SIZE];
    static {
        DECIMATION_MASKS_BY_BIT[0] = 1;
        for (int i = 1, j = 0; i < Long.SIZE; i++, j++) {
            DECIMATION_MASKS_BY_BIT[i] = DECIMATION_MASKS_BY_BIT[j] << 1;
            DECIMATION_MASKS_BY_BIT[j] = ~DECIMATION_MASKS_BY_BIT[j];
        }
        DECIMATION_MASKS_BY_BIT[Long.SIZE -1] = ~DECIMATION_MASKS_BY_BIT[Long.SIZE - 1];
    }

    private DragstrLocusCases downsample(final DragstrLocusCases in, final int minDecimationBit, final int maxCount) {
        final int inSize = in.size();
        if (inSize <= maxCount) {
            return in;
        } else {
            final int[] countByFirstDecimatingBit = new int[Long.SIZE - minDecimationBit];
            for (final DragstrLocus locus : in.loci) {
                long mask = locus.getMask();
                for (int j = minDecimationBit; mask != 0 && j < 64;  j++) {
                    final long newMask = mask & DECIMATION_MASKS_BY_BIT[j];
                    if (newMask != mask) {
                         countByFirstDecimatingBit[j]++;
                         break;
                     } else {
                         mask = newMask;
                     }
                }
            }
            int finalSize = inSize;
            IntList progressiveSizes = new IntArrayList(10);
            long filterMask = 0;
            progressiveSizes.add(finalSize);
            for (int j = minDecimationBit; finalSize > maxCount && j < Long.SIZE; j++) {
                finalSize -= countByFirstDecimatingBit[j];
                filterMask |= ~ DECIMATION_MASKS_BY_BIT[j];
                progressiveSizes.add(finalSize);
            }
            final DragstrLocusCases out = new DragstrLocusCases(finalSize);
            for (int i = 0; i < inSize; i++) {
                final long mask = in.loci.get(i).getMask();
                if ((mask & filterMask) == filterMask) {
                    out.loci.add(in.loci.get(i));
                    out.depth.add(in.depth.get(i));
                    out.nonref.add(in.nonref.get(i));
                }
            }
            progressiveSizes.add(out.size());
            if (logger.isDebugEnabled()) {
                logger.debug("" + in.loci.get(0).getPeriod() + " "  + in.loci.get(0).getRepeats() + " " + Arrays.toString(progressiveSizes.toArray()));
            }
            return out;
        }
    }

    private void logCounts(StratifiedDragstrLocusCases allSites) {
        if (logger.isDebugEnabled()) {
            logger.debug("Number of usable cases per period and repeat length:\n");
            final int[] columnWidths = IntStream.range(1, config.maxPeriod + 1).map(period -> {
                final int max = IntStream.range(1, config.maxRepeats).map(repeat -> allSites.get(period,repeat).size())
                        .max().orElse(0);
                return (int) Math.ceil(Math.log10(max)) + 1; }).toArray();
            logger.debug("      " + IntStream.range(0, config.maxPeriod).mapToObj(i -> String.format("-%" + columnWidths[i] + "s", (i + 1))).collect(Collectors.joining()));
            for (int i = 1; i <= config.maxRepeats; i++) {
                final int repeat = i;
                logger.debug(String.format("%-4s", repeat) + "  " + IntStream.range(1, config.maxPeriod + 1)
                        .mapToObj(period -> String.format("%-" + columnWidths[period - 1] + "s",
                                allSites.get(period, repeat).size())).collect(Collectors.joining("")));
            }
        }
    }

    private StratifiedDragstrLocusCases traverseSerial(List<SimpleInterval> intervals, final File lociDir) {
        return intervals.stream()
                .map(interval -> { final DragstrLocusCases result = traverse(interval, directlyAccessEngineReadsDataSource(), directlyAccessEngineReferenceDataSource(), lociDir);
                                   if (!result.loci.isEmpty()) {
                                       progressMeter.update(result.loci.get(result.loci.size() - 1).getStartInterval(dictionary, 0), result.loci.size());
                                   }
                                   return result;
                })
                .collect(Stratificator.make(config.maxPeriod, config.maxRepeats));
    }

    private StratifiedDragstrLocusCases traverseParallel(final List<SimpleInterval> intervals, final File lociDir, final int shardSize) {
        final List<SimpleInterval> shards = shardIntervals(intervals, shardSize);
        try (final ReadsDataSourcePool readsDataSourcePool = new ReadsDataSourcePool(readArguments.getReadPaths());
             final ReferenceDataSourcePool referenceDataSourcePool = new ReferenceDataSourcePool(referenceArguments.getReferencePath())) {
             final EstimateDragstrParameters tool = this;

            final AtomicLong numberBasesProcessed = new AtomicLong(0);
            return Utils.runInParallel(threads, () ->
                        shards.stream().parallel()
                                .map(shard -> {
                                    try (ReadsDataSourcePool.AutoReturn readsSource = readsDataSourcePool.borrowAutoReturn();
                                         ReferenceDataSourcePool.AutoReturn referenceSource = referenceDataSourcePool.borrowAutoReturn()) {
                                         final DragstrLocusCases result = tool.traverse(shard, readsSource.get(), referenceSource.get(), lociDir);
                                         final int resultSize = result.size();
                                         synchronized (numberBasesProcessed) {
                                             final long processed = numberBasesProcessed.updateAndGet(l -> l + shard.size());
                                             if (resultSize > 0) {
                                                 progressMeter.update(absoluteCoordinates.toSimpleInterval(processed, 1), resultSize);
                                             }
                                         }
                                         return result;
                                    }
                                }).collect(Stratificator.make(config.maxPeriod, config.maxRepeats))
            );
        }
    }

    private List<SimpleInterval> shardIntervals(final List<SimpleInterval> unshared, final int shardSize) {
        final List<SimpleInterval> output = new ArrayList<>(10 + unshared.size() * 10);
        final int shardingSizeThreshold = (int) Math.max(100, Math.round(shardSize * 1.5));
        for (final SimpleInterval in : unshared) {
            if (in.size() < shardingSizeThreshold) {
                output.add(in);
            } else {
                int start = in.getStart();
                final int inEnd = in.getEnd();
                final int stop = in.getEnd() - shardingSizeThreshold + 1;
                while (start < stop) {
                    final int end = start + shardSize - 1;
                    output.add(new SimpleInterval(in.getContig(), start, end));
                    start = end + 1;
                }
                if (start <= inEnd) {
                    output.add(new SimpleInterval(in.getContig(), start, inEnd));
                }
            }
        }
        return output;
    }

    private List<SimpleInterval> sortAndMergeOverlappingIntervals(final List<SimpleInterval> userIntervals) {
        final List<SimpleInterval> sorted = userIntervals.stream()
                .sorted(Comparator.comparing(SimpleInterval::getContig).thenComparing(SimpleInterval::getStart)
                        .thenComparing(SimpleInterval::getEnd))
                .collect(Collectors.toList());
        if (sorted.isEmpty()) {
            return Collections.emptyList();
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
            return new ArrayList<>(combined);
        }
    }

    private AutoCloseableReference<File> stageLociData() {
        final File lociDir;
        try {
            lociDir = File.createTempFile("loci", ".dir");
        } catch (final IOException e) {
            throw new UserException.CouldNotReadInputFile(lociPath, "could not unzip it", e);
        }
        lociDir.delete();
        lociDir.mkdir();
        ZipUtils.unzip(lociPath, lociDir);
        return AutoCloseableReference.of(lociDir, Utils::deleteFileTree);
    }

    public static class DragstrLocusCases {
        public List<DragstrLocus> loci;
        public IntList depth;
        public IntList nonref;

        public DragstrLocusCases() {
            this(100);
        }

        public DragstrLocusCases(int initialCapacity) {
            loci = new ArrayList<>(initialCapacity);
            depth = new IntArrayList(initialCapacity);
            nonref = new IntArrayList(initialCapacity);
        }

        public void add(final DragstrLocus locus, final int total, final int nonref) {
            loci.add(locus);
            depth.add(total);
            this.nonref.add(nonref);
        }

        public void addAll(final DragstrLocusCases other) {
            loci.addAll(other.loci);
            depth.addAll(other.depth);
            nonref.addAll(other.nonref);
        }

        public int size() {
            return loci.size();
        }
    }

    private static class Stratificator implements Collector<DragstrLocusCases, StratifiedDragstrLocusCases, StratifiedDragstrLocusCases> {


        private final int maxPeriod;
        private final int maxRepeats;

        private static Stratificator make(final int maxPeriod, final int maxRepeats) {
            return new Stratificator(maxPeriod, maxRepeats);
        }

        private Stratificator(final int maxPeriod, final int maxRepeats) {
            this.maxPeriod = maxPeriod;
            this.maxRepeats = maxRepeats;
        }

        @Override
        public Supplier<StratifiedDragstrLocusCases> supplier() {
            return () -> new StratifiedDragstrLocusCases(maxPeriod, maxRepeats);
        }

        @Override
        public BiConsumer<StratifiedDragstrLocusCases, DragstrLocusCases> accumulator() {
            return StratifiedDragstrLocusCases::addAll;
        }

        @Override
        public BinaryOperator<StratifiedDragstrLocusCases> combiner() {
            return StratifiedDragstrLocusCases::addAll;
        }

        @Override
        public Function<StratifiedDragstrLocusCases, StratifiedDragstrLocusCases> finisher() {
            return a -> a;
        }

        @Override
        public Set<Characteristics> characteristics() {
            return EnumSet.of(Characteristics.IDENTITY_FINISH, Characteristics.UNORDERED);
        }
    }

    public static class StratifiedDragstrLocusCases {
        public DragstrLocusCases[][] perPeriodAndRepeat;

        public StratifiedDragstrLocusCases(final int maxPeriod, final int maxRepeats) {
            perPeriodAndRepeat = new DragstrLocusCases[maxPeriod][maxRepeats];
            for (int i = 0; i < maxPeriod; i++) {
                for (int j = 0; j < maxRepeats; j++) {
                    perPeriodAndRepeat[i][j] = new DragstrLocusCases();
                }
            }
        }

        public StratifiedDragstrLocusCases addAll(final StratifiedDragstrLocusCases other) {
            Utils.validate(perPeriodAndRepeat.length == other.perPeriodAndRepeat.length, "incompatible dimensions");
            for (int i = 0; i < perPeriodAndRepeat.length; i++) {
                Utils.validate(perPeriodAndRepeat[i].length == other.perPeriodAndRepeat[i].length, "invalid dimensions");
                for (int j = 0; j < perPeriodAndRepeat[i].length; j++) {
                    perPeriodAndRepeat[i][j].addAll(other.perPeriodAndRepeat[i][j]);
                }
            }
            return this;
        }


        public DragstrLocusCases get(final int period, final int repeat) {
            return perPeriodAndRepeat[period - 1][repeat - 1];
        }

        public void addAll(final DragstrLocusCases traversalResult) {
            final int size = traversalResult.loci.size();
            for (int i = 0; i < size; i++) {
                final DragstrLocus locus = traversalResult.loci.get(i);
                final int depth = traversalResult.depth.getInt(i);
                final int nonref = traversalResult.nonref.getInt(i);
                final DragstrLocusCases[] periodResults = perPeriodAndRepeat[Math.min(perPeriodAndRepeat.length - 1, locus.getPeriod() - 1)];
                final DragstrLocusCases results = periodResults[Math.min(periodResults.length - 1, locus.getRepeats() - 1)];
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

    private DragstrLocusCases traverse(final SimpleInterval interval, final ReadsDataSource readsDataSource, final ReferenceDataSource referenceDataSource, final File lociDir) {

            final Stream<GATKRead> readStream = getTransformedReadStream(readsDataSource, interval, makeReadFilter());
            final Stream<DragstrLocus> lociStream = getSubjectStream(interval, lociDir);
            final Iterator<GATKRead> readIterator = readStream.iterator();
            final Iterator<DragstrLocus> subjectIterator = lociStream.iterator();
            final IntervalBuffer<GATKRead> readBuffer = new IntervalBuffer<>();
            final DragstrLocusCases result = new DragstrLocusCases();
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
                final LocusResult locusResult = evaluateLocus(lastSubject, reads, referenceDataSource);
                if (locusResult.qualifies) {
                    if (locusResult.nonrefDepth > 0) okCaseWithNonref++;
                    refReads += locusResult.depth - locusResult.nonrefDepth;
                    nonRefReads += locusResult.nonrefDepth;
                    result.add(lastSubject, locusResult.depth, locusResult.nonrefDepth);
                }
                caseCounts[locusResult.code.ordinal()]++;
                lastSubject = subjectIterator.hasNext() ? subjectIterator.next() : null;
            }
            //if (logger.isDebugEnabled()) {
            //    logger.debug("Summary data collection " + interval + ":\n"
            //            + "    " + Arrays.stream(LocusResultCode.values()).map(c -> "" + c + " = " + caseCounts[c.ordinal()]).collect(Collectors.joining("\n    "))
            //            + "\n"
            //            + "    ref/non-ref reads = " + refReads + "/" + nonRefReads + "\n"
            //            + "    ref-only/non-ref-containing reads cases = " + (caseCounts[0] - okCaseWithNonref) + "/" + okCaseWithNonref + "\n");
            //}
            return result;
    }

    protected Stream<GATKRead> getTransformedReadStream(final ReadsDataSource source, final SimpleInterval interval, final ReadFilter filter) {
        // if has reads, return an transformed/filtered/transformed stream
        if (hasReads()) {
            final ReadTransformer preTransformer = makePreReadFilterTransformer();
            final ReadTransformer postTransformer = makePostReadFilterTransformer();

            return interval == null ? Utils.stream(source) : Utils.stream(source.query(interval))
                    .map(preTransformer)
                    .filter(filter)
                    .map(postTransformer);
        }
        // returns an empty Stream if there are no reads
        return Stream.empty();
    }


    private int readEnd(final GATKRead read) {
        return read.isUnmapped() ? read.getAssignedStart() : read.getEnd();
    }
    private Stream<DragstrLocus> getSubjectStream(final SimpleInterval interval, final File lociDir) {
        try {
            final DragstrLocus.BinaryTableIndex lociIndex = DragstrLocus.BinaryTableIndex.load(new File(lociDir, "all.idx").toString());
            return DragstrLocus.binaryReader(new File(lociDir, "all.bin").toString(), lociIndex,
                    dictionary.getSequenceIndex(interval.getContig()), interval.getStart(), interval.getEnd()).stream();
        } catch (final IOException ex) {
            throw new UserException.CouldNotReadInputFile(new File(lociDir, "all.idx"), ex);
        }
    }

    private LocusResult evaluateLocus(DragstrLocus locus, List<GATKRead> rawReads, final ReferenceDataSource referenceDataSource) {
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
            byte[] referenceSequence = referenceDataSource.queryAndPrefetch(repeatInterval).getBases();
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

    private static class PeriodRepeatCombo {
        public final int period;
        public final int repeat;

        private PeriodRepeatCombo(final int period, final int repeat) {
            this.period = period;
            this.repeat = repeat;
        }

        private static PeriodRepeatCombo of(final int period, final int repeat) {
            return new PeriodRepeatCombo(period, repeat);
        }
    }
}
