package org.broadinstitute.hellbender.tools.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.OverlapDetector;
import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CoverageAnalysisProgramGroup;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.CalledContigPloidyCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.LocatableCopyNumberPosteriorDistributionCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.LocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SimpleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CalledContigPloidy;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.LocatableCopyNumberPosteriorDistribution;
import org.broadinstitute.hellbender.tools.copynumber.gcnv.GermlineCNVNamingConstants;
import org.broadinstitute.hellbender.tools.copynumber.gcnv.IntegerCopyNumberState;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.tools.sv.collections.StructuralVariantDepthQualityCollection;
import org.broadinstitute.hellbender.tools.sv.records.StructuralVariantDepthQuality;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.utils.tsv.TableUtils;

import java.io.*;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;
import java.util.stream.Stream;

/**
 * Integrates copy number posteriors for a given set of structural variant records.
 *
 * <h3>Inputs</h3>
 *
 * <ul>
 *     <li>
 *         Table of structural variants
 *     </li>
 *     <li>
 *         Germline copy number posteriors
 *     </li>
 * </ul>
 *
 * <h3>Output</h3>
 *
 * <ul>
 *     <li>
 *         Table of structural variant copy number probabilities
 *     </li>
 * </ul>
 *
 * <h3>Usage example</h3>
 *
 * <pre>
 *     gatk EventCopyNumberPosteriors
 * </pre>
 *
 * @author Andrey Smirnov &lt;asmirnov@broadinstitute.org&gt;
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Collects read counts at specified intervals",
        oneLineSummary = "Collects read counts at specified intervals",
        programGroup = CoverageAnalysisProgramGroup.class
)
@DocumentedFeature
public final class EventCopyNumberPosteriors extends GATKTool {
    public static final String STRUCTURAL_VARIANT_TABLE_LONG_NAME = "table";
    public static final String COPY_NUMBER_POSTERIORS_LONG_NAME = "posteriors";
    public static final String CONTIG_PLOIDY_CALLS_LONG_NAME = "ploidy-calls-dir";
    public static final String MIN_SIZE_LONG_NAME = "min-size";
    public static final String MAX_SIZE_LONG_NAME = "max-size";
    public static final String COPY_NEUTRAL_PRIOR_LONG_NAME = "copy-neutral-prior";

    @Argument(
            doc = "Structural variants table",
            fullName = STRUCTURAL_VARIANT_TABLE_LONG_NAME
    )
    private File structuralVariantTableFile;

    @Argument(
            doc = "Germline copy number posteriors VCF",
            fullName = COPY_NUMBER_POSTERIORS_LONG_NAME
    )
    private List<File> copyNumberPosteriorsFiles;

    @Argument(
            doc = "Contig ploidy path file",
            fullName = CONTIG_PLOIDY_CALLS_LONG_NAME
    )
    private File contigPloidyCallsDir;

    @Argument(
            doc = "Output event probabilities table",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private File outputFile;

    @Argument(
            doc = "Min event size",
            fullName = MIN_SIZE_LONG_NAME,
            minValue = 0,
            maxValue = Integer.MAX_VALUE,
            optional = true
    )
    private int minEventSize = 0;

    @Argument(
            doc = "Max event size",
            fullName = MAX_SIZE_LONG_NAME,
            minValue = 0,
            maxValue = Integer.MAX_VALUE,
            optional = true
    )
    private int maxEventSize = Integer.MAX_VALUE;

    @Argument(
            doc = "Prior probability of neutral copy number in regions where copy number calls are not available.",
            fullName = COPY_NEUTRAL_PRIOR_LONG_NAME,
            minValue = 0,
            maxValue = 1
    )
    private double copyNeutralPrior = 0.99;

    @Override
    public void traverse() {
        validateArguments();
        final Collection<EventRecord> eventsCollection = EventRecordReader.readEventRecords(structuralVariantTableFile)
                .filter(e -> isValidSize(e, minEventSize, maxEventSize)).collect(Collectors.toList());
        final Map<String,Collection<EventRecord>> sampleIdEventMap = getSampleIdEventMap(eventsCollection);
        final Map<String,CalledContigPloidyCollection> samplePloidyCallsMap = readCohortPloidyCalls().stream()
                .collect(Collectors.toMap(c -> c.getMetadata().getSampleName(), c -> c));
        final Map<String, SitePosterior> eventIdPosteriorsMap = new HashMap<>();

        final Set<String> processedSamples = new HashSet<>(SVUtils.hashMapCapacity(copyNumberPosteriorsFiles.size()));
        for (final File file : copyNumberPosteriorsFiles) {
            logger.info("Loading copy number file: " + file.getAbsolutePath());
            final LocatableCopyNumberPosteriorDistributionCollection samplePosteriors = new LocatableCopyNumberPosteriorDistributionCollection(file);
            final int intervalSize = validatePosteriors(samplePosteriors);
            if (processedSamples.contains(samplePosteriors.getSampleName())) {
                throw new UserException("Duplicate sample \"" + samplePosteriors.getSampleName() + "\" found in file " + file.getAbsolutePath());
            }
            if (sampleIdEventMap.containsKey(samplePosteriors.getSampleName())) {
                logger.info("Adding sample " + samplePosteriors.getSampleName() + " to site posteriors...");
                addSampleToSitePosteriors(eventIdPosteriorsMap, samplePloidyCallsMap, eventsCollection, samplePosteriors, intervalSize);
            } else {
                logger.warn("Posteriors were provided for sample " + samplePosteriors.getSampleName() + " but no events were found in this sample.");
            }
            processedSamples.add(samplePosteriors.getSampleName());
        }

        logger.info("Computing site qualities...");
        final StructuralVariantDepthQualityCollection siteQualities = getSiteQualities(eventIdPosteriorsMap.values(), samplePloidyCallsMap, getBestAvailableSequenceDictionary());
        writeSiteQualities(siteQualities);
    }

    private boolean isValidSize(final Locatable loc, final int min, final int max) {
        final int size = loc.getLengthOnReference();
        return size >= min && size <= max;
    }

    private int validatePosteriors(final LocatableCopyNumberPosteriorDistributionCollection posteriors) {
        if (posteriors.getRecords().isEmpty()) {
            throw new UserException.BadInput("Copy number posteriors collection for sample " + posteriors.getSampleName() + " is empty");
        }
        final int size = posteriors.getRecords().stream().map(LocatableCopyNumberPosteriorDistribution::getLengthOnReference)
                .collect(Collectors.groupingBy(i -> i)).entrySet().stream()
                .sorted(Comparator.comparingInt(e -> -e.getValue().size())).map(Map.Entry::getKey)
                .collect(Collectors.toList()).get(0);
        logger.info("Automatically determined interval size to be: " + size);
        final long nonUniformIntervals = posteriors.getRecords().stream().filter(r -> r.getLengthOnReference() != size).count();
        if (nonUniformIntervals > 0) {
            logger.warn("There were " + nonUniformIntervals + " copy number intervals for sample " + posteriors.getSampleName() + " not of size " + size);
        }
        return size;
    }

    private void writeSiteQualities(final StructuralVariantDepthQualityCollection siteQualities) {
        siteQualities.write(outputFile);
    }

    private static StructuralVariantDepthQualityCollection getSiteQualities(final Collection<SitePosterior> sitePosteriors,
                                                                            final Map<String,CalledContigPloidyCollection> samplePloidyCallsMap,
                                                                            final SAMSequenceDictionary dictionary) {
        final List<StructuralVariantDepthQuality> qualities = sitePosteriors.stream().map(p -> getSiteQuality(p, samplePloidyCallsMap)).collect(Collectors.toList());
        final LocatableMetadata metadata = new SimpleLocatableMetadata(dictionary);
        return new StructuralVariantDepthQualityCollection(metadata, qualities);
    }

    // Gets probability of depth support in at least one sample
    private static StructuralVariantDepthQuality getSiteQuality(final SitePosterior sitePosterior,
                                                                final Map<String,CalledContigPloidyCollection> samplePloidyCallsMap) {
        final Collection<LocatableCopyNumberPosteriorDistribution> samplePosteriors = sitePosterior.getPosteriors();
        final String contig = sitePosterior.getEvent().getInterval().getContig();
        double siteP = 1;
        for (final LocatableCopyNumberPosteriorDistribution samplePosterior : samplePosteriors) {
            final String sample = samplePosterior.getSample();
            final int ploidy = getSamplePloidy(sample, contig, samplePloidyCallsMap);
            double p = 0;
            for (final IntegerCopyNumberState state : samplePosterior.getIntegerCopyNumberStateList()) {
                if (isEventState(state, sitePosterior.getEvent(), ploidy)) {
                    p += FastMath.exp(samplePosterior.getCopyNumberPosterior(state));
                }
            }
            siteP = siteP * (1.0 - p);
        }
        // Probability that there is no supporting depth signal in any of the samples
        siteP = Math.max(siteP, Double.MIN_VALUE);
        final double qual =  Math.min(-10.0 * Math.log10(siteP), 99);
        return new StructuralVariantDepthQuality(sitePosterior.getEvent(), qual);
    }

    private static int getSamplePloidy(final String sample, final String contig,
                                       final Map<String,CalledContigPloidyCollection> samplePloidyCallsMap) {
        if (!samplePloidyCallsMap.containsKey(sample)) {
            throw new UserException("Could not find ploidy calls for sample: " + sample);
        }
        final Optional<CalledContigPloidy> ploidy = samplePloidyCallsMap.get(sample).getRecords().stream()
                .filter(record -> record.getContig().equals(contig)).findFirst();
        if (!ploidy.isPresent()) {
            throw new UserException("Sample " + sample + " does not have a ploidy call for contig: " + contig);
        }
        return ploidy.get().getPloidy();
    }

    private static boolean isEventState(final IntegerCopyNumberState state, final EventRecord event, final int ploidy) {
        if (event.getType().equals("DEL")) {
            return state.getCopyNumber() < ploidy;
        }
        if (event.getType().equals("DUP")) {
            return state.getCopyNumber() > ploidy;
        }
        return false;
    }

    private Collection<CalledContigPloidyCollection> readCohortPloidyCalls() {
        final File[] subdirs = contigPloidyCallsDir.listFiles(File::isDirectory);
        final Collection<CalledContigPloidyCollection> ploidyCollections = new ArrayList<>();
        for (final File subdir : subdirs) {
            if (subdir.getName().startsWith(GermlineCNVNamingConstants.SAMPLE_PREFIX)) {
                ploidyCollections.add(readPloidyCalls(subdir));
            }
        }
        return ploidyCollections;
    }

    private static CalledContigPloidyCollection readPloidyCalls(final File dir) {
        final Path ploidyPath = Paths.get(dir.getAbsolutePath(), GermlineCNVNamingConstants.CONTIG_PLOIDY_FILE_NAME);
        final File ploidyFile = new File(ploidyPath.toString());
        return new CalledContigPloidyCollection(ploidyFile);
    }

    private void addSampleToSitePosteriors(final Map<String, SitePosterior> eventIdPosteriorsMap,
                                           final Map<String,CalledContigPloidyCollection> samplePloidyCallsMap,
                                           final Collection<EventRecord> eventsToAdd,
                                           final LocatableCopyNumberPosteriorDistributionCollection samplePosteriors,
                                           final int intervalSize) {
        final String sample = samplePosteriors.getSampleName();
        if (samplePosteriors.getRecords().isEmpty()) {
            throw new IllegalArgumentException("Sample posteriors collection was empty");
        }
        final List<IntegerCopyNumberState> copyStates = samplePosteriors.getRecords().get(0).getIntegerCopyNumberStateList();
        if (copyStates.isEmpty()) {
            throw new IllegalArgumentException("The first posteriors record had no copy number states");
        }
        final int maxCopyState = copyStates.stream().mapToInt(IntegerCopyNumberState::getCopyNumber).max().getAsInt();
        final OverlapDetector<LocatableCopyNumberPosteriorDistribution> overlapDetector = samplePosteriors.createNewOverlapDetector();
        for (final EventRecord event : eventsToAdd) {
            if (event.isInSample(sample)) {
                final int eventPloidy = getSamplePloidy(sample, event.getInterval().getContig(), samplePloidyCallsMap);
                final LocatableCopyNumberPosteriorDistribution eventPosterior = getEventPosterior(event, sample, overlapDetector, maxCopyState, intervalSize, eventPloidy);
                eventIdPosteriorsMap.putIfAbsent(event.getId(), new SitePosterior(event));
                eventIdPosteriorsMap.get(event.getId()).getPosteriors().add(eventPosterior);
            }
        }
    }

    private LocatableCopyNumberPosteriorDistribution getEventPosterior(final EventRecord event,
                                                                       final String sample,
                                                                       final OverlapDetector<LocatableCopyNumberPosteriorDistribution> overlapDetector,
                                                                       final int maxCopyState,
                                                                       final int posteriorsIntervalSize,
                                                                       final int ploidy) {
        Utils.validateArg(maxCopyState >= 0, "Maximum copy state must be non-negative");
        final SimpleInterval eventInterval = event.getInterval();
        final Set<LocatableCopyNumberPosteriorDistribution> overlappers = overlapDetector.getOverlaps(eventInterval);
        final double[] copyStateSums = new double[maxCopyState + 1];
        Arrays.fill(copyStateSums, Double.MIN_VALUE);
        int overlapSize = 0;
        for (final LocatableCopyNumberPosteriorDistribution overlapper : overlappers) {
            final int overlap = overlapper.getInterval().intersect(eventInterval).size();
            // TODO requires 50% overlap with interval; make adjustable
            if (overlap >= 0.5 * overlapper.getInterval().size()) {
                overlapSize += overlap;
                for (final IntegerCopyNumberState state : overlapper.getIntegerCopyNumberStateList()) {
                    copyStateSums[state.getCopyNumber()] += overlapper.getCopyNumberPosterior(state);
                }
            }
        }

        // Fill in missing copy number posterior intervals with a prior
        final double unsupportedIntervals = (event.getInterval().size() - overlapSize) / (double) posteriorsIntervalSize;
        if (unsupportedIntervals > 0) {
            final double logNeutralProb = FastMath.log(copyNeutralPrior);
            final double logNonNeutralProb = FastMath.log((1.0 - copyNeutralPrior) / (copyStateSums.length - 1));
            for (int i = 0; i < copyStateSums.length; i++) {
                if (i != ploidy) {
                    copyStateSums[i] += logNonNeutralProb * unsupportedIntervals;
                } else {
                    copyStateSums[i] += logNeutralProb * unsupportedIntervals;
                }
            }
        }

        double denom = 0;
        final double maxStateSum = DoubleStream.of(copyStateSums).max().getAsDouble();
        for (int i = 0; i < copyStateSums.length; i++) {
            // Normalize to avoid underflow error
            copyStateSums[i] -= maxStateSum;
            denom += FastMath.exp(copyStateSums[i]);
        }
        final double logDenom = Math.log(denom);
        final Map<IntegerCopyNumberState,Double> eventPosterior = new HashMap<>(SVUtils.hashMapCapacity(maxCopyState + 1));
        for (int i = 0; i < copyStateSums.length; i++) {
            final Double p = copyStateSums[i] - logDenom;
            eventPosterior.put(new IntegerCopyNumberState(i), p);
        }
        return new LocatableCopyNumberPosteriorDistribution(eventPosterior, sample, event.getInterval());
    }

    private static Map<String,Collection<EventRecord>> getSampleIdEventMap(final Collection<EventRecord> events) {
        final Map<String,Collection<EventRecord>> map = new HashMap<>();
        for (final EventRecord record : events) {
            for (final String sample : record.getSamples()) {
                map.putIfAbsent(sample, new ArrayList<>());
                map.get(sample).add(record);
            }
        }
        return map;
    }

    private static Map<String,Collection<String>> getEventIdSamplesMap(final Collection<EventRecord> events) {
        final Map<String,Collection<String>> map = new HashMap<>(SVUtils.hashMapCapacity(events.size()));
        for (final EventRecord record : events) {
            if (map.containsKey(record.getId())) {
                throw new UserException("Duplicate event id found: " + record.getId());
            }
            map.put(record.getId(), record.getSamples());
        }
        return map;
    }

    private static class SitePosterior {
        private EventRecord event;
        private Collection<LocatableCopyNumberPosteriorDistribution> posteriors;

        public SitePosterior(final EventRecord event) {
            this.event = event;
            this.posteriors = new ArrayList<>(event.getSamples().size());
        }

        public EventRecord getEvent() {
            return event;
        }

        public Collection<LocatableCopyNumberPosteriorDistribution> getPosteriors() {
            return posteriors;
        }
    }

    private final static class EventRecordReader {

        private enum TABLE_COLUMNS {
            CHROM,
            START,
            END,
            ID,
            SAMPLES,
            TYPE
        }

        protected static Stream<EventRecord> readEventRecords(final File file) {
            try {
                final Reader reader = new FileReader(file);
                final TableReader<EventRecord> tableReader = TableUtils.reader(reader, EventRecordReader::eventRecordFactory);
                return tableReader.stream();
            } catch (final FileNotFoundException e) {
                throw new UserException(e.getMessage());
            } catch (final IOException e) {
                throw new UserException(e.getMessage());
            }
        }

        private static Function<DataLine, EventRecord> eventRecordFactory(final TableColumnCollection columns,
                                                                         final Function<String, RuntimeException> errorFunction) {
            return EventRecordReader::eventRecordParser;
        }

        private static EventRecord eventRecordParser(final DataLine data) {
            final String chrom = data.get(TABLE_COLUMNS.CHROM);
            final int start = Integer.valueOf(data.get(TABLE_COLUMNS.START));
            final int end = Integer.valueOf(data.get(TABLE_COLUMNS.END));
            final SimpleInterval interval = new SimpleInterval(chrom, start, end);
            final String[] samples = data.get(TABLE_COLUMNS.SAMPLES).split(",");
            return new EventRecord(data.get(TABLE_COLUMNS.ID), data.get(TABLE_COLUMNS.TYPE), interval, samples);
        }
    }

    private void validateArguments() {
        Utils.validateArg(minEventSize <= maxEventSize,
                "Minimum event size cannot be larger than maximum event size.");
        Utils.validateArg(getBestAvailableSequenceDictionary() != null, "Sequence dictionary is required.");
    }

    public final static class EventRecord implements Locatable {
        private final String id;
        private final String type;
        private final SimpleInterval interval;
        private final Set<String> samples;

        public EventRecord(final String id, final String type, final SimpleInterval interval, final String[] samples) {
            this.id = id;
            this.type = type;
            this.interval = interval;
            this.samples = new HashSet<>(Arrays.asList(samples));
        }

        public String getId() {
            return id;
        }

        public String getType() {
            return type;
        }

        public SimpleInterval getInterval() {
            return interval;
        }

        public Collection<String> getSamples() {
            return samples;
        }

        public boolean isInSample(final String id) {
            return samples.contains(id);
        }

        @Override
        public String getContig() {
            return interval.getContig();
        }

        @Override
        public int getStart() {
            return interval.getStart();
        }

        @Override
        public int getEnd() {
            return interval.getEnd();
        }
    }
}
