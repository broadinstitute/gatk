package org.broadinstitute.hellbender.tools.sv;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CoverageAnalysisProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.CopyNumberPosteriorDistributionCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleMetadata;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.utils.tsv.TableUtils;
import scala.Tuple2;

import java.io.*;
import java.util.*;
import java.util.function.Function;
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
public final class EventCopyNumberPosteriors extends CommandLineProgram {
    public static final String STRUCTURAL_VARIANT_TABLE_LONG_NAME = "table";
    public static final String COPY_NUMBER_POSTERIORS_LONG_NAME = "posteriors";
    public static final String CONTIG_PLOIDY_CALLS_LONG_NAME = "ploidy-calls";
    public static final String MIN_SIZE_LONG_NAME = "min-size";
    public static final String MAX_SIZE_LONG_NAME = "max-size";

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
            doc = "Contig ploidy calls",
            fullName = CONTIG_PLOIDY_CALLS_LONG_NAME
    )
    private File contigPloidyCallsFile;

    @Argument(
            doc = "Output event probabilities table",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private File outputPath;

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

    @Override
    public Object doWork() {
        validateArguments();
        final Collection<EventRecord> eventsCollection = EventRecordReader.readEventRecords(structuralVariantTableFile);
        final Map<String,Collection<EventRecord>> sampleIdEventMap = getSampleIdEventMap(eventsCollection);
        final Map<String,Collection<String>> eventIdSamplesMap = getEventIdSamplesMap(eventsCollection);

        final Set<String> samples = sampleIdEventMap.keySet();
        final Map<String,Collection<String>> eventIdSampleMap = new HashMap<>();

        final Set<String> processedSamples = new HashSet<>(SVUtils.hashMapCapacity(copyNumberPosteriorsFiles.size()));
        for (final File file : copyNumberPosteriorsFiles) {
            final CopyNumberPosteriorDistributionCollection posteriors = new CopyNumberPosteriorDistributionCollection(file);
            final SampleMetadata metadata = posteriors.getMetadata();
            if (processedSamples.contains(metadata.getSampleName())) {
                throw new UserException("Duplicate sample \"" + metadata.getSampleName() + "\" found in file " + file.getAbsolutePath());
            }
            if (sampleIdEventMap.containsKey(metadata.getSampleName())) {

            } else {
                logger.warn("Posteriors were provided for sample " + metadata.getSampleName() + " but no events for this sample were found.");
            }
            processedSamples.add(metadata.getSampleName());
        }
        return null;
    }

    private static void addSampleToCohortPosteriors(final Map<String,CohortPosterior> cohortPosteriors,
                                                    Collection<EventRecord> eventsToAdd,
                                                    final CopyNumberPosteriorDistributionCollection samplePosteriors) {
        final String sample = samplePosteriors.getMetadata().getSampleName();
        for (final EventRecord event : eventsToAdd) {
            final SimpleInterval interval = event.getInterval();
            samplePosteriors.getRecords().get(0).
        }
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
            map.put(record.getId(), Arrays.asList(record.getSamples()));
        }
        return map;
    }

    private static class CohortPosterior {
        private String eventId;
        private Collection<SamplePosterior> posteriors;

        public CohortPosterior(final String eventId, final int numSamples) {
            this.eventId = eventId;
            this.posteriors = new ArrayList<>(numSamples);
        }

        public String getEventId() {
            return eventId;
        }

        public Collection<SamplePosterior> getPosteriors() {
            return posteriors;
        }
    }

    private static class SamplePosterior {
        private final String sample;
        private final double pDeletion;
        private final double pNeutral;

        public SamplePosterior(final String sample, final double pDeletion, final double pNeutral) {
            this.sample = sample;
            this.pDeletion = pDeletion;
            this.pNeutral = pNeutral;
        }

        public String getSample() {
            return sample;
        }

        public double getDeletionProbability() {
            return pDeletion;
        }

        public double getNeutralProbability() {
            return pNeutral;
        }

        public double getDuplicationProbability() {
            return 1.0 - pDeletion - pNeutral;
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

        protected static Collection<EventRecord> readEventRecords(final File file) {
            try {
                final Reader reader = new FileReader(file);
                final TableReader<EventRecord> tableReader = TableUtils.reader(reader, EventRecordReader::eventRecordFactory);
                return tableReader.toList();
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
    }

    private final static class EventRecord {
        private final String id;
        private final String type;
        private final SimpleInterval interval;
        private final String[] samples;

        public EventRecord(final String id, final String type, final SimpleInterval interval, final String[] samples) {
            this.id = id;
            this.type = type;
            this.interval = interval;
            this.samples = samples;
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

        public String[] getSamples() {
            return samples;
        }

        public Stream<Tuple2<String,EventRecord>> sampleTuples() {
            return Arrays.stream(samples).map(s -> new Tuple2<>(s, this));
        }
    }
}
