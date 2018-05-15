package org.broadinstitute.hellbender.utils.read.markduplicates;

import htsjdk.samtools.*;
import htsjdk.samtools.DuplicateScoringStrategy.ScoringStrategy;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Histogram;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import picard.sam.markduplicates.util.OpticalDuplicateFinder;
import picard.sam.util.PhysicalLocation;

import java.io.File;
import java.util.*;

/**
 * Abstract class that holds parameters and methods common to classes that perform duplicate
 * detection and/or marking within SAM/BAM/CRAM files.
 *
 * @author Nils Homer
 */
public abstract class AbstractMarkDuplicatesCommandLineProgram extends AbstractOpticalDuplicateFinderCommandLineProgram {

    @Argument(fullName = StandardArgumentDefinitions.INPUT_LONG_NAME, shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME,
            doc = "One or more input SAM/BAM/CRAM files to analyze. Must be coordinate sorted.")
    public List<File> INPUT;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "The output file to write marked records to")
    public File OUTPUT;

    @Argument(shortName = StandardArgumentDefinitions.METRICS_FILE_SHORT_NAME,
            fullName = StandardArgumentDefinitions.METRICS_FILE_LONG_NAME,
            doc = "File to write duplication metrics to")
    public File METRICS_FILE;

    @Argument(shortName = StandardArgumentDefinitions.PROGRAM_RECORD_ID_SHORT_NAME,
            doc = "The program record ID for the @PG record(s) created by this program. Set to null to disable " +
                    "PG record creation.  This string may have a suffix appended to avoid collision with other " +
                    "program record IDs.",
            optional = true)
    public String PROGRAM_RECORD_ID = "MarkDuplicatesGATK";

    @Argument(shortName = "PG_VERSION",
            doc = "Value of VN tag of PG record to be created. If not specified, the version will be detected automatically.",
            optional = true)
    public String PROGRAM_GROUP_VERSION;

    @Argument(shortName = "PG_COMMAND",
            doc = "Value of CL tag of PG record to be created. If not supplied the command line will be detected automatically.",
            optional = true)
    public String PROGRAM_GROUP_COMMAND_LINE;

    @Argument(shortName = "PG_NAME",
            doc = "Value of PN tag of PG record to be created.")
    public String PROGRAM_GROUP_NAME = getClass().getSimpleName();

    @Argument(shortName = "CO",
            doc = "Comment(s) to include in the output file's header.",
            optional = true)
    public List<String> COMMENT = new ArrayList<>();

    @Argument(doc = "If true do not write duplicates to the output file instead of writing them with appropriate flags set.")
    public boolean REMOVE_DUPLICATES = false;

    @Argument(shortName = StandardArgumentDefinitions.ASSUME_SORTED_SHORT_NAME,
            doc = "If true, assume that the input file is coordinate sorted even if the header says otherwise.")
    public boolean ASSUME_SORTED = false;

    @Argument(shortName = StandardArgumentDefinitions.DUPLICATE_SCORING_STRATEGY_SHORT_NAME,
            fullName = StandardArgumentDefinitions.DUPLICATE_SCORING_STRATEGY_LONG_NAME,
            doc = "The scoring strategy for choosing the non-duplicate among candidates.")
    public DuplicateScoringStrategy.ScoringStrategy DUPLICATE_SCORING_STRATEGY = ScoringStrategy.TOTAL_MAPPED_REFERENCE_LENGTH;

    /** The program groups that have been seen during the course of examining the input records. */
    protected final Set<String> pgIdsSeen = new LinkedHashSet<>();

    /**
     * We have to re-chain the program groups based on this algorithm.  This returns the map from existing program group ID
     * to new program group ID.
     */
    protected Map<String, String> getChainedPgIds(final SAMFileHeader outputHeader) {
        final Map<String, String> chainedPgIds;
        // Generate new PG record(s)
        if (PROGRAM_RECORD_ID != null) {
            final PgIdGenerator pgIdGenerator = new PgIdGenerator(outputHeader);
            if (PROGRAM_GROUP_VERSION == null) {
                PROGRAM_GROUP_VERSION = this.getVersion();
            }
            if (PROGRAM_GROUP_COMMAND_LINE == null) {
                PROGRAM_GROUP_COMMAND_LINE = this.getCommandLine();
            }
            chainedPgIds = new LinkedHashMap<>();
            for (final String existingId : this.pgIdsSeen) {
                final String newPgId = pgIdGenerator.getNonCollidingId(PROGRAM_RECORD_ID);
                chainedPgIds.put(existingId, newPgId);
                final SAMProgramRecord programRecord = new SAMProgramRecord(newPgId);
                programRecord.setProgramVersion(PROGRAM_GROUP_VERSION);
                programRecord.setCommandLine(PROGRAM_GROUP_COMMAND_LINE);
                programRecord.setProgramName(PROGRAM_GROUP_NAME);
                programRecord.setPreviousProgramGroupId(existingId);
                outputHeader.addProgramRecord(programRecord);
            }
        } else {
            chainedPgIds = null;
        }
        return chainedPgIds;
    }

    /**
     * Writes the metrics given by the libraryIdGenerator to the METRICS_FILE.
     *
     * @param libraryIdGenerator
     */
    protected void finalizeAndWriteMetrics(final LibraryIdGenerator libraryIdGenerator) {
        //We want to sort libraries by name
        final SortedMap<String, DuplicationMetrics> metricsByLibrary = new TreeMap<>(Utils.COMPARE_STRINGS_NULLS_FIRST);
        metricsByLibrary.putAll(libraryIdGenerator.getMetricsByLibraryMap());

        final Histogram<Short> opticalDuplicatesByLibraryId = libraryIdGenerator.getOpticalDuplicatesByLibraryIdMap();
        final Map<String, Short> libraryIds = libraryIdGenerator.getLibraryIdsMap();

        // Write out the metrics
        final MetricsFile<DuplicationMetrics, Double> file = getMetricsFile();
        for (final Map.Entry<String, DuplicationMetrics> entry : metricsByLibrary.entrySet()) {
            final String libraryName = entry.getKey();
            final DuplicationMetrics metrics = entry.getValue();

            metrics.READ_PAIRS_EXAMINED = metrics.READ_PAIRS_EXAMINED / 2;
            metrics.READ_PAIR_DUPLICATES = metrics.READ_PAIR_DUPLICATES / 2;

            // Add the optical dupes to the metrics
            final Short libraryId = libraryIds.get(libraryName);
            if (libraryId != null) {
                @SuppressWarnings("unchecked")//type-checker being annoying here
                final Histogram.Bin<Short> bin = opticalDuplicatesByLibraryId.get(libraryId);
                if (bin != null) {
                    metrics.READ_PAIR_OPTICAL_DUPLICATES = (long) bin.getValue();
                }
            }
            metrics.calculateDerivedMetrics();
            file.addMetric(metrics);
        }

        if (metricsByLibrary.size() == 1) {
            file.setHistogram(metricsByLibrary.values().iterator().next().calculateRoiHistogram());
        }

        file.write(METRICS_FILE);
    }

    /** Little class to generate program group IDs */
    static class PgIdGenerator {
        private int recordCounter;

        private final Set<String> idsThatAreAlreadyTaken = new LinkedHashSet<>();

        PgIdGenerator(final SAMFileHeader header) {
            for (final SAMProgramRecord pgRecord : header.getProgramRecords()) {
                idsThatAreAlreadyTaken.add(pgRecord.getProgramGroupId());
            }
            recordCounter = idsThatAreAlreadyTaken.size();
        }

        String getNonCollidingId(final String recordId) {
            if (!idsThatAreAlreadyTaken.contains(recordId)) {
                // don't remap 1st record. If there are more records
                // with this id, they will be remapped in the 'else'.
                idsThatAreAlreadyTaken.add(recordId);
                ++recordCounter;
                return recordId;
            } else {
                String newId;
                // Below we tack on one of roughly 1.7 million possible 4 digit base36 at random. We do this because
                // our old process of just counting from 0 upward and adding that to the previous id led to 1000s of
                // calls idsThatAreAlreadyTaken.contains() just to resolve 1 collision when merging 1000s of similarly
                // processed bams.
                while (idsThatAreAlreadyTaken.contains(newId = recordId + "." + SamFileHeaderMerger.positiveFourDigitBase36Str(recordCounter++)))
                    ;

                idsThatAreAlreadyTaken.add(newId);
                return newId;
            }

        }
    }

    /** Little class used to package up a header and an iterable/iterator. */
    public static final class SamHeaderAndIterator implements AutoCloseable{
        public final SAMFileHeader header;
        public final CloseableIterator<SAMRecord> iterator;
        private final List<SamReader> readers;

        public SamHeaderAndIterator(final SAMFileHeader header, final CloseableIterator<SAMRecord> iterator, final List<SamReader> readers) {
            this.header = header;
            this.iterator = iterator;
            this.readers = readers;
        }

        @Override
        public void close(){
            CloserUtil.close(readers);
            CloserUtil.close(iterator);
        }
    }

    /**
     * Since this may read it's inputs more than once this method does all the opening
     * and checking of the inputs.
     */
    protected SamHeaderAndIterator openInputs() {
        final List<SAMFileHeader> headers = new ArrayList<>(INPUT.size());
        final List<SamReader> readers = new ArrayList<>(INPUT.size());

        for (final File f : INPUT) {
            final SamReader reader = SamReaderFactory.makeDefault()
                    .enable(SamReaderFactory.Option.EAGERLY_DECODE)
                    .referenceSequence(REFERENCE_SEQUENCE)
                    .open(f); // eager decode
            final SAMFileHeader header = reader.getFileHeader();

            if (!ASSUME_SORTED && header.getSortOrder() != SAMFileHeader.SortOrder.coordinate) {
                throw new UserException("Input file " + f.getAbsolutePath() + " is not coordinate sorted.");
            }

            headers.add(header);
            readers.add(reader);
        }

        if (headers.size() == 1) {
            return new SamHeaderAndIterator(headers.get(0), readers.get(0).iterator(), readers);
        } else {
            final SamFileHeaderMerger headerMerger = new SamFileHeaderMerger(SAMFileHeader.SortOrder.coordinate, headers, false);
            final MergingSamRecordIterator iterator = new MergingSamRecordIterator(headerMerger, readers, ASSUME_SORTED);
            return new SamHeaderAndIterator(headerMerger.getMergedHeader(), iterator, readers);
        }
    }

    /**
     * Looks through the set of reads and identifies how many of the duplicates are
     * in fact optical duplicates, and stores the data in the instance level histogram.
     */
    public static void trackOpticalDuplicates(List<? extends ReadEnds> ends,
                                              final ReadEnds keeper,
                                              final OpticalDuplicateFinder opticalDuplicateFinder,
                                              final LibraryIdGenerator libraryIdGenerator) {
        boolean hasFR = false, hasRF = false;

        // Check to see if we have a mixture of FR/RF
        for (final ReadEnds end : ends) {
            if (ReadEnds.FR == end.orientationForOpticalDuplicates) {
                hasFR = true;
            } else if (ReadEnds.RF == end.orientationForOpticalDuplicates) {
                hasRF = true;
            }
        }

        // Check if we need to partition since the orientations could have changed
        if (hasFR && hasRF) { // need to track them independently
            // Variables used for optical duplicate detection and tracking
            final List<ReadEnds> trackOpticalDuplicatesF = new ArrayList<>();
            final List<ReadEnds> trackOpticalDuplicatesR = new ArrayList<>();

            // Split into two lists: first of pairs and second of pairs, since they must have orientation and same starting end
            for (final ReadEnds end : ends) {
                if (ReadEnds.FR == end.orientationForOpticalDuplicates) {
                    trackOpticalDuplicatesF.add(end);
                } else if (ReadEnds.RF == end.orientationForOpticalDuplicates) {
                    trackOpticalDuplicatesR.add(end);
                } else {
                    throw new UserException("Found an unexpected orientation: " + end.orientation);
                }
            }

            // track the duplicates
            trackOpticalDuplicates(trackOpticalDuplicatesF, opticalDuplicateFinder, libraryIdGenerator.getOpticalDuplicatesByLibraryIdMap(), keeper);
            trackOpticalDuplicates(trackOpticalDuplicatesR, opticalDuplicateFinder, libraryIdGenerator.getOpticalDuplicatesByLibraryIdMap(), keeper);
        } else { // No need to partition
            AbstractMarkDuplicatesCommandLineProgram.trackOpticalDuplicates(ends, opticalDuplicateFinder, libraryIdGenerator.getOpticalDuplicatesByLibraryIdMap(), keeper);
        }
    }

    /**
     * Looks through the set of reads and identifies how many of the duplicates are
     * in fact optical duplicates, and stores the data in the instance level histogram.
     */
    private static void trackOpticalDuplicates(final List<? extends PhysicalLocation> list,
                                               final OpticalDuplicateFinder opticalDuplicateFinder,
                                               final Histogram<Short> opticalDuplicatesByLibraryId,
                                               final ReadEnds keeper) {
        final boolean[] opticalDuplicateFlags = opticalDuplicateFinder.findOpticalDuplicates(list, keeper);

        int opticalDuplicates = 0;
        for (final boolean b : opticalDuplicateFlags) if (b) ++opticalDuplicates;
        if (opticalDuplicates > 0) {
            opticalDuplicatesByLibraryId.increment(list.get(0).getLibraryId(), opticalDuplicates);
        }
    }
}
