package org.broadinstitute.hellbender.tools.picard.markduplicates;

import htsjdk.samtools.*;
import htsjdk.samtools.util.*;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;
import org.broadinstitute.hellbender.utils.iterators.MarkDuplicatesWithMateCigarIterator;
import org.broadinstitute.hellbender.utils.sam.markduplicates.AbstractMarkDuplicatesCommandLineProgram;

import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * An even better duplication marking algorithm that handles all cases including clipped
 * and gapped alignments.
 * <p/>
 * This tool differs with MarkDuplicates as it may break ties differently.  Furthermore,
 * as it is a one-pass algorithm, it cannot know the program records contained in the file
 * that should be chained in advance.  Therefore it will only be able to examine the header
 * to attempt to infer those program group records that have no associated previous program
 * group record. If a read is encountered without a program record, or not one as previously
 * defined, it will not be updated.
 * <p/>
 * This tool will also not work with alignments that have large gaps or skips, such as those
 * from RNA-seq data.  This is due to the need to buffer small genomic windows to ensure
 * integrity of the duplicate marking, while large skips (ex. skipping introns) in the
 * alignment records would force making that window very large, thus exhausting memory.
 *
 * @author Nils Homer
 */
@CommandLineProgramProperties(
        usage = "Examines aligned records in the supplied SAM or BAM file to locate duplicate molecules. " +
                "All records are then written to the output file with the duplicate records flagged.",
        usageShort = "Examines aligned records in the supplied SAM or BAM file to locate duplicate molecules.",
        programGroup = ReadProgramGroup.class
)
public class MarkDuplicatesWithMateCigar extends AbstractMarkDuplicatesCommandLineProgram {
    private final Log log = Log.getInstance(MarkDuplicatesWithMateCigar.class);

    @Argument(doc = "The minimum distance to buffer records to account for clipping on the 5' end of the records." +
            "Set this number to -1 to use twice the first read's read length (or 100, whichever is smaller).", optional = true)
    public int MINIMUM_DISTANCE = -1;

    @Argument(doc = "Skip record pairs with no mate cigar and include them in the output.")
    boolean SKIP_PAIRS_WITH_NO_MATE_CIGAR = true;

    @Argument(doc = "The block size for use in the coordinate-sorted record buffer.", optional = true)
    public int BLOCK_SIZE = 100000;

    /** Warnings that will only be emitted once */
    private boolean warnedNullProgramRecords = false;
    private boolean warnedMissingProgramRecords = false;

    /**
     * Main work method.
     */
    protected Object doWork() {
        IOUtil.assertFilesAreReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);
        IOUtil.assertFileIsWritable(METRICS_FILE);

        // Open the inputs
        final SamHeaderAndIterator headerAndIterator = openInputs();
        final SAMFileHeader header = headerAndIterator.header;

        // Create the output header
        final SAMFileHeader outputHeader = header.clone();
        outputHeader.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        for (final String comment : COMMENT) outputHeader.addComment(comment);

        // Since this is one-pass, unlike MarkDuplicates, we cannot only chain together program
        // group records we have seen, we have to assume all of them may be seen.  We can perhaps
        // filter out any program groups which have been referenced previously.
        setPGIdsSeen(outputHeader);
        // Key: previous PG ID on a SAM Record (or null).  Value: New PG ID to replace it.
        final Map<String, String> chainedPgIds = getChainedPgIds(outputHeader);

        // Open the output
        final SAMFileWriter out = new SAMFileWriterFactory().makeSAMOrBAMWriter(outputHeader,
                true,
                OUTPUT);

        // Create the mark duplicate iterator.  The duplicate marking is handled by the iterator, conveniently.
        final MarkDuplicatesWithMateCigarIterator iterator = new MarkDuplicatesWithMateCigarIterator(headerAndIterator.header,
                headerAndIterator.iterator,
                this.opticalDuplicateFinder,
                this.DUPLICATE_SCORING_STRATEGY,
                this.MINIMUM_DISTANCE,
                this.REMOVE_DUPLICATES,
                this.SKIP_PAIRS_WITH_NO_MATE_CIGAR,
                this.MAX_RECORDS_IN_RAM,
                this.BLOCK_SIZE,
                this.TMP_DIR);

        // progress logger!
        final ProgressLogger progress = new ProgressLogger(log, (int) 1e6, "Read");

        // Go through the records
        for (final SAMRecord record : new IterableAdapter<SAMRecord>(iterator)) {
            if (progress.record(record)) {
                iterator.logMemoryStats(log);
            }

            // Update the program record if necessary
            updateProgramRecord(record, chainedPgIds);

            // Write the alignment
            out.addAlignment(record);
        }

        // remember to close the inputs
        iterator.close();

        out.close();

        // For convenience to reference
        final Histogram<Short> opticalDupesByLibraryId = iterator.getOpticalDupesByLibraryId();

        // Log info
        log.info("Processed " + progress.getCount() + " records");
        log.info("Found " + iterator.getNumRecordsWithNoMateCigar() + " records with no mate cigar optional tag.");
        log.info("Marking " + iterator.getNumDuplicates() + " records as duplicates.");
        log.info("Found " + ((long) opticalDupesByLibraryId.getSumOfValues()) + " optical duplicate clusters."); // cast as long due to returning a double

        // Write out the metrics
        finalizeAndWriteMetrics(iterator.getLibraryIdGenerator());

        return null;
    }

    /**
     * Updates the program record if necessary.
     */
    private void updateProgramRecord(final SAMRecord record, final Map<String, String> chainedPgIds) {
        if (PROGRAM_RECORD_ID != null) {
            final String pgId = record.getStringAttribute(SAMTag.PG.name());
            if (null == pgId) {
                if (!warnedNullProgramRecords) {
                    warnedNullProgramRecords = true;
                    log.warn("Encountered a record with no program record, program group chaining will not occur for this read: " + record);
                } // else already warned!
            } else if (!chainedPgIds.containsKey(pgId)) {
                if (!warnedMissingProgramRecords) {
                    warnedMissingProgramRecords = true;
                    log.warn("Encountered a record with an intermediate program record, program group chaining will not occur for this read: " + record);
                } // else already warned!
            } else {
                record.setAttribute(SAMTag.PG.name(), chainedPgIds.get(pgId));
            }
        }
    }

    /**
     * Generate the list of program records seen in the SAM file, approximating this with those in the header that were not
     * themselves mentioned elsewhere.
     */
    private void setPGIdsSeen(final SAMFileHeader header) {
        final Set<String> pgIdsSeenAsPrevious = new HashSet<String>();

        // get all program record ids that are mentioned as previously seen
        for (final SAMProgramRecord samProgramRecord : header.getProgramRecords()) {
            final String previousProgramGroupID = samProgramRecord.getPreviousProgramGroupId();
            if (null != previousProgramGroupID) pgIdsSeenAsPrevious.add(previousProgramGroupID);
        }

        // ignore those that were previously seen
        for (final SAMProgramRecord samProgramRecord : header.getProgramRecords()) {
            final String pgId = samProgramRecord.getId();
            if (!pgIdsSeenAsPrevious.contains(pgId)) this.pgIdsSeen.add(pgId);
        }
    }
}
