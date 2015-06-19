package org.broadinstitute.hellbender.tools.picard.sam;

import htsjdk.samtools.*;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import org.broadinstitute.hellbender.cmdline.*;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

/**
 * Reads a SAM or BAM file and combines the output to one file
 *
 * @author Tim Fennell
 */
@CommandLineProgramProperties(
        usage = "Merges multiple SAM/BAM files into one file.",
        usageShort = "Merges multiple SAM or BAM files into one file",
        programGroup = ReadProgramGroup.class
)
public final class MergeSamFiles extends PicardCommandLineProgram {
    private static final Log log = Log.getInstance(MergeSamFiles.class);

    @Argument(shortName = "I", doc = "SAM or BAM input file", optional=false)
    public List<File> INPUT = new ArrayList<>();

    @Argument(shortName = "O", doc = "SAM or BAM file to write merged result to")
    public File OUTPUT;

    @Argument(shortName = StandardArgumentDefinitions.SORT_ORDER_SHORT_NAME, doc = "Sort order of output file", optional = true)
    public SAMFileHeader.SortOrder SORT_ORDER = SAMFileHeader.SortOrder.coordinate;

    @Argument(doc = "If true, assume that the input files are in the same sort order as the requested output sort order, even if their headers say otherwise.",
            shortName = StandardArgumentDefinitions.ASSUME_SORTED_SHORT_NAME)
    public boolean ASSUME_SORTED = false;

    @Argument(shortName = "MSD", doc = "Merge the sequence dictionaries", optional = true)
    public boolean MERGE_SEQUENCE_DICTIONARIES = false;

    @Argument(doc = "Option to create a background thread to encode, " +
            "compress and write to disk the output file. The threaded version uses about 20% more CPU and decreases " +
            "runtime by ~20% when writing out a compressed BAM file.")
    public boolean USE_THREADING = false;

    @Argument(doc = "Comment(s) to include in the merged output file's header.", optional = true, shortName = "CO")
    public List<String> COMMENT = new ArrayList<>();

    private static final int PROGRESS_INTERVAL = 1000000;

    /** Combines multiple SAM/BAM files into one. */
    @Override
    protected Object doWork() {
        boolean matchedSortOrders = true;

        // Open the files for reading and writing
        final List<SamReader> readers = new ArrayList<>();
        final List<SAMFileHeader> headers = new ArrayList<>();
        {
            SAMSequenceDictionary dict = null; // Used to try and reduce redundant SDs in memory

            for (final File inFile : INPUT) {
                IOUtil.assertFileIsReadable(inFile);
                final SamReader in = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(inFile);
                readers.add(in);
                headers.add(in.getFileHeader());

                // A slightly hackish attempt to keep memory consumption down when merging multiple files with
                // large sequence dictionaries (10,000s of sequences). If the dictionaries are identical, then
                // replace the duplicate copies with a single dictionary to reduce the memory footprint.
                if (dict == null) {
                    dict = in.getFileHeader().getSequenceDictionary();
                } else if (dict.equals(in.getFileHeader().getSequenceDictionary())) {
                    in.getFileHeader().setSequenceDictionary(dict);
                }

                matchedSortOrders = matchedSortOrders && in.getFileHeader().getSortOrder() == SORT_ORDER;
            }
        }

        // If all the input sort orders match the output sort order then just merge them and
        // write on the fly, otherwise setup to merge and sort before writing out the final file
        IOUtil.assertFileIsWritable(OUTPUT);
        final boolean presorted;
        final SAMFileHeader.SortOrder headerMergerSortOrder;
        final boolean mergingSamRecordIteratorAssumeSorted;

        if (matchedSortOrders || SORT_ORDER == SAMFileHeader.SortOrder.unsorted || ASSUME_SORTED) {
            log.info("Input files are in same order as output so sorting to temp directory is not needed.");
            headerMergerSortOrder = SORT_ORDER;
            mergingSamRecordIteratorAssumeSorted = ASSUME_SORTED;
            presorted = true;
        } else {
            log.info("Sorting input files using temp directory " + TMP_DIR);
            headerMergerSortOrder = SAMFileHeader.SortOrder.unsorted;
            mergingSamRecordIteratorAssumeSorted = false;
            presorted = false;
        }
        final SamFileHeaderMerger headerMerger = new SamFileHeaderMerger(headerMergerSortOrder, headers, MERGE_SEQUENCE_DICTIONARIES);
        final MergingSamRecordIterator iterator = new MergingSamRecordIterator(headerMerger, readers, mergingSamRecordIteratorAssumeSorted);
        final SAMFileHeader header = headerMerger.getMergedHeader();
        for (final String comment : COMMENT) {
            header.addComment(comment);
        }
        header.setSortOrder(SORT_ORDER);
        final SAMFileWriterFactory samFileWriterFactory = new SAMFileWriterFactory();
        if (USE_THREADING) {
            samFileWriterFactory.setUseAsyncIo(true);
        }
        try (final SAMFileWriter out = samFileWriterFactory.makeSAMOrBAMWriter(header, presorted, OUTPUT)) {

            // Lastly loop through and write out the records
            final ProgressLogger progress = new ProgressLogger(log, PROGRESS_INTERVAL);
            while (iterator.hasNext()) {
                final SAMRecord record = iterator.next();
                out.addAlignment(record);
                progress.record(record);
            }

            log.info("Finished reading inputs.");
            CloserUtil.close(readers);
        }
        return null;
    }

    @Override
    protected String[] customCommandLineValidation() {
        if (CREATE_INDEX && SORT_ORDER != SAMFileHeader.SortOrder.coordinate) {
            return new String[]{"Can't CREATE_INDEX unless SORT_ORDER is coordinate"};
        }
        return null;
    }

}
