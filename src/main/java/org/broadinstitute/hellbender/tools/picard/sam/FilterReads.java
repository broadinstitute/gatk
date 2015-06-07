package org.broadinstitute.hellbender.tools.picard.sam;

import htsjdk.samtools.*;
import htsjdk.samtools.filter.AlignedFilter;
import htsjdk.samtools.filter.FilteringIterator;
import htsjdk.samtools.filter.ReadNameFilter;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import org.broadinstitute.hellbender.cmdline.*;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;

/**
 * From a SAM or BAM file, produce a new SAM or BAM by filtering aligned reads or a list of read
 * names provided in a file (one readname per line)
 * <p/>
 * $Id$
 */
@CommandLineProgramProperties(
        usage = "Produces a new SAM or BAM file by including or excluding aligned reads " +
                "or a list of reads names supplied in the READ_LIST_FILE from the INPUT SAM or BAM file.\n",
        usageShort = "Creates a new SAM or BAM file by including or excluding aligned reads",
        programGroup = ReadProgramGroup.class
)
public final class FilterReads extends PicardCommandLineProgram {

    private static final Log log = Log.getInstance(FilterReads.class);

    private static enum Filter {
        includeAligned("OUTPUT SAM/BAM will contain aligned reads only. INPUT SAM/BAM must be in queryname SortOrder. (Note that *both* first and second of paired reads must be aligned to be included in the OUTPUT SAM or BAM)"),
        excludeAligned("OUTPUT SAM/BAM will contain un-mapped reads only. INPUT SAM/BAM must be in queryname SortOrder. (Note that *both* first and second of pair must be aligned to be excluded from the OUTPUT SAM or BAM)"),
        includeReadList("OUTPUT SAM/BAM will contain reads that are supplied in the READ_LIST_FILE file"),
        excludeReadList("OUTPUT bam will contain reads that are *not* supplied in the READ_LIST_FILE file");

        private final String description;

        Filter(final String description) {
            this.description = description;
        }

        @Override
        public String toString() {
            return this.name() + " [" + description + "]";
        }
    }

    @Argument(doc = "The SAM or BAM file that will be filtered.",
            optional = false,
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME)
    public File INPUT;

    @Argument(doc = "Filter.", optional = false)
    public Filter FILTER = null;

    @Argument(doc = "Read List File containing reads that will be included or excluded from the OUTPUT SAM or BAM file.",
            optional = true,
            shortName = "RLF")
    public File READ_LIST_FILE;

    @Argument(
            doc = "SortOrder of the OUTPUT SAM or BAM file, otherwise use the SortOrder of the INPUT file.",
            optional = true, shortName = "SO")
    public SAMFileHeader.SortOrder SORT_ORDER;

    @Argument(
            doc = "Create .reads files (for debugging purposes)",
            optional = true)
    public boolean WRITE_READS_FILES = true;

    @Argument(doc = "SAM or BAM file to write read excluded results to",
            optional = false, shortName = "O")
    public File OUTPUT;

    private void filterReads(final FilteringIterator filteringIterator) {

        // get OUTPUT header from INPUT and overwrite it if necessary
        final SAMFileHeader fileHeader = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).getFileHeader(INPUT);
        final SAMFileHeader.SortOrder inputSortOrder = fileHeader.getSortOrder();
        if (SORT_ORDER != null) {
            fileHeader.setSortOrder(SORT_ORDER);
        }
        final boolean presorted = inputSortOrder.equals(fileHeader.getSortOrder());
        log.info("Filtering [presorted=" + presorted + "] " + INPUT.getName() + " -> OUTPUT=" +
                OUTPUT.getName() + " [sortorder=" + fileHeader.getSortOrder().name() + "]");

        final ProgressLogger progress = new ProgressLogger(log, (int) 1e6, "Written");

        // create OUTPUT file
        try (final SAMFileWriter outputWriter = new SAMFileWriterFactory().makeSAMOrBAMWriter(fileHeader, presorted, OUTPUT)) {

           while (filteringIterator.hasNext()) {
               final SAMRecord rec = filteringIterator.next();
               outputWriter.addAlignment(rec);
               progress.record(rec);
           }

           filteringIterator.close();
       }
       log.info(new DecimalFormat("#,###").format(progress.getCount()) + " SAMRecords written to " + OUTPUT.getName());
    }

    /**
     * Write out a file of read names for debugging purposes.
     *
     * @param samOrBamFile The SAM or BAM file for which we are going to write out a file of its
     *                     containing read names
     */
    private void writeReadsFile(final File samOrBamFile) throws IOException {
        File readsFile = null;
        try (final SamReader reader = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(samOrBamFile)) {
            readsFile = new File(OUTPUT.getParentFile(), IOUtil.basename(samOrBamFile) + ".reads");
            IOUtil.assertFileIsWritable(readsFile);
            try (final BufferedWriter bw = IOUtil.openFileForBufferedWriting(readsFile, false)) {

                for (final SAMRecord rec : reader) {
                    bw.write(rec.toString() + "\n");
                }
            }
            IOUtil.assertFileIsReadable(readsFile);
        }
    }

    @Override
    protected Object doWork() {
        try {
            IOUtil.assertFileIsReadable(INPUT);
            IOUtil.assertFileIsWritable(OUTPUT);
            if (WRITE_READS_FILES) writeReadsFile(INPUT);

            switch (FILTER) {
                case includeAligned:
                    filterReads(new FilteringIterator(SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(INPUT).iterator(),
                            new AlignedFilter(true), true));
                    break;
                case excludeAligned:
                    filterReads(new FilteringIterator(SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(INPUT).iterator(),
                            new AlignedFilter(false), true));
                    break;
                case includeReadList:
                    filterReads(new FilteringIterator(SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(INPUT).iterator(),
                            new ReadNameFilter(READ_LIST_FILE, true)));
                    break;
                case excludeReadList:
                    filterReads(new FilteringIterator(SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(INPUT).iterator(),
                            new ReadNameFilter(READ_LIST_FILE, false)));
                    break;
                default:
                    throw new UnsupportedOperationException(FILTER.name() + " has not been implemented!");
            }

            IOUtil.assertFileIsReadable(OUTPUT);
            if (WRITE_READS_FILES) writeReadsFile(OUTPUT);

        } catch (IOException e) {
            if (OUTPUT.exists() && !OUTPUT.delete()) {
                throw new UserException("Failed to delete existing output: " + OUTPUT.getAbsolutePath());
            } else {
                throw new UserException("Failed to filter " + INPUT.getName());
            }
        }

        return null;
    }

    @Override
    protected String[] customCommandLineValidation() {
        if (INPUT.equals(OUTPUT)) {
            return new String[]{"INPUT file and OUTPUT file must differ!"};
        }

        if ((FILTER.equals(Filter.includeReadList) ||
                FILTER.equals(Filter.excludeReadList)) &&
                READ_LIST_FILE == null) {
            return new String[]{"A READ_LIST_FILE must be specified when using the " + FILTER.name() + " option"};

        }

        return super.customCommandLineValidation();
    }
}
