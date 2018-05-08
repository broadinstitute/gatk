package org.broadinstitute.hellbender.utils.test.testers;

import htsjdk.samtools.*;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.utils.test.CommandLineProgramTester;

import java.io.File;
import java.util.*;

/**
 * Abstract class for doing basic on the fly SAM file testing.
 */
public abstract class SamFileTester implements CommandLineProgramTester {

    private final SAMRecordSetBuilder samRecordSetBuilder;
    protected final Map<String, Boolean> duplicateFlags = new LinkedHashMap<>();
    private File outputDir;
    private File output;
    private int readNameCounter = 0;
    private boolean noMateCigars = false;
    private boolean deleteOnExit = true;
    protected final List<String> args = new ArrayList<>();

    public SamFileTester(final int readLength, final boolean deleteOnExit, final int defaultChromosomeLength, final DuplicateScoringStrategy.ScoringStrategy duplicateScoringStrategy) {
        this.deleteOnExit = deleteOnExit;
        this.samRecordSetBuilder = new SAMRecordSetBuilder(true, SAMFileHeader.SortOrder.coordinate, true, defaultChromosomeLength, duplicateScoringStrategy);
        samRecordSetBuilder.setReadLength(readLength);
        setOutputDir();
    }

    public SamFileTester(final int readLength, final boolean deleteOnExit, final int defaultChromosomeLength) {
        this(readLength, deleteOnExit, defaultChromosomeLength, SAMRecordSetBuilder.DEFAULT_DUPLICATE_SCORING_STRATEGY);
    }

    public void setHeader(final SAMFileHeader header) {
        this.samRecordSetBuilder.setHeader(header);
    }

    public SAMFileHeader getHeader() {
        return this.samRecordSetBuilder.getHeader();
    }

    public void addRecord(final SAMRecord record) {
        this.duplicateFlags.put(samRecordToDuplicatesFlagsKey(record), record.getDuplicateReadFlag());
        this.samRecordSetBuilder.addRecord(record);
    }

    public int getNumberOfRecords() {
        return this.samRecordSetBuilder.size();
    }

    public CloseableIterator<SAMRecord> getRecordIterator() {
        return this.samRecordSetBuilder.iterator();
    }

    public File getOutput() {
        return output;
    }

    public void setOutput(final File output) {
        this.output = output;
    }

    public void addArg(final String key) {
        args.add(key);
    }

    public void addArg(final String key, final String value) {
        args.add(key); args.add(value);
    }

    public List<String> getArgs() {
        return args;
    }

    public File getOutputDir() {
        return outputDir;
    }

    private void setOutputDir() {
        this.outputDir = IOUtil.createTempDir(this.getClass().getSimpleName() + ".", ".tmp");
        if (deleteOnExit) {
          outputDir.deleteOnExit();
        }
    }

    public void setNoMateCigars(final boolean value) {
        this.noMateCigars = value;
    }

    public boolean getDeleteOnExit() {
        return deleteOnExit;
    }

    protected String samRecordToDuplicatesFlagsKey(final SAMRecord record) {
        String readName = record.getReadName()
                + "-"
                + record.getReadPairedFlag()
                + "-";
        if (record.getReadPairedFlag()) {
            readName += record.getFirstOfPairFlag()
                    + "-"
                    + record.getSecondOfPairFlag();
        } else {
            readName += "false-false";
        }
        return readName;
    }

    // Below are a bunch of utility methods for adding records to the SAMRecordSetBuilder
    public void addUnmappedFragment(final int referenceSequenceIndex,
                                    final int defaultQualityScore) {
        addFragment(referenceSequenceIndex, -1, true, false, null, null, defaultQualityScore, false);
    }

    public void addUnmappedPair(final int referenceSequenceIndex,
                                final int defaultQualityScore) {
        addMatePair(referenceSequenceIndex, -1, -1, true, true, false, false, null, null, false, false, false, false, false, defaultQualityScore);
    }

    public void addMappedFragment(final int referenceSequenceIndex, final int alignmentStart, final boolean isDuplicate,
                                  final int defaultQualityScore) {
        addFragment(referenceSequenceIndex, alignmentStart, false, isDuplicate, null, null, defaultQualityScore, false);
    }

    public void addMappedFragment(final int referenceSequenceIndex, final int alignmentStart, final boolean isDuplicate,
                                  final int defaultQualityScore, final boolean isSecondary) {
        addFragment(referenceSequenceIndex, alignmentStart, false, isDuplicate, null, null, defaultQualityScore, isSecondary);
    }

    public void addMappedFragment(final int referenceSequenceIndex, final int alignmentStart, final boolean isDuplicate, final String cigar,
                                  final int defaultQualityScore) {
        addFragment(referenceSequenceIndex, alignmentStart, false, isDuplicate, cigar, null, defaultQualityScore, false);
    }

    public void addMappedPair(final int referenceSequenceIndex,
                              final int alignmentStart1,
                              final int alignmentStart2,
                              final boolean isDuplicate1,
                              final boolean isDuplicate2,
                              final int defaultQualityScore) {
        addMappedPair(referenceSequenceIndex, alignmentStart1, alignmentStart2, isDuplicate1, isDuplicate2, null, null,
                false, defaultQualityScore);
    }

    public void addMappedPair(final int referenceSequenceIndex,
                              final int alignmentStart1,
                              final int alignmentStart2,
                              final boolean isDuplicate1,
                              final boolean isDuplicate2,
                              final String cigar1,
                              final String cigar2,
                              final boolean firstOnly,
                              final int defaultQualityScore) {
        addMappedPair(referenceSequenceIndex, alignmentStart1, alignmentStart2, isDuplicate1, isDuplicate2, cigar1,
                cigar2, false, true, firstOnly, defaultQualityScore);
    }

    public void addMappedPair(final int referenceSequenceIndex,
                              final int alignmentStart1,
                              final int alignmentStart2,
                              final boolean isDuplicate1,
                              final boolean isDuplicate2,
                              final String cigar1,
                              final String cigar2,
                              final boolean strand1,
                              final boolean strand2,
                              final boolean firstOnly,
                              final int defaultQualityScore) {
        addMatePair(referenceSequenceIndex, alignmentStart1, alignmentStart2, false, false, isDuplicate1, isDuplicate2, cigar1, cigar2,
                strand1, strand2, firstOnly, false, false, defaultQualityScore);
    }

    public void addMatePair(final int referenceSequenceIndex,
                            final int alignmentStart1,
                            final int alignmentStart2,
                            final boolean record1Unmapped,
                            final boolean record2Unmapped,
                            final boolean isDuplicate1,
                            final boolean isDuplicate2,
                            final String cigar1,
                            final String cigar2,
                            final boolean strand1,
                            final boolean strand2,
                            final boolean firstOnly,
                            final boolean record1NonPrimary,
                            final boolean record2NonPrimary,
                            final int defaultQualityScore) {
        addMatePair("READ" + readNameCounter++, referenceSequenceIndex, alignmentStart1, alignmentStart2, record1Unmapped, record2Unmapped,
                isDuplicate1, isDuplicate2, cigar1, cigar2, strand1, strand2, firstOnly, record1NonPrimary, record2NonPrimary,
                defaultQualityScore);
    }

    private void addFragment(final int referenceSequenceIndex, final int alignmentStart, final boolean recordUnmapped, final boolean isDuplicate, final String cigar,
                             final String qualityString, final int defaultQualityScore, final boolean isSecondary) {
        final SAMRecord record = samRecordSetBuilder.addFrag("READ" + readNameCounter++, referenceSequenceIndex, alignmentStart, false,
                recordUnmapped, cigar, qualityString, defaultQualityScore, isSecondary);

        this.duplicateFlags.put(samRecordToDuplicatesFlagsKey(record), isDuplicate);
    }

    public void addMatePair(final String readName,
                            final int referenceSequenceIndex1,
                            final int referenceSequenceIndex2,
                            final int alignmentStart1,
                            final int alignmentStart2,
                            final boolean record1Unmapped,
                            final boolean record2Unmapped,
                            final boolean isDuplicate1,
                            final boolean isDuplicate2,
                            final String cigar1,
                            final String cigar2,
                            final boolean strand1,
                            final boolean strand2,
                            final boolean firstOnly,
                            final boolean record1NonPrimary,
                            final boolean record2NonPrimary,
                            final int defaultQuality,
                            final String readGorup) {
        final List<SAMRecord> samRecordList = samRecordSetBuilder.addPair(readName, referenceSequenceIndex1, referenceSequenceIndex2, alignmentStart1, alignmentStart2,
                record1Unmapped, record2Unmapped, cigar1, cigar2, strand1, strand2, record1NonPrimary, record2NonPrimary, defaultQuality);

        final SAMRecord record1 = samRecordList.get(0);
        final SAMRecord record2 = samRecordList.get(1);

        if (this.noMateCigars) {
            record1.setAttribute("MC", null);
            record2.setAttribute("MC", null);
        }

        if (readGorup!=null) {
            record1.setAttribute("RG", readGorup);
            record2.setAttribute("RG", readGorup);
        }

        if (firstOnly) {
            samRecordSetBuilder.getRecords().remove(record2);
        }

        this.duplicateFlags.put(samRecordToDuplicatesFlagsKey(record1), isDuplicate1);
        this.duplicateFlags.put(samRecordToDuplicatesFlagsKey(record2), isDuplicate2);
    }

    public void addMatePair(final String readName,
                            final int referenceSequenceIndex,
                            final int alignmentStart1,
                            final int alignmentStart2,
                            final boolean record1Unmapped,
                            final boolean record2Unmapped,
                            final boolean isDuplicate1,
                            final boolean isDuplicate2,
                            final String cigar1,
                            final String cigar2,
                            final boolean strand1,
                            final boolean strand2,
                            final boolean firstOnly,
                            final boolean record1NonPrimary,
                            final boolean record2NonPrimary,
                            final int defaultQuality) {
        addMatePair(readName, referenceSequenceIndex,referenceSequenceIndex, alignmentStart1, alignmentStart2, record1Unmapped, record2Unmapped,
                isDuplicate1, isDuplicate2, cigar1, cigar2, strand1, strand2, firstOnly, record1NonPrimary, record2NonPrimary, defaultQuality, null);
    }

    protected abstract void test();

    /**
     * Sets up the basic command line arguments for input and output and runs instanceMain.
     */
    public void runTest() {
        final File input = createInputFile();
        output = new File(outputDir, "output.sam");
        addArg("--" + StandardArgumentDefinitions.INPUT_LONG_NAME, input.getAbsolutePath());
        addArg("--" + StandardArgumentDefinitions.OUTPUT_LONG_NAME, output.getAbsolutePath());
        addArgs();
        runCommandLine(args);
        test();
    }

    protected void addArgs() {
        // subclasses may override to add more arguments
    }

    private File createInputFile() {
        // Create the input file
        final File input = new File(outputDir, "input.bam");
        try (final SAMFileWriter writer = new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(samRecordSetBuilder.getHeader(), true, input)) {
            for (final SAMRecord record : samRecordSetBuilder.getRecords()) {
                writer.addAlignment(record);
            }
        }
        return input;
    }

    public SamReader getInput() {
        return samRecordSetBuilder.getSamReader();
    }
}
