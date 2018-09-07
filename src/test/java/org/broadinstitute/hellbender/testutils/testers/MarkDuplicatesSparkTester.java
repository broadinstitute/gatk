package org.broadinstitute.hellbender.testutils.testers;

import htsjdk.samtools.*;
import htsjdk.samtools.DuplicateScoringStrategy.ScoringStrategy;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.FormatUtil;
import htsjdk.samtools.util.TestUtil;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.MarkDuplicatesSparkArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.transforms.markduplicates.MarkDuplicatesSpark;
import org.testng.Assert;
import picard.sam.DuplicationMetrics;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;

/**
 * This class is an extension of SamFileTester used to test AbstractMarkDuplicatesCommandLineProgram's with SAM files generated on the fly.
 * This performs the underlying tests defined by classes such as AbstractMarkDuplicatesCommandLineProgramTest.
 */
public class MarkDuplicatesSparkTester extends SamFileTester {

    final private File metricsFile;
    final DuplicationMetrics expectedMetrics;
    boolean markUnmappedReads = false;

    public MarkDuplicatesSparkTester(final ScoringStrategy duplicateScoringStrategy, final SAMFileHeader.SortOrder sortOrder) {
        this(duplicateScoringStrategy, sortOrder, true);
    }

    public MarkDuplicatesSparkTester(final ScoringStrategy duplicateScoringStrategy, final SAMFileHeader.SortOrder sortOrder, final boolean recordNeedSorting) {
        super(50, true, SAMRecordSetBuilder.DEFAULT_CHROMOSOME_LENGTH, duplicateScoringStrategy, sortOrder, recordNeedSorting);

        expectedMetrics = new DuplicationMetrics();
        expectedMetrics.READ_PAIR_OPTICAL_DUPLICATES = 0;

        metricsFile = new File(getOutputDir(), "metrics.txt");
        addArg("--"+StandardArgumentDefinitions.METRICS_FILE_LONG_NAME, metricsFile.getAbsolutePath());
        addArg("--"+StandardArgumentDefinitions.DUPLICATE_SCORING_STRATEGY_LONG_NAME, duplicateScoringStrategy.name());
    }

    public MarkDuplicatesSparkTester(final ScoringStrategy duplicateScoringStrategy, boolean markUnmappedReads) {
        this(duplicateScoringStrategy, SAMFileHeader.SortOrder.coordinate);
        this.markUnmappedReads = markUnmappedReads;
    }

    public MarkDuplicatesSparkTester(boolean markUnmappedReads) {
        this(DuplicateScoringStrategy.ScoringStrategy.TOTAL_MAPPED_REFERENCE_LENGTH, markUnmappedReads);
    }

    public MarkDuplicatesSparkTester() {
        this(false);
    }

    @Override
    public final String getTestedToolName() { return getProgram().getClass().getSimpleName(); }

    /**
     * Fill in expected duplication metrics directly from the input records given to this tester
     */
    private void updateExpectedDuplicationMetrics() {

        final FormatUtil formatter = new FormatUtil();

        final CloseableIterator<SAMRecord> inputRecordIterator = this.getRecordIterator();
        while (inputRecordIterator.hasNext()) {
            final SAMRecord record = inputRecordIterator.next();
            if (record.isSecondaryOrSupplementary()) {
                ++expectedMetrics.SECONDARY_OR_SUPPLEMENTARY_RDS;
            } else {
                final String key = samRecordToDuplicatesFlagsKey(record);
                if (!this.duplicateFlags.containsKey(key)) {
                    System.err.println("DOES NOT CONTAIN KEY: " + key);
                }
                final boolean isDuplicate = this.duplicateFlags.get(key);

                // First bring the simple metricsFile up to date
                if (record.getReadUnmappedFlag()) {
                    ++expectedMetrics.UNMAPPED_READS;
                } else if (!record.getReadPairedFlag() || record.getMateUnmappedFlag()) {
                    ++expectedMetrics.UNPAIRED_READS_EXAMINED;
                    if (isDuplicate) ++expectedMetrics.UNPAIRED_READ_DUPLICATES;
                } else {
                    ++expectedMetrics.READ_PAIRS_EXAMINED; // will need to be divided by 2 at the end
                    if (isDuplicate) ++expectedMetrics.READ_PAIR_DUPLICATES; // will need to be divided by 2 at the end
                }
            }
        }
        expectedMetrics.READ_PAIR_DUPLICATES = expectedMetrics.READ_PAIR_DUPLICATES / 2;
        expectedMetrics.READ_PAIRS_EXAMINED = expectedMetrics.READ_PAIRS_EXAMINED / 2;
        expectedMetrics.calculateDerivedFields();

        // Have to run this Double value through the same format/parsing operations as during a file write/read
        expectedMetrics.PERCENT_DUPLICATION = formatter.parseDouble(formatter.format(expectedMetrics.PERCENT_DUPLICATION));
    }

    public void setExpectedOpticalDuplicate(final int expectedOpticalDuplicatePairs) {
        expectedMetrics.READ_PAIR_OPTICAL_DUPLICATES = expectedOpticalDuplicatePairs;
    }

    @Override
    public void test() {
        try {
            updateExpectedDuplicationMetrics();
            // Read the output and check the duplicate flag
            int outputRecords = 0;
            final SamReader reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.LENIENT).open(getOutput());
            for (final SAMRecord record : reader) {
                outputRecords++;
                final String key = samRecordToDuplicatesFlagsKey(record);
                if (!this.duplicateFlags.containsKey(key)) {
                    System.err.println("DOES NOT CONTAIN KEY: " + key);
                }
                Assert.assertTrue(this.duplicateFlags.containsKey(key));
                final boolean value = this.duplicateFlags.get(key);
                this.duplicateFlags.remove(key);
                if (value != record.getDuplicateReadFlag()) {
                    System.err.println("Mismatching read:");
                    System.err.print(record.getSAMString());
                }
                Assert.assertEquals(record.getDuplicateReadFlag(), value);
            }
            CloserUtil.close(reader);

            // Ensure the program output the same number of records as were read in
            Assert.assertEquals(outputRecords, this.getNumberOfRecords(), ("saw " + outputRecords + " output records, vs. " + this.getNumberOfRecords() + " input records"));

            // Check the values written to metrics.txt against our input expectations
            final MetricsFile<DuplicationMetrics, Comparable<?>> metricsOutput = new MetricsFile<>();
            try{
                metricsOutput.read(new FileReader(metricsFile));
            }
            catch (final FileNotFoundException ex) {
                System.err.println("Metrics file not found: " + ex);
            }

            Assert.assertEquals(metricsOutput.getMetrics().size(), 1);
            final DuplicationMetrics observedMetrics = metricsOutput.getMetrics().get(0);
            Assert.assertEquals(observedMetrics.UNPAIRED_READS_EXAMINED, expectedMetrics.UNPAIRED_READS_EXAMINED, "UNPAIRED_READS_EXAMINED does not match expected");
            Assert.assertEquals(observedMetrics.READ_PAIRS_EXAMINED, expectedMetrics.READ_PAIRS_EXAMINED, "READ_PAIRS_EXAMINED does not match expected");
            Assert.assertEquals(observedMetrics.UNMAPPED_READS, expectedMetrics.UNMAPPED_READS, "UNMAPPED_READS does not match expected");
            Assert.assertEquals(observedMetrics.UNPAIRED_READ_DUPLICATES, expectedMetrics.UNPAIRED_READ_DUPLICATES, "UNPAIRED_READ_DUPLICATES does not match expected");
            Assert.assertEquals(observedMetrics.READ_PAIR_DUPLICATES, expectedMetrics.READ_PAIR_DUPLICATES, "READ_PAIR_DUPLICATES does not match expected");
            Assert.assertEquals(observedMetrics.READ_PAIR_OPTICAL_DUPLICATES, expectedMetrics.READ_PAIR_OPTICAL_DUPLICATES, "READ_PAIR_OPTICAL_DUPLICATES does not match expected");
            Assert.assertEquals(observedMetrics.PERCENT_DUPLICATION, expectedMetrics.PERCENT_DUPLICATION, "PERCENT_DUPLICATION does not match expected");
            // The dataflow version outputs 0 instead of null because it is
            // coded and needs to have real values for each field.
            if (observedMetrics.ESTIMATED_LIBRARY_SIZE == null) {
              observedMetrics.ESTIMATED_LIBRARY_SIZE = 0L;
            }
            if (expectedMetrics.ESTIMATED_LIBRARY_SIZE == null) {
              expectedMetrics.ESTIMATED_LIBRARY_SIZE = 0L;
            }
            Assert.assertEquals(observedMetrics.ESTIMATED_LIBRARY_SIZE, expectedMetrics.ESTIMATED_LIBRARY_SIZE, "ESTIMATED_LIBRARY_SIZE does not match expected");
            Assert.assertEquals(observedMetrics.SECONDARY_OR_SUPPLEMENTARY_RDS, expectedMetrics.SECONDARY_OR_SUPPLEMENTARY_RDS, "SECONDARY_OR_SUPPLEMENTARY_RDS does not match expected");
        } finally {
            TestUtil.recursiveDelete(getOutputDir());
        }
    }

    protected CommandLineProgram getProgram() { return new MarkDuplicatesSpark(); }

    @Override
    protected void addArgs() {
        if (!markUnmappedReads) {
            addArg("--" + MarkDuplicatesSparkArgumentCollection.DO_NOT_MARK_UNMAPPED_MATES_LONG_NAME);
        }
    }
}
