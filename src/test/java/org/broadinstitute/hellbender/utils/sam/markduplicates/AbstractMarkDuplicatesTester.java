/*
 * The MIT License
 *
 * Copyright (c) 2014 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package org.broadinstitute.hellbender.utils.sam.markduplicates;

import htsjdk.samtools.DuplicateScoringStrategy.ScoringStrategy;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordSetBuilder;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.FormatUtil;
import htsjdk.samtools.util.TestUtil;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.metrics.DuplicationMetrics;
import org.broadinstitute.hellbender.utils.sam.testers.SamFileTester;
import org.testng.Assert;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;

/**
 * This class is an extension of SamFileTester used to test AbstractMarkDuplicatesCommandLineProgram's with SAM files generated on the fly.
 * This performs the underlying tests defined by classes such as AbstractMarkDuplicatesCommandLineProgramTest.
 */
abstract public class AbstractMarkDuplicatesTester extends SamFileTester {

    final private File metricsFile;
    final DuplicationMetrics expectedMetrics;

    public AbstractMarkDuplicatesTester(final ScoringStrategy duplicateScoringStrategy) {
        super(50, true, SAMRecordSetBuilder.DEFAULT_CHROMOSOME_LENGTH, duplicateScoringStrategy);

        expectedMetrics = new DuplicationMetrics();
        expectedMetrics.READ_PAIR_OPTICAL_DUPLICATES = 0;

        metricsFile = new File(getOutputDir(), "metrics.txt");
        addArg("METRICS_FILE=" + metricsFile);
        addArg("DUPLICATE_SCORING_STRATEGY=" + duplicateScoringStrategy.name());
    }


    public AbstractMarkDuplicatesTester() {
        this(SAMRecordSetBuilder.DEFAULT_DUPLICATE_SCORING_STRATEGY);
    }

    @Override
    public String getCommandLineProgramName() { return getProgram().getClass().getSimpleName(); }

    /**
     * Fill in expected duplication metrics directly from the input records given to this tester
     */
    private void updateExpectedDuplicationMetrics() {

        final FormatUtil formatter = new FormatUtil();

        final CloseableIterator<SAMRecord> inputRecordIterator = this.getRecordIterator();
        while (inputRecordIterator.hasNext()) {
            final SAMRecord record = inputRecordIterator.next();
            if (!record.isSecondaryOrSupplementary()) {
                final String key = samRecordToDuplicatesFlagsKey(record);
                if (!this.duplicateFlags.containsKey(key)) {
                    System.err.println("DOES NOT CONTAIN KEY: " + key);
                }
                final boolean isDuplicate = this.duplicateFlags.get(key);

                // First bring the simple metricsFile up to date
                if (record.getReadUnmappedFlag()) {
                    ++expectedMetrics.UNMAPPED_READS;
                }
                else if (!record.getReadPairedFlag() || record.getMateUnmappedFlag()) {
                    ++expectedMetrics.UNPAIRED_READS_EXAMINED;
                    if (isDuplicate) ++expectedMetrics.UNPAIRED_READ_DUPLICATES;
                }
                else {
                    ++expectedMetrics.READ_PAIRS_EXAMINED; // will need to be divided by 2 at the end
                    if (isDuplicate) ++expectedMetrics.READ_PAIR_DUPLICATES; // will need to be divided by 2 at the end
                }
            }
        }
        expectedMetrics.READ_PAIR_DUPLICATES = expectedMetrics.READ_PAIR_DUPLICATES / 2;
        expectedMetrics.READ_PAIRS_EXAMINED = expectedMetrics.READ_PAIRS_EXAMINED / 2;
        expectedMetrics.calculateDerivedMetrics();

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
            final SamReader reader = SamReaderFactory.makeDefault().open(getOutput());
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
            final MetricsFile<DuplicationMetrics, Comparable<?>> metricsOutput = new MetricsFile<DuplicationMetrics, Comparable<?>>();
            try{
                metricsOutput.read(new FileReader(metricsFile));
            }
            catch (final FileNotFoundException ex) {
                System.err.println("Metrics file not found: " + ex);
            }
            // NB: Test writes an initial metrics line with a null entry for LIBRARY and 0 values for all metrics. Why?
            final DuplicationMetrics observedMetrics = metricsOutput.getMetrics().get(metricsOutput.getMetrics().size() - 1);
            Assert.assertEquals(observedMetrics.UNPAIRED_READS_EXAMINED, expectedMetrics.UNPAIRED_READS_EXAMINED, "UNPAIRED_READS_EXAMINED does not match expected");
            Assert.assertEquals(observedMetrics.READ_PAIRS_EXAMINED, expectedMetrics.READ_PAIRS_EXAMINED, "READ_PAIRS_EXAMINED does not match expected");
            Assert.assertEquals(observedMetrics.UNMAPPED_READS, expectedMetrics.UNMAPPED_READS, "UNMAPPED_READS does not match expected");
            Assert.assertEquals(observedMetrics.UNPAIRED_READ_DUPLICATES, expectedMetrics.UNPAIRED_READ_DUPLICATES, "UNPAIRED_READ_DUPLICATES does not match expected");
            Assert.assertEquals(observedMetrics.READ_PAIR_DUPLICATES, expectedMetrics.READ_PAIR_DUPLICATES, "READ_PAIR_DUPLICATES does not match expected");
            Assert.assertEquals(observedMetrics.READ_PAIR_OPTICAL_DUPLICATES, expectedMetrics.READ_PAIR_OPTICAL_DUPLICATES, "READ_PAIR_OPTICAL_DUPLICATES does not match expected");
            Assert.assertEquals(observedMetrics.PERCENT_DUPLICATION, expectedMetrics.PERCENT_DUPLICATION, "PERCENT_DUPLICATION does not match expected");
            Assert.assertEquals(observedMetrics.ESTIMATED_LIBRARY_SIZE, expectedMetrics.ESTIMATED_LIBRARY_SIZE, "ESTIMATED_LIBRARY_SIZE does not match expected");
        } finally {
            TestUtil.recursiveDelete(getOutputDir());
        }
    }

    abstract protected CommandLineProgram getProgram();
}

