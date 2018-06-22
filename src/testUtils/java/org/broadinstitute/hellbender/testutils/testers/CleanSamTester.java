package org.broadinstitute.hellbender.testutils.testers;

import htsjdk.samtools.*;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.TestUtil;
import picard.sam.CleanSam;
import org.testng.Assert;

import java.io.PrintWriter;
import java.util.*;

/**
 * This class is the extension of the SamFileTester to test CleanSam with SAM files generated on the fly.
 */
// TODO: this should live in the tests sources, because is testing just one class
// TODO: and it is not really part of the testing framework (which should be helpers for re-use with different classes)
public final class CleanSamTester extends SamFileTester {
    private final String expectedCigar;

    @Override
    public String getTestedToolName() { return CleanSam.class.getSimpleName(); }

    public CleanSamTester(final String expectedCigar, final int readLength, final int defaultChromosomeLength) {
        super(readLength, true, defaultChromosomeLength);
        this.expectedCigar = expectedCigar;
    }


    @Override
    protected void test() {
        try {
            final SamFileValidator validator = new SamFileValidator(new PrintWriter(System.out), 8000);

            // Validate it has the expected cigar
            validator.setIgnoreWarnings(true);
            validator.setVerbose(true, 1000);
            validator.setErrorsToIgnore(Collections.singletonList(SAMValidationError.Type.MISSING_READ_GROUP));
            SamReaderFactory factory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.LENIENT);
            SamReader samReader = factory.open(getOutput());
            final SAMRecordIterator iterator = samReader.iterator();
            while (iterator.hasNext()) {
                final SAMRecord rec = iterator.next();
                Assert.assertEquals(rec.getCigarString(), expectedCigar);
                if (SAMUtils.hasMateCigar(rec)) {
                    Assert.assertEquals(SAMUtils.getMateCigarString(rec), expectedCigar);
                }
            }
            CloserUtil.close(samReader);

            // Run validation on the output file
            samReader = factory.open(getOutput());
            final boolean validated = validator.validateSamFileVerbose(samReader, null);
            CloserUtil.close(samReader);

            Assert.assertTrue(validated, "ValidateSamFile failed");
        } finally {
            TestUtil.recursiveDelete(getOutputDir());
        }
    }
}
