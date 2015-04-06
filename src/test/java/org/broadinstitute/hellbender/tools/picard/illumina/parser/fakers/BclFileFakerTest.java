package org.broadinstitute.hellbender.tools.picard.illumina.parser.fakers;

import org.broadinstitute.hellbender.tools.picard.illumina.parser.readers.BclQualityEvaluationStrategy;
import org.broadinstitute.hellbender.tools.picard.illumina.parser.readers.BclReader;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;

public class BclFileFakerTest {

    @Test
    public void testFileLengthMatchesHeaderLength() throws Exception {
        final File fakeFile = File.createTempFile("BclFileFakerTest", ".bcl");
        fakeFile.deleteOnExit();

        new BclFileFaker().fakeFile(fakeFile, 100000);
        // .make() has a number of checks for the file
        final BclReader bclReader = new BclReader(
            fakeFile,
            new BclQualityEvaluationStrategy(BclQualityEvaluationStrategy.ILLUMINA_ALLEGED_MINIMUM_QUALITY), false);
        Assert.assertEquals(100000, BclReader.getNumberOfClusters(fakeFile));
        Assert.assertEquals(BclReader.getNumberOfClusters(fakeFile), fakeFile.length() - 4);
    }

    @Test
    public void testGZFileIsActuallyGZipped() throws Exception {
        final File fakeFile = File.createTempFile("BclFileFakerTest", ".bcl.gz");
        fakeFile.deleteOnExit();

        new BclFileFaker().fakeFile(fakeFile, 100000);
        new BclReader(
                fakeFile,
                new BclQualityEvaluationStrategy(BclQualityEvaluationStrategy.ILLUMINA_ALLEGED_MINIMUM_QUALITY), false);
    }
}
