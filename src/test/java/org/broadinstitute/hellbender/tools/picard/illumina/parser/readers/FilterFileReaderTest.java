package org.broadinstitute.hellbender.tools.picard.illumina.parser.readers;

import htsjdk.samtools.SAMException;
import org.broadinstitute.hellbender.tools.picard.illumina.parser.fakers.FilterFileFaker;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.*;

public final class FilterFileReaderTest {
    public static File TEST_DATA_DIR = new File("src/test/resources/org/broadinstitute/hellbender/tools/picard/illumina/readerTests");
    public static final File PASSING_FILTER_FILE = new File(TEST_DATA_DIR, "pf_passing.filter");

    public static final boolean [] expectedPfs = {
            false, false, false, false,     true,  false, false, false,      true,  true,  false, false,   true,  false, true,  false,
            true,  false, false, true,      true,  true,  true,  false,      true,  true,  false, true,    true,  true,  true,  false,
            true,  true,  true,  true,      false, true,  false, false,      false, true,  true,  false,   false, true,  false, true,
            false, true,  true,  true,      false, false, true,  false,      false, false, true,  true,    false, false, false, true
    };

    @Test
    public void readValidFile() {
        final FilterFileReader reader = new FilterFileReader(PASSING_FILTER_FILE);
        Assert.assertEquals(reader.numClusters, expectedPfs.length);
        for(int i = 0; i < reader.numClusters; i++) {
            Assert.assertEquals(reader.hasNext(), true);
            Assert.assertEquals(reader.next().booleanValue(), expectedPfs[i]);
        }

        Assert.assertEquals(false, reader.hasNext());
    }

    @Test
    void readFakedFile() throws Exception {
        final File fakeFile = File.createTempFile("FilterFileFakerTest", ".filter");
        fakeFile.deleteOnExit();

        new FilterFileFaker().fakeFile(fakeFile, 100);
        final FilterFileReader reader = new FilterFileReader(fakeFile);
        Assert.assertEquals(reader.numClusters, 1);         // A faked file has one cluster - with PF = false.
        Assert.assertTrue(reader.hasNext());
        Assert.assertFalse(reader.next());
        Assert.assertFalse(reader.hasNext());
    }

    @Test(expectedExceptions = NoSuchElementException.class)
    public void readPastEnd() {
        final FilterFileReader reader = new FilterFileReader(PASSING_FILTER_FILE);
        for(int i = 0; i < reader.numClusters; i++) {
            reader.next();
        }

        Assert.assertEquals(false, reader.hasNext());
        reader.next();
    }

    @DataProvider(name="failingFilesForSAMException")
    public Object[][] failingFilesForSAMException() {
        return new Object[][] {
                {"pf_nonExistentFile.filter"}
        };
    }

    @DataProvider(name="failingFilesForReaderException")
    public Object[][] failingFilesForReaderException() {
        return new Object[][] {
            {"pf_failing1.filter"},
            {"pf_failing2.filter"},
            {"pf_tooLarge.filter"},
            {"pf_tooShort.filter"},
            {"pf_badOpeningBytes.filter"},
            {"pf_badVersionBytes.filter"}
        };
    }

    @Test(dataProvider = "failingFilesForSAMException", expectedExceptions = SAMException.class)
    public void readInvalidValuesForSAMException(final String failingFile) {
        final FilterFileReader reader = new FilterFileReader(new File(TEST_DATA_DIR, failingFile));
        while(reader.hasNext()) {
            reader.next();
        }
    }

    @Test(dataProvider = "failingFilesForReaderException", expectedExceptions = IlluminaReaderException.class)
    public void readInvalidValuesForPicardException(final String failingFile) {
        final FilterFileReader reader = new FilterFileReader(new File(TEST_DATA_DIR, failingFile));
        while(reader.hasNext()) {
            reader.next();
        }
    }
}