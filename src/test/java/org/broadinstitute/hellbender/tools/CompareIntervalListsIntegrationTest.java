package org.broadinstitute.hellbender.tools;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;

public class CompareIntervalListsIntegrationTest extends CommandLineProgramTest {

    private static final String INTERVALS_1 = "intervals1.intervals";
    private static final String INTERVALS_2 = "intervals2.intervals";
    private static final String INTERVALS_3 = "intervals3.intervals"; // same territory as
    private static final String INTERVALS_4 = "intervals4.intervals"; // empty file

    @DataProvider
    public Object[][] getIntervalFiles(){
        return new Object[][]{
                {getTestFile(INTERVALS_1), getTestFile(INTERVALS_1), true},
                {getTestFile(INTERVALS_2), getTestFile(INTERVALS_1), false},
                {getTestFile(INTERVALS_1), getTestFile(INTERVALS_2), false},
                {getTestFile(INTERVALS_2), getTestFile(INTERVALS_2), true},
                {getTestFile(INTERVALS_3), getTestFile(INTERVALS_1), true},
                {getTestFile(INTERVALS_3), getTestFile(INTERVALS_2), false},
                {getTestFile(INTERVALS_1), getTestFile(INTERVALS_3), true},
                {getTestFile(INTERVALS_1), getTestFile(INTERVALS_4), false},
                {getTestFile(INTERVALS_4), getTestFile(INTERVALS_1), false},
                {getTestFile(INTERVALS_4), getTestFile(INTERVALS_4), false},
        };
    }

    @Test(dataProvider = "getIntervalFiles")
    public void testCompareIntervals(File left, File right, boolean expectedMatch) {
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addReference(new File(hg19MiniReference))
                .add("L", left)
                .add("L2", right);

        if(expectedMatch) {
            runCommandLine(args);
        } else {
            Assert.expectThrows(UserException.class, () -> runCommandLine(args));
        }

    }
}