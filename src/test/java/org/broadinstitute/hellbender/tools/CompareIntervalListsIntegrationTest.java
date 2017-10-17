package org.broadinstitute.hellbender.tools;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;

public class CompareIntervalListsIntegrationTest extends CommandLineProgramTest{

    private static final String INTERVALS_1 = "intervals1.intervals";
    private static final String INTERVLS_2  = "intervals2.intervals";
    private static final String INTERVALS_3 = "intervals3.intervals"; //same territory as

    @DataProvider
    public Object[][] getIntervalFiles(){
        return new Object[][]{
                {getTestFile(INTERVALS_1), getTestFile(INTERVALS_1), true},
                {getTestFile(INTERVLS_2), getTestFile(INTERVALS_1), false},
                {getTestFile(INTERVALS_1), getTestFile(INTERVLS_2), false},
                {getTestFile(INTERVLS_2), getTestFile(INTERVLS_2), true},
                {getTestFile(INTERVALS_3), getTestFile(INTERVALS_1), true},
                {getTestFile(INTERVALS_3), getTestFile(INTERVLS_2), false},
                {getTestFile(INTERVALS_1), getTestFile(INTERVALS_3), true},
        };
    }

    @Test(dataProvider = "getIntervalFiles")
    public void testCompareIntervals(File left, File right, boolean expectedMatch) {
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addReference(new File(hg19MiniReference))
                .addFileArgument("L", left)
                .addFileArgument("L2", right);

        if(expectedMatch) {
            runCommandLine(args);
        } else {
            Assert.expectThrows(UserException.class, () -> runCommandLine(args));
        }

    }
}