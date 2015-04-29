package org.broadinstitute.hellbender.tools.picard.illumina.parser;

import org.testng.annotations.DataProvider;

import java.io.File;
import java.util.*;

import static htsjdk.samtools.util.CollectionUtil.makeList;

/**
* @author jburke@broadinstitute.org
*/
public class CycleIlluminaFileMapTest {
    //TODO: REVAMP THIS
    private static final File TEST_DATA_DIR = new File("src/test/resources/org/broadinstitute/hellbender/tools/picard/illumina/25T8B25T/Data/Intensities/BaseCalls/L001");
    private static final File ZERO_LENGTH_TEST_DATA_DIR = new File("src/test/resources/org/broadinstitute/hellbender/tools/picard/illumina/25T8B25T/Data/Intensities/BaseCalls/L002");
    private static final int [] ALL_CYCLES = {1,2,3,4};

    private static String laneToDir(int lane) {
        String outStr = String.valueOf(lane);
        while(outStr.length() < 3) {
            outStr = "0" + outStr;
        }
        return "L" + outStr;
    }

    @DataProvider(name = "iteratorTestData")
    public Object[][] iteratorTestData() {
        return new Object[][] {
            new Object[] {
                TEST_DATA_DIR, 1, 1101, ".bcl",
                    makeList(new File(TEST_DATA_DIR + "/C1.1", "s_1_1101.bcl"),
                             new File(TEST_DATA_DIR + "/C2.1", "s_1_1101.bcl"),
                             new File(TEST_DATA_DIR + "/C3.1", "s_1_1101.bcl"),
                             new File(TEST_DATA_DIR + "/C4.1", "s_1_1101.bcl")),
                    ALL_CYCLES
            },
            new Object[] {
                ZERO_LENGTH_TEST_DATA_DIR, 1, 1101, ".bcl", new ArrayList<File>(), new int[]{},

            },
            new Object[] {
                TEST_DATA_DIR, 1, 1201, ".bcl", new ArrayList<File>(), new int[]{}
            },

            new Object[] {
                TEST_DATA_DIR, 2, 1101, ".bcl", new ArrayList<File>(), new int[]{}
            },
        };
    }
}
