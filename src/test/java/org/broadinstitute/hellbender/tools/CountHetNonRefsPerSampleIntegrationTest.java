package org.broadinstitute.hellbender.tools;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Collections;
import java.util.Map;

public class CountHetNonRefsPerSampleIntegrationTest extends CommandLineProgramTest {

    @Test
    public void testCountHetsEmptyResults() throws Exception {
        String fileIn = "no_genotypes.vcf";
        final File ORIG_FILE = new File(getToolDataDir(), fileIn);
        final String[] args = new String[]{
            "--variant" ,  ORIG_FILE.getAbsolutePath(),

        };
        final Map<String, Long> res = (Map<String, Long>)this.runCommandLine(args);
        Assert.assertEquals(res, Collections.emptyMap());
    }

}
