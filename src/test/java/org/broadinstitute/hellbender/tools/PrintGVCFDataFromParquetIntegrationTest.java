package org.broadinstitute.hellbender.tools;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

/**
 * Created by akiezun on 2/10/16.
 */
public class PrintGVCFDataFromParquetIntegrationTest extends CommandLineProgramTest {
    @Test
    public void test() throws IOException {
        final File input = new File(".", "gvcfs.ANKEF1.parquet");
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--file");
        args.add(input.getCanonicalPath());

        runCommandLine(args.getArgsList());
    }
}
