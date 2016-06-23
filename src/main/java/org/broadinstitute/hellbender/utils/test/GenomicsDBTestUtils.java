package org.broadinstitute.hellbender.utils.test;

import org.apache.commons.io.FileUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.runtime.ProcessController;
import org.broadinstitute.hellbender.utils.runtime.ProcessOutput;
import org.broadinstitute.hellbender.utils.runtime.ProcessSettings;
import org.testng.Assert;
import org.testng.SkipException;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;


/**
 * In order to run tests using this class, the environment variable {@link GenomicsDBTestUtils#GENOMICSDB_BIN_PATH_PROPERTY_NAME} must be set
 * and point to a location containing the genomics db binaries {@link GenomicsDBTestUtils#CREATE_TILEDB_WORKSPACE_BINARY_NAME} and
 * {@link GenomicsDBTestUtils#VCF2TILEDB_BINARY_NAME}
 *
 * For information about how to build these binaries for your system
 * see https://github.com/Intel-HLS/GenomicsDB/wiki/Compiling-GenomicsDB
 */
public final class GenomicsDBTestUtils {
    private GenomicsDBTestUtils() {};

    /**
     * environment variable that is used to lookup the path to GenomicsDB binaries
     */
    public static final String GENOMICSDB_BIN_PATH_PROPERTY_NAME = "GATK_GENOMICSDB_BIN";
    /**
     * names of the two binaries that are used to setup a GenomicsDB array
     */
    public static final String CREATE_TILEDB_WORKSPACE_BINARY_NAME = "create_tiledb_workspace";
    public static final String VCF2TILEDB_BINARY_NAME = "vcf2tiledb";

    /**
     * creates a new GenomicsDB workspace and loads data into it before running a test
     *
     * In order to run tests using this method, the environment variable {@link GenomicsDBTestUtils#GENOMICSDB_BIN_PATH_PROPERTY_NAME}
     * must be set and point to a location containing the genomics db binaries {@link GenomicsDBTestUtils#CREATE_TILEDB_WORKSPACE_BINARY_NAME} and
     * {@link GenomicsDBTestUtils#VCF2TILEDB_BINARY_NAME}
     *
     * For information about how to build these binaries for your system
     * see https://github.com/Intel-HLS/GenomicsDB/wiki/Compiling-GenomicsDB
     *
     * if {@link GenomicsDBTestUtils#GENOMICSDB_BIN_PATH_PROPERTY_NAME} is not set, then any tests using this method will be skipped*
     *
     * @param test this will be run after the data has been loaded successfully
     * @param workspace where to create the workspace, if something exists here already it will be deleted
     * @param loader loader.json file to use to load data from, the workspace in this file must match workspace
     */
    public synchronized static void runOnGenomicsDBArray(GenomicsDBTestAction test, File workspace, File loader) throws IOException {
        //delete any existing workspace
        FileUtils.deleteQuietly(workspace);

        final String gendbBinariesPath = System.getenv(GENOMICSDB_BIN_PATH_PROPERTY_NAME);
        if (gendbBinariesPath == null) {
            throw new SkipException("Skipping GenomicsDB test because " + GENOMICSDB_BIN_PATH_PROPERTY_NAME + " wasn't specified");
        }

        IOUtils.deleteRecursivelyOnExit(workspace);
        final ProcessController processController = new ProcessController();

        runProcess(processController, new String[]{gendbBinariesPath + "/" + CREATE_TILEDB_WORKSPACE_BINARY_NAME, workspace.getAbsolutePath()});
        Assert.assertTrue(workspace.exists(), "Workspace was not created at " + workspace.getAbsolutePath());
        runProcess(processController, new String[]{gendbBinariesPath + "/" + VCF2TILEDB_BINARY_NAME, loader.getAbsolutePath()});
        test.run();
    }

    private static void runProcess(ProcessController processController, String[] command) {
        final ProcessSettings prs = new ProcessSettings(command);
        prs.getStderrSettings().printStandard(true);
        prs.getStdoutSettings().printStandard(true);
        final ProcessOutput output = processController.exec(prs);
        Assert.assertEquals(output.getExitValue(), 0, "Process exited with non-zero value. Command: "+ Arrays.toString(command) + "\n");
    }

    /**
     * test action which is allowed to throw IOException since Runnable doesn't allow checked exceptions
     */
    @FunctionalInterface
    public interface GenomicsDBTestAction {

        void run() throws IOException;
    }
}
