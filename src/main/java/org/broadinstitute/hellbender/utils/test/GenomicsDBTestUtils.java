package org.broadinstitute.hellbender.utils.test;

import genomicsdb.VCF2TileDB;
import org.apache.commons.io.FileUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.File;
import java.io.IOException;


/**
 * Utils for running tests using TileDB
 */
public final class GenomicsDBTestUtils {
    private GenomicsDBTestUtils() {};


    /**
     * creates a new GenomicsDB workspace and loads data into it before running a test
     *
     * @param test this will be run after the data has been loaded successfully
     * @param workspace where to create the workspace, if something exists here already it will be deleted
     * @param loader loader.json file to use to load data from, the workspace in this file must match workspace
     */
    public synchronized static void runOnGenomicsDBArray(GenomicsDBTestAction test, File workspace, File loader) throws IOException {
        //delete any existing workspace
        FileUtils.deleteQuietly(workspace);

        IOUtils.deleteRecursivelyOnExit(workspace);
        final VCF2TileDB vcf2TileDB = new VCF2TileDB(loader.getAbsolutePath());
        vcf2TileDB.write();
        test.run();
    }

    /**
     * test action which is allowed to throw IOException since Runnable doesn't allow checked exceptions
     */
    @FunctionalInterface
    public interface GenomicsDBTestAction {

        void run() throws IOException;
    }
}
