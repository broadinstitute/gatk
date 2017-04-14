package org.broadinstitute.hellbender.utils.test;

import com.intel.genomicsdb.GenomicsDBImporter;
import org.apache.commons.io.FileUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;
import java.util.function.BiConsumer;


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
        final GenomicsDBImporter importer = new GenomicsDBImporter(loader.getAbsolutePath());
        importer.write();
        test.run();
    }

    /**
     * test action which is allowed to throw IOException since Runnable doesn't allow checked exceptions
     */
    @FunctionalInterface
    public interface GenomicsDBTestAction {

        void run() throws IOException;
    }

    public static <T> void assertCondition(Iterable<T> actual, Iterable<T> expected, BiConsumer<T,T> assertion){
        final Iterator<T> iterActual = actual.iterator();
        final Iterator<T> iterExpected = expected.iterator();
        while(iterActual.hasNext() && iterExpected.hasNext()){
            assertion.accept(iterActual.next(), iterExpected.next());
        }
        if (iterActual.hasNext()){
            Assert.fail("actual is longer than expected with at least one additional element: " + iterActual.next());
        }
        if (iterExpected.hasNext()){
            Assert.fail("actual is shorter than expected, missing at least one element: " + iterExpected.next());
        }

    }
}
