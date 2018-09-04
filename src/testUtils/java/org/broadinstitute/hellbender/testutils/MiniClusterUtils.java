package org.broadinstitute.hellbender.testutils;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.hdfs.MiniDFSCluster;
import org.testng.Assert;

import java.io.File;
import java.io.IOException;
import java.util.UUID;

/**
 * Utilities to help manage a {@link MiniDFSCluster} for use in tests.
 */
public final class MiniClusterUtils {

    /**
     * @return a new empty cluster with the given Configuration
     * @throws IOException
     */
    public static MiniDFSCluster getMiniCluster(Configuration otherConf) throws IOException {
        final File baseDir = BaseTest.createTempDir("minicluster_storage");
        final Configuration configuration;
        if (otherConf != null) {
            configuration = new Configuration(otherConf);
        } else {
            configuration = new Configuration();
        }
        configuration.set(MiniDFSCluster.HDFS_MINIDFS_BASEDIR, baseDir.getAbsolutePath());
        return new MiniDFSCluster.Builder(configuration).build();
    }

    /**
     * @return a new empty cluster
     * @throws IOException
     */
    public static MiniDFSCluster getMiniCluster() throws IOException {
        return getMiniCluster(null);
    }

    /**
     * @return a unique path in the filesystem on the given cluster
     * @throws IOException
     */
    public static Path getTempPath(MiniDFSCluster cluster, String prefix, String extension) throws IOException {
        return new Path(cluster.getFileSystem().getWorkingDirectory(), prefix + "_" + UUID.randomUUID() + extension);
    }

    /**
     * @return path to the working directory in the given cluster
     * @throws IOException
     */
    public static Path getWorkingDir(MiniDFSCluster cluster) throws IOException {
        return cluster.getFileSystem().getWorkingDirectory();
    }

    /**
     * shut down the cluster
     */
    public static void stopCluster(MiniDFSCluster cluster){
        if(cluster != null){
            cluster.shutdown(true);
        }
    }

    /**
     * Create a new isolated {@link MiniDFSCluster}, run a {@link MiniClusterTest} on it and then shut it down.
     * @param test a function to run on the cluster, this should indicate success or failure by using {@link Assert}
     * @throws Exception
     */
    public static void runOnIsolatedMiniCluster(MiniClusterTest test) throws Exception {
        MiniDFSCluster cluster = null;
        try {
            cluster = getMiniCluster();
            test.test(cluster);
        } finally {
            stopCluster(cluster);
        }
    }

    /**
     * An interface for writing tests that run on a minicluster.
     * Implementers should write test to make use of the passed in cluster
     */
    @FunctionalInterface
    public interface MiniClusterTest {
        /**
         * A test to be run using an hdfs MiniCluster
         * It's alright for this to make destructive changes to the cluster since it is given it's own isolated setup.
         *
         * Test failure is indicated use standard {@link Assert} methods
         *
         * @param cluster an isolated MiniDFSCluster instance
         * @throws Exception in order to allow any checked exception to be thrown
         */
        void test(MiniDFSCluster cluster) throws Exception;
    }


}
