package org.broadinstitute.hellbender.testutils;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.runtime.ProcessController;

import java.util.*;

/**
 * utilities to help run tests on a Dataproc cluster
 *
 * These require you to have a working installation of the google-cloud-sdk (gcloud)
 * They also require that you have logged in to gcloud with a user that has appropriate permissions to create dataproc
 * clusters in your project and launch jobs on dataproc clusters
 *
 * See README.md for instructions on how to authenticate gcloud so it works with gatk
 * See https://cloud.google.com/dataproc/docs/concepts/iam for detailed information about the necessary permissions
 *
 */
public final class DataprocTestUtils {
    private static final Logger logger = LogManager.getLogger(DataprocTestUtils.class);

    /**
     * environment variable that specifies where to find gcloud utilities
     */
    public static final String CLOUD_HOME_ENVIRONMENT_VARIABLE_NAME = "GCLOUD_HOME";

    /**
     * environment variable used to override cluster creation and use one with the given name instead
     */
    public static final String CLUSTER_NAME_ENVIRONMENT_VARIABLE_NAME = "HELLBENDER_DATAPROC_CLUSTER_NAME";

    private DataprocTestUtils(){}

    /**
     * @return the value of the GCLOUD_HOME environment variable if present, or "gcloud"
     */
    private static String getGCloudPath(){
        final String gcloudBinDir = System.getenv(CLOUD_HOME_ENVIRONMENT_VARIABLE_NAME);
        return gcloudBinDir == null ? "gcloud" : gcloudBinDir + "/gcloud";
    }

    /**
     *  get the name of a dataproc cluster to dispatch jobs to
     *  checks if the {@link #CLUSTER_NAME_ENVIRONMENT_VARIABLE_NAME} environment variable is set, if not, it will
     *  create a new short duration cluster to use
     *
     *  $$$$$  This method costs you money!  $$$$$
     *
     * @return the name of the cluster
     */
    public static String getTestCluster(){
        final String clusterName = System.getenv(CLUSTER_NAME_ENVIRONMENT_VARIABLE_NAME);
        return clusterName != null ? clusterName : createTestCluster();
    }

    /**
     * create a minimal dataproc cluster which will shutdown automatically after 10 minutes of idling or 30 minutes
     * this is accomplished by shelling out to gcloud
     *
     * $$$$$  This method costs you money!  $$$$$
     *
     * @return the name of the cluster
     */
    private static String createTestCluster() {
        logger.info("Starting Dataproc cluster creation.");
        final String clusterName = "gatk-test-" + UUID.randomUUID();
        final String[] command = new String[]{
            getGCloudPath(), "beta", "dataproc", "clusters", "create",
            "--max-idle", "10m",
            "--max-age", "30m",
            "--num-workers", "2",
            "--master-machine-type", "n1-highmem-2",
            "--worker-machine-type", "n1-highmem-2",
            clusterName
        };
        BaseTest.runProcess(ProcessController.getThreadLocal(), command, "Couldn't create dataproc cluster");
        logger.info("Finished Dataproc cluster creation: cluster named" + clusterName);
        return clusterName;
    }

    /**
     * run a gatk spark tool on a dataproc cluster and throw if it fails
     *
     * this makes use of the gatk script and gcloud and will not work if those are not available
     * projects depending on GATK as a library will probably have to implement an equivalent method instead of using this one
     *
     * @param tool name of the tool
     * @param args arguments to the tool
     * @param clusterName name of the cluster
     */
    public static void launchGatkTool(final String tool, final List<String> args, final String clusterName){
        //set the jar staging directory so jar caching works
        final Map<String, String> env = new HashMap<>(System.getenv());
        final String gcpTestStaging = BaseTest.getGCPTestStaging();
        env.put("GATK_GCS_STAGING", gcpTestStaging);

        final List<String> command = new ArrayList<>();
        command.add("./gatk");
        command.add(tool);
        command.addAll(args);
        command.add("--");
        command.add("--spark-runner"); command.add("GCS");
        command.add("--cluster"); command.add(clusterName);
        final String[] commandArray = command.toArray(new String[command.size()]);
        BaseTest.runProcess(ProcessController.getThreadLocal(), commandArray, env, "gatk run on dataproc failed");
    }
}
