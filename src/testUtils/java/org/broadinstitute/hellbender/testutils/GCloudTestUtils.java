package org.broadinstitute.hellbender.testutils;

import com.google.common.base.Strings;

public final class GCloudTestUtils {
    /**
     *  This is a public requester pays bucket owned by the broad-gatk-test project.
     *  It must be owned by a different project than the service account doing the testing or the test may fail because it can access the
     *  file directly through alternative permissions.
     */
    public static final String REQUESTER_PAYS_BUCKET_DEFAULT = "gs://hellbender-requester-pays-test/";

    public static final String TEST_INPUTS_DEFAULT = "gs://hellbender/test/resources/";
    public static final String TEST_STAGING_DEFAULT = "gs://hellbender-test-logs/staging/";
    public static final String TEST_PROJECT_DEFAULT = "broad-dsde-dev";


    /**
     * A publicly readable GCS bucket set as requester pays, this should not be owned by the same project that is set
     * as {@link #getTestProject()} or the tests for requester pays access may be invalid.
     *
     * @return HELLBENDER_REQUESTER_PAYS_BUCKET env. var if defined, {@value GCloudTestUtils#REQUESTER_PAYS_BUCKET_DEFAULT}.
     */
    public static String getRequesterPaysBucket() {
        return getEnvironmentVariable("HELLBENDER_REQUESTER_PAYS_BUCKET", REQUESTER_PAYS_BUCKET_DEFAULT);
    }

    private static String getEnvironmentVariable(final String variableName, final String defaultValue) {
        final String valueFromEnvironment = System.getProperty(variableName);
        return valueFromEnvironment == null || valueFromEnvironment.isEmpty()? defaultValue : valueFromEnvironment;
    }

    /**
     * name of the google cloud project that stores the data and will run the code
     *
     * @return HELLBENDER_TEST_PROJECT env. var if defined or {@value #TEST_PROJECT_DEFAULT}
     */
    public static String getTestProject() {
        return getEnvironmentVariable("HELLBENDER_TEST_PROJECT", TEST_PROJECT_DEFAULT);
    }

    /**
     * A writable GCS path where java files can be cached and temporary test files can be written,
     * of the form gs://bucket/, or gs://bucket/path/.
     *
     * @return HELLBENDER_TEST_STAGING env. var if defined, or {@value #TEST_STAGING_DEFAULT}
     */
    public static String getTestStaging() {
        return getEnvironmentVariable("HELLBENDER_TEST_STAGING", TEST_STAGING_DEFAULT);
    }

    /**
     * A GCS path where the test inputs are stored.
     * <p>
     * The value of HELLBENDER_TEST_INPUTS should end in a "/" (for example, "gs://hellbender/test/resources/")
     *
     * @return HELLBENDER_TEST_INPUTS env. var if defined or {@value #TEST_INPUTS_DEFAULT}.
     */
    public static String getTestInputPath() {
        return getEnvironmentVariable("HELLBENDER_TEST_INPUTS", TEST_INPUTS_DEFAULT);
    }

}
