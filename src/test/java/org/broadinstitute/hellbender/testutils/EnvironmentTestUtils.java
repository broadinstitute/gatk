package org.broadinstitute.hellbender.testutils;

import org.broadinstitute.hellbender.utils.runtime.ScriptExecutor;

import java.util.Map;
import java.util.Properties;

/**
 * Utility class for running tests with specific system properties set.
 */
public class EnvironmentTestUtils {
    /**
     * Runs the given Runnable with the GATK Lite Docker property set to true.
     *
     * @param toRun      The Runnable to execute.
     */
    public static void checkWithGATKDockerPropertySet(Runnable toRun) {
        runWithSystemProperties(toRun, Map.of(ScriptExecutor.GATK_LITE_DOCKER_ENV_VAR, "true"));
    }
   

    /**
     * Runs the given Runnable with the specified system properties set.
     *
     * @param toRun      The Runnable to execute.
     * @param properties A map of system properties to set before running the Runnable.
     */
    public static void runWithSystemProperties(Runnable toRun, Map<String, String> properties) {
        final Properties originalProperties = System.getProperties();
        try {
            properties.forEach(System::setProperty);
            toRun.run();
        } finally {
            properties.keySet().forEach(System::clearProperty);
            System.setProperties(originalProperties);
        }
    }
}
