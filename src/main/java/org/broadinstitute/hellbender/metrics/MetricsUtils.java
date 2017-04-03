package org.broadinstitute.hellbender.metrics;

import htsjdk.samtools.SAMException;
import htsjdk.samtools.metrics.MetricsFile;
import org.broadinstitute.hellbender.engine.AuthHolder;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;

import java.io.PrintWriter;

/**
 * Utility methods for dealing with {@link MetricsFile} and related classes.
 */
public final class MetricsUtils {

    private MetricsUtils(){} //don't instantiate this utility class

    /**
     * Write a {@link MetricsFile} to the given path, can be any destination supported by {@link BucketUtils#createFile(String)}
     * @param metricsFile a {@link MetricsFile} object to write to disk
     * @param metricsOutputPath the path (or uri) to write the metrics to
     * @param authHolder authentication for remote paths, can be null if the path is not remote
     */
    public static void saveMetrics(final MetricsFile<?, ?> metricsFile, String metricsOutputPath, AuthHolder authHolder) {
        try(PrintWriter out = new PrintWriter(BucketUtils.createFile(metricsOutputPath))) {
            metricsFile.write(out);
        } catch (SAMException e ){
            throw new UserException.CouldNotCreateOutputFile("Could not write metrics to file: " + metricsOutputPath, e);
        }
    }

}
