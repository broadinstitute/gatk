package org.broadinstitute.hellbender.metrics;

import htsjdk.samtools.SAMException;
import htsjdk.samtools.metrics.MetricsFile;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.*;
import java.nio.file.Path;

/**
 * Utility methods for dealing with {@link MetricsFile} and related classes.
 */
public final class MetricsUtils {

    private MetricsUtils(){} //don't instantiate this utility class

    /**
     * Write a {@link MetricsFile} to the given path, can be any destination supported by {@link IOUtils#openOutputStream(Path)} )}
     * @param metricsFile a {@link MetricsFile} object to write to disk
     * @param metricsOutputPath the path (or uri) to write the metrics to
     */
    public static void saveMetrics(final MetricsFile<?, ?> metricsFile, String metricsOutputPath) {
        try(Writer out = new BufferedWriter(new OutputStreamWriter(IOUtils.openOutputStream(IOUtils.getPath(metricsOutputPath))))) {
            metricsFile.write(out);
        } catch (IOException | SAMException e ){
            throw new UserException.CouldNotCreateOutputFile("Could not write metrics to file: " + metricsOutputPath, e);
        }
    }

}
