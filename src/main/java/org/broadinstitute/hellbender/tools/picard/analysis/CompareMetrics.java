package org.broadinstitute.hellbender.tools.picard.analysis;

import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.PicardCommandLineProgram;
import org.broadinstitute.hellbender.cmdline.PositionalArguments;
import org.broadinstitute.hellbender.cmdline.programgroups.QCProgramGroup;
import org.broadinstitute.hellbender.exceptions.GATKException;

import java.io.File;
import java.io.FileReader;
import java.util.*;

/**
 * Compare two metrics files.
 */
@CommandLineProgramProperties(
        usage = CompareMetrics.USAGE,
        usageShort = CompareMetrics.USAGE,
        programGroup = QCProgramGroup.class
)
public final class CompareMetrics extends PicardCommandLineProgram {

    static final String USAGE = "Compare two metrics files";

    @PositionalArguments(minElements = 2, maxElements = 2)
    public List<File> metricsFiles;

    private static final Log log = Log.getInstance(CompareMetrics.class);

    @Override
    protected Object doWork() {
        IOUtil.assertFilesAreReadable(metricsFiles);
        final MetricsFile<?, ?> metricsA = new MetricsFile<>();
        final MetricsFile<?, ?> metricsB = new MetricsFile<>();
        try {
            metricsA.read(new FileReader(metricsFiles.get(0)));
            metricsB.read(new FileReader(metricsFiles.get(1)));
            final boolean areEqual = metricsA.areMetricsEqual(metricsB) && metricsA.areHistogramsEqual(metricsB);
            final String status = areEqual ? "EQUAL" : "NOT EQUAL";
            log.info("Files " + metricsFiles.get(0) + " and " + metricsFiles.get(1) + "are " + status);
        } catch (final Exception e) {
            throw new GATKException(e.getMessage());
        }
        return null;
    }
}
