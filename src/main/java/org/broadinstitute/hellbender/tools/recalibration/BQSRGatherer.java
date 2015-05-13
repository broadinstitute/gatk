package org.broadinstitute.hellbender.tools.recalibration;

import org.apache.commons.collections4.CollectionUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.commandline.Gatherer;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.report.GATKReport;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.*;

/**
 * Merges BQSR Recaliration tables.
 */
public final class BQSRGatherer extends Gatherer {

    private static final Logger logger = LogManager.getLogger(BQSRGatherer.class);
    private static final String EMPTY_INPUT_LIST = "list of inputs files is empty or there is no usable data in any input file";
    private static final String MISSING_READ_GROUPS = "Missing read group(s)";

    @Override
    public void gather(final List<File> inputs, final File output) throws IOException {
        try (final PrintStream outputFile = IOUtils.makePrintStreamMaybeGzipped(output)) {
            final GATKReport report = gatherReport(inputs);
            report.print(outputFile);
        }
    }

    /**
     * Gathers the input recalibration reports into a single report.
     *
     * @param inputs Input recalibration GATK reports
     * @return gathered recalibration GATK report
     */
    public static GATKReport gatherReport(final List<File> inputs) {
        final SortedSet<String> allReadGroups = new TreeSet<>();
        final Map<File, Set<String>> inputReadGroups = new LinkedHashMap<>();

        // Get the read groups from each input report
        for (final File input : inputs) {
            final Set<String> readGroups = RecalibrationReport.getReadGroups(input);
            inputReadGroups.put(input, readGroups);
            allReadGroups.addAll(readGroups);
        }

        // Log the read groups that are missing from specific inputs
        for (Map.Entry<File, Set<String>> entry: inputReadGroups.entrySet()) {
            final File input = entry.getKey();
            final Set<String> readGroups = entry.getValue();
            if (allReadGroups.size() != readGroups.size()) {
                // Since this is not completely unexpected, more than debug, but less than a proper warning.
                logger.info(MISSING_READ_GROUPS + ": " + input.getAbsolutePath());
                for (final Object readGroup: CollectionUtils.subtract(allReadGroups, readGroups)) {
                    logger.info("  " + readGroup);
                }
            }
        }

        RecalibrationReport generalReport = null;
        for (File input : inputs) {
            final RecalibrationReport inputReport = new RecalibrationReport(input, allReadGroups);
            if( inputReport.isEmpty() ) {
                continue;
            }

            if (generalReport == null) {
                generalReport = inputReport;
            } else {
                generalReport.combine(inputReport);
            }
        }
        if (generalReport == null) {
            throw new GATKException(EMPTY_INPUT_LIST);
        }

        generalReport.calculateQuantizedQualities();

        return generalReport.createGATKReport();
    }
}
