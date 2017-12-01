package org.broadinstitute.hellbender.tools.exome.conversion.titanconversion;


import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollection;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollectionUtils;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCountCollection;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.File;
import java.io.IOException;

/**
 * This class is a helper class for CLIs to support TITAN
 *
 * http://genome.cshlp.org/content/early/2014/07/24/gr.180281.114.abstract
 */
public class TitanFileConverter {
    private static final Logger logger = LogManager.getLogger(TitanFileConverter.class);

    private TitanFileConverter() {}

    /**
     * Create a het pulldown file that is compatible with TITAN.
     *
     * @param hetPulldown Readable file from one of the het pulldown tools
     * @param outputFile Not {@code null}
     */
    public static void convertHetPulldownToTitanHetFile(final File hetPulldown, final File outputFile) {
        IOUtils.canReadFile(hetPulldown);
        try {
            final AllelicCountCollection acc = new AllelicCountCollection(hetPulldown);
            final TitanAllelicCountWriter titanAllelicCountWriter = new TitanAllelicCountWriter(outputFile);
            titanAllelicCountWriter.writeAllRecords(acc.getCounts());
            titanAllelicCountWriter.close();
        } catch (final IOException ioe) {
            throw new UserException.BadInput("Bad output file: " + outputFile);
        } catch (final IllegalArgumentException iae) {
            logger.warn("Cannot convert pulldown with basic verbosity (i.e., missing ref/alt nucleotide information) to TITAN-compatible file.  An empty file will be output.");
        }
    }

    /**
     * Create a target file that is compatible with TITAN.
     *
     * @param tnFile Readable file from {@link org.broadinstitute.hellbender.tools.exome.NormalizeSomaticReadCounts}
     * @param outputFile Not {@code null}
     */
    public static void convertCRToTitanCovFile(final File tnFile, final File outputFile) {
        IOUtils.canReadFile(tnFile);
        try {
            final ReadCountCollection rcc = ReadCountCollectionUtils.parse(tnFile);
            final TitanCopyRatioEstimateWriter titanCopyRatioEstimateWriter = new TitanCopyRatioEstimateWriter(outputFile);
            titanCopyRatioEstimateWriter.writeAllRecords(rcc.records());
            titanCopyRatioEstimateWriter.close();
        } catch (final IOException ioe) {
            throw new UserException.BadInput("Bad output file: " + outputFile);
        }
    }
}
