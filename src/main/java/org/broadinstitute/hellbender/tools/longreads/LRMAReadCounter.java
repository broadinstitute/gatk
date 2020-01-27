package org.broadinstitute.hellbender.tools.longreads;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.hellbender.cmdline.programgroups.CoverageAnalysisProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

@CommandLineProgramProperties(
        summary = "Counts up all reads with specific sets of flags set in one pass.",
        oneLineSummary = "Counts up all reads with specific sets of flags set in one pass.",
        programGroup = CoverageAnalysisProgramGroup.class
)
@ExperimentalFeature
public class LRMAReadCounter extends ReadWalker {
    private static final Logger logger = LogManager.getLogger(LRMAReadCounter.class);

    //==================================================================================================================
    // Public Static Members:

    //==================================================================================================================
    // Private Static Members:

    //==================================================================================================================
    // Public Members:
    @Argument(
            fullName  = "output-csv-file",
            optional = true,
            doc = "Output the results to the given CSV file name.")
    private String outputCsvFileName = "";

    //==================================================================================================================
    // Private Members:
    private long totalReadsCount = 0;
    private long unmappedReadsCount = 0;
    private long mappedReadsCount = 0;
    private long primaryReadsCount = 0;
    private long secondaryReadsCount = 0;
    private long supplementaryReadsCount = 0;

    static final Pattern PAC_BIO_READ_NAME_PATTERN = Pattern.compile("(m[0-9]*_[0-9]*_[0-9]*)/([0-9]*)/(.*)");
    private long numUniqueZmws = 0;
    private HashSet<Integer> zmwSet = new HashSet<>();


    //==================================================================================================================
    // Constructors:

    //==================================================================================================================
    // Override Methods:

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        return Collections.emptyList();
    }

    @Override
    public void apply(final GATKRead read, final ReferenceContext referenceContext, final FeatureContext featureContext) {
        accountForZmw(read);
        accountForMappingFlags(read);
    }

    @Override
    public Object onTraversalSuccess() {
        printStats();

        if ( !outputCsvFileName.isEmpty() ) {
            writeStatsToCsvFile();
        }

        return null;
    }

    //==================================================================================================================
    // Static Methods:

    //==================================================================================================================
    // Instance Methods:

    private void writeStatsToCsvFile() {
        final Path path = Paths.get(outputCsvFileName);

        logger.info("Writing CSV file to: {}", path.toUri().toString());

        final StringBuilder sb = new StringBuilder();

        // Write field names:
        sb.append("Num Unique ZMWs"); sb.append(',');
        sb.append("Unmapped Reads"); sb.append(',');
        sb.append("Secondary Aligned Reads"); sb.append(',');
        sb.append("Supplementary Aligned Reads"); sb.append(',');
        sb.append("Primary Aligned Reads"); sb.append(',');
        sb.append("Mapped Reads"); sb.append(',');
        sb.append("Total Reads");

        sb.append('\n');

        // Write values:
        sb.append(numUniqueZmws); sb.append(',');
        sb.append(unmappedReadsCount); sb.append(',');
        sb.append(secondaryReadsCount); sb.append(',');
        sb.append(supplementaryReadsCount); sb.append(',');
        sb.append(primaryReadsCount); sb.append(',');
        sb.append(mappedReadsCount); sb.append(',');
        sb.append(totalReadsCount);

        sb.append('\n');

        try {
            Files.write(path, sb.toString().getBytes());
        }
        catch (final IOException ex ) {
            throw new UserException("Could not write out stats csv!", ex);
        }
    }

    private void printStats() {
        // Print out our stats:
        logger.info("----");
        logger.info("");
        logger.info("Num Unique ZMWs:             {}", numUniqueZmws);
        logger.info("");
        logger.info("Unmapped Reads:              {}", unmappedReadsCount);
        logger.info("Secondary Aligned Reads:     {}", secondaryReadsCount);
        logger.info("Supplementary Aligned Reads: {}", supplementaryReadsCount);
        logger.info("Primary Aligned Reads:       {}", primaryReadsCount);
        logger.info("");
        logger.info("Mapped Reads:                {}", mappedReadsCount);
        logger.info("");
        logger.info("Total Reads:                 {}", totalReadsCount);
        logger.info("");
        logger.info("----");
    }

    private void accountForZmw(final GATKRead read) {
        final Matcher matcher = PAC_BIO_READ_NAME_PATTERN.matcher(read.getName());
        if ( matcher.matches() ) {
            final Integer zmw = Integer.valueOf(matcher.group(2));
            if (!zmwSet.contains(zmw)) {
                zmwSet.add(zmw);
                ++numUniqueZmws;
            }
        }
    }

    private void accountForMappingFlags(final GATKRead read) {
        ++totalReadsCount;

        if ( read.isUnmapped() ) {
            ++unmappedReadsCount;
        }
        else {
            ++mappedReadsCount;

            if (read.isSecondaryAlignment() &&  !read.isSupplementaryAlignment()) {
                ++secondaryReadsCount;
            }
            else if (!read.isSecondaryAlignment() &&  read.isSupplementaryAlignment()) {
                ++supplementaryReadsCount;
            }
            else {
                ++primaryReadsCount;
            }
        }
    }

    //==================================================================================================================
    // Helper Data Types:

}
