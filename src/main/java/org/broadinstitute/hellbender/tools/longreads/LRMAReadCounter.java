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
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

@CommandLineProgramProperties(
        summary = "Counts up all PacBio reads with specific sets of flags set in one pass.",
        oneLineSummary = "Counts up all PacBio reads with specific sets of flags set in one pass.",
        programGroup = CoverageAnalysisProgramGroup.class
)
@ExperimentalFeature
public class LRMAReadCounter extends ReadWalker {
    private static final Logger logger = LogManager.getLogger(LRMAReadCounter.class);

    //==================================================================================================================
    // Public Static Members:

    //==================================================================================================================
    // Private Static Members:

    /**
     * Simple data class to hold counts and supplementary for stats that we need on a per-movie basis.
     * NOTE: Pac-Bio specific.
     */
    private static class ReadStatHolder {
        private long         totalReadsCount = 0;
        private long         unmappedReadsCount = 0;
        private long         mappedReadsCount = 0;
        private long         primaryReadsCount = 0;
        private long         secondaryReadsCount = 0;
        private long         supplementaryReadsCount = 0;

        ReadStatHolder() {}

        long getTotalReadsCount() {
            return totalReadsCount;
        }

        long getUnmappedReadsCount() {
            return unmappedReadsCount;
        }

        long getMappedReadsCount() {
            return mappedReadsCount;
        }

        long getPrimaryReadsCount() {
            return primaryReadsCount;
        }

        long getSecondaryReadsCount() {
            return secondaryReadsCount;
        }

        long getSupplementaryReadsCount() {
            return supplementaryReadsCount;
        }

        String toCsvRow() {
            return unmappedReadsCount + "," + secondaryReadsCount + "," + supplementaryReadsCount + "," +
                    primaryReadsCount + "," + mappedReadsCount + "," + totalReadsCount;
        }

//        void incrementTotalReadsCount() { ++totalReadsCount; }
//        void incrementUnmappedReadsCount() { ++unmappedReadsCount; }
//        void incrementMappedReadsCount() { ++mappedReadsCount; }
//        void incrementPrimaryReadsCount() { ++primaryReadsCount; }
//        void incrementSecondaryReadsCount() { ++secondaryReadsCount; }
//        void incrementSupplementaryReadsCount() { ++supplementaryReadsCount; }

        void accountForMappingFlags(final GATKRead read) {
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
    }

    private static class PacBioMovieStatHolder extends ReadStatHolder {

        private final HashSet<Integer> zmwSet = new HashSet<>();
        private final String movieName;

        PacBioMovieStatHolder(final String movieName) {
            super();
            this.movieName = movieName;
        }

        String getMovieName() {
            return movieName;
        }

        void addZmw(final Integer zmw) {
            zmwSet.add(zmw);
        }

        int getNumUniqueZmws() {
            return zmwSet.size();
        }

        @Override
        String toCsvRow() {
            return super.toCsvRow() + "," + zmwSet.size();
        }
    }


    //==================================================================================================================
    // Public Members:
    @Argument(
            fullName  = "output-csv-file",
            optional = true,
            doc = "Output the results to the given CSV file name.")
    private String outputCsvFileName = "";

    //==================================================================================================================
    // Private Members:

    private static final Pattern                               PAC_BIO_READ_NAME_PATTERN = Pattern.compile("(m[0-9]*_[0-9]*_[0-9]*)/([0-9]*)/(.*)");
    private final LinkedHashMap<String, PacBioMovieStatHolder> movieNameToStatHolderMap  = new LinkedHashMap<>();
    private final ReadStatHolder                               overallReadStatHolder     = new ReadStatHolder();

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

        // PacBio reads are special:
        final Matcher matcher = PAC_BIO_READ_NAME_PATTERN.matcher(read.getName());
        if ( matcher.matches() ) {

            final String movieName = matcher.group(1);
            final Integer zmw = Integer.valueOf(matcher.group(2));
            if ( movieNameToStatHolderMap.get(movieName) == null ) {
                movieNameToStatHolderMap.put(movieName, new PacBioMovieStatHolder(movieName));
            }

            movieNameToStatHolderMap.get(movieName).addZmw(zmw);
            movieNameToStatHolderMap.get(movieName).accountForMappingFlags(read);
        }

        overallReadStatHolder.accountForMappingFlags(read);
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
        sb.append("Movie Name / Sample Name"); sb.append(',');
        sb.append("Unmapped Reads"); sb.append(',');
        sb.append("Secondary Aligned Reads"); sb.append(',');
        sb.append("Supplementary Aligned Reads"); sb.append(',');
        sb.append("Primary Aligned Reads"); sb.append(',');
        sb.append("Mapped Reads"); sb.append(',');
        sb.append("Total Reads"); sb.append(',');
        sb.append("Num Unique ZMWs"); sb.append(',');

        sb.append('\n');

        // Write overall stats first:
        sb.append("OVERALL");
        sb.append(',');
        sb.append(overallReadStatHolder.toCsvRow());
        sb.append(',');
        sb.append("N/A");
        sb.append('\n');

        // Write PacBio Stats:
        for ( final Map.Entry<String, PacBioMovieStatHolder> entry : movieNameToStatHolderMap.entrySet() ) {
            sb.append(entry.getKey());
            sb.append(',');
            sb.append(entry.getValue().toCsvRow());
            sb.append('\n');
        }

        try {
            Files.write(path, sb.toString().getBytes());
        }
        catch (final IOException ex ) {
            throw new UserException("Could not write out stats csv!", ex);
        }
    }

    private void logStats(final String name, final ReadStatHolder statHolder) {
        logger.info("----");
        logger.info("");
        logger.info("{}", name);
        logger.info("");

        logger.info("Unmapped Reads:              {}", statHolder.getUnmappedReadsCount());
        logger.info("Secondary Aligned Reads:     {}", statHolder.getSecondaryReadsCount());
        logger.info("Supplementary Aligned Reads: {}", statHolder.getSupplementaryReadsCount());
        logger.info("Primary Aligned Reads:       {}", statHolder.getPrimaryReadsCount());
        logger.info("");
        logger.info("Mapped Reads:                {}", statHolder.getMappedReadsCount());
        logger.info("");
        logger.info("Total Reads:                 {}", statHolder.getTotalReadsCount());
        logger.info("");
        if ( statHolder instanceof PacBioMovieStatHolder) {
            logger.info("Num Unique ZMWs:             {}", ((PacBioMovieStatHolder)statHolder).getNumUniqueZmws());
            logger.info("");
        }
        logger.info("----");
        logger.info("");
    }

    private void printStats() {
        // Print out our stats:

        logger.info("");
        logStats("OVERALL", overallReadStatHolder);

        logger.info("===========================================");

        for ( final Map.Entry<String, PacBioMovieStatHolder> entry : movieNameToStatHolderMap.entrySet() ) {
            logStats(entry.getKey(), entry.getValue());
        }
    }

    //==================================================================================================================
    // Helper Data Types:

}

