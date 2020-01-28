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
import java.util.stream.Collectors;

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

    private static enum MapStatus {
        UNMAPPED("UNMAPPED"),
        UNPLACED("UNPLACED"),
        MAPPED("MAPPED");

        private String name;

        MapStatus(final String name) {
            this.name = name;
        }

        public String getName() { return name; }
        public String toString() { return name; }
    }

    /**
     * Simple data class to hold counts and supplementary for stats that we need on a per-movie basis.
     * NOTE: Pac-Bio specific.
     */
    private static class ReadStatHolder {
        private long         totalReadsCount = 0;

        private long         unmappedReadsCount = 0;
        private long         unplacedReadsCount = 0;
        private long         mappedReadsCount = 0;

        private long         primaryReadsCount = 0;
        private long         secondaryReadsCount = 0;
        private long         supplementaryReadsCount = 0;

        private long adapterAnnotatedReadCount    = 0;
        // Post-correction:
        private long barcodeAnnotatedReadCount    = 0;
        // Pre-correction:
        private long rawBarcodeAnnotatedReadCount = 0;
        private long umiAnnotatedReadCount        = 0;

        private HashMap<MapStatus, Integer> mapStatusPrimaryCountMap = new HashMap<>();
        private HashMap<MapStatus, Integer> mapStatusSecondaryCountMap = new HashMap<>();
        private HashMap<MapStatus, Integer> mapStatusSupplementaryCountMap = new HashMap<>();

        ReadStatHolder() {
            Arrays.stream(MapStatus.values()).forEach(
                    mapStatus -> {
                        mapStatusPrimaryCountMap.put(mapStatus, 0);
                        mapStatusSecondaryCountMap.put(mapStatus, 0);
                        mapStatusSupplementaryCountMap.put(mapStatus, 0);
                    }
            );
        }

        public long getAdapterAnnotatedReadCount() {
            return adapterAnnotatedReadCount;
        }

        public long getBarcodeAnnotatedReadCount() {
            return barcodeAnnotatedReadCount;
        }

        public long getRawBarcodeAnnotatedReadCount() {
            return rawBarcodeAnnotatedReadCount;
        }

        public long getUmiAnnotatedReadCount() {
            return umiAnnotatedReadCount;
        }

        long getTotalReadsCount() {
            return totalReadsCount;
        }

        long getUnmappedReadsCount() {
            return unmappedReadsCount;
        }
        long getUnplacedReadsCount() {
            return unplacedReadsCount;
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
            return unmappedReadsCount + "," + unplacedReadsCount + "," + secondaryReadsCount + "," + supplementaryReadsCount + "," +
                    primaryReadsCount + "," + mappedReadsCount + "," + totalReadsCount + "," +
                    adapterAnnotatedReadCount + "," +
                    umiAnnotatedReadCount + "," +
                    rawBarcodeAnnotatedReadCount + "," +
                    barcodeAnnotatedReadCount + "," +
            createMapStatusAlignmentCsvData();
        }

        String createMapStatusAlignmentCsvData() {
            final StringBuilder sb = new StringBuilder();
            for (final MapStatus status : MapStatus.values()) {
                sb.append(mapStatusPrimaryCountMap.get(status));
                sb.append(',');
            }
            for (final MapStatus status : MapStatus.values()) {
                sb.append(mapStatusSecondaryCountMap.get(status));
                sb.append(',');
            }
            for (final MapStatus status : MapStatus.values()) {
                sb.append(mapStatusSupplementaryCountMap.get(status));
                sb.append(',');
            }

            return sb.toString();
        }

        void incrementCountMap( final MapStatus mapStatus, final HashMap<MapStatus, Integer> mapStatusCountMap ) {
            final Integer count = mapStatusCountMap.get(mapStatus) + 1;
            mapStatusCountMap.put(mapStatus, count);
        }

        void accountForTranscriptome10xFlags(final GATKRead read ) {
            if ( read.getAttributeAsString("ZA") != null ) {
                ++adapterAnnotatedReadCount;
            }
            if ( read.getAttributeAsString("CB") != null ) {
                ++barcodeAnnotatedReadCount;
            }
            if ( read.getAttributeAsString("CR") != null ) {
                ++rawBarcodeAnnotatedReadCount;
            }
            if ( read.getAttributeAsString("ZU") != null ) {
                ++umiAnnotatedReadCount;
            }
        }

        void accountForMappingFlags(final GATKRead read) {
            ++totalReadsCount;

            // Mapping:
            if ( read.isUnmapped() ) {
                ++unmappedReadsCount;
            }
            else if ( read.isUnplaced() ) {
                ++unplacedReadsCount;
            }
            else {
                ++mappedReadsCount;
            }

            // Alignment:
            if (read.isSecondaryAlignment() &&  !read.isSupplementaryAlignment()) {
                ++secondaryReadsCount;

                if ( read.isUnmapped() ) {
                    incrementCountMap(MapStatus.UNMAPPED, mapStatusSecondaryCountMap);
                }
                else if ( read.isUnplaced() ) {
                    incrementCountMap(MapStatus.UNPLACED, mapStatusSecondaryCountMap);
                }
                else {
                    incrementCountMap(MapStatus.MAPPED, mapStatusSecondaryCountMap);
                }
            }
            else if (!read.isSecondaryAlignment() &&  read.isSupplementaryAlignment()) {
                ++supplementaryReadsCount;

                if ( read.isUnmapped() ) {
                    incrementCountMap(MapStatus.UNMAPPED, mapStatusSupplementaryCountMap);
                }
                else if ( read.isUnplaced() ) {
                    incrementCountMap(MapStatus.UNPLACED, mapStatusSupplementaryCountMap);
                }
                else {
                    incrementCountMap(MapStatus.MAPPED, mapStatusSupplementaryCountMap);
                }
            }
            else {
                ++primaryReadsCount;

                if ( read.isUnmapped() ) {
                    incrementCountMap(MapStatus.UNMAPPED, mapStatusPrimaryCountMap);
                }
                else if ( read.isUnplaced() ) {
                    incrementCountMap(MapStatus.UNPLACED, mapStatusPrimaryCountMap);
                }
                else {
                    incrementCountMap(MapStatus.MAPPED, mapStatusPrimaryCountMap);
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
            return super.toCsvRow() + zmwSet.size();
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
            movieNameToStatHolderMap.get(movieName).accountForTranscriptome10xFlags(read);
        }

        overallReadStatHolder.accountForMappingFlags(read);
        overallReadStatHolder.accountForTranscriptome10xFlags(read);
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
        sb.append("Unplaced Reads"); sb.append(',');
        sb.append("Secondary Aligned Reads"); sb.append(',');
        sb.append("Supplementary Aligned Reads"); sb.append(',');
        sb.append("Primary Aligned Reads"); sb.append(',');
        sb.append("Mapped Reads"); sb.append(',');
        sb.append("Total Reads"); sb.append(',');
        sb.append("Adapter Annotated Reads"); sb.append(',');
        sb.append("UMI Annotated Reads"); sb.append(',');
        sb.append("Raw Barcode Annotated Reads"); sb.append(',');
        sb.append("Barcode Annotated Reads"); sb.append(',');

        for ( final String alignmentType : Arrays.asList("Primary Aligned", "Secondary Aligned", "Supplementary Aligned")) {
            for ( final MapStatus status : MapStatus.values() ) {
                sb.append(status.name);
                sb.append(' ');
                sb.append(alignmentType);
                sb.append(',');
            }
        }

        sb.append("Num Unique ZMWs"); sb.append(',');

        sb.append('\n');

        // Write overall stats first:
        sb.append("OVERALL");
        sb.append(',');
        sb.append(overallReadStatHolder.toCsvRow());
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

    private void logMapAlignmentMatrix(final ReadStatHolder statHolder) {
        final String header = Arrays.stream(MapStatus.values()).map(Object::toString).collect(Collectors.joining("\t"));
        logger.info("                 " + header);

        final String primary = Arrays.stream(MapStatus.values())
                .map(status -> statHolder.mapStatusPrimaryCountMap.get(status))
                .map(Object::toString).collect(Collectors.joining("\t"));
        logger.info("Primary          " + primary);

        final String secondary = Arrays.stream(MapStatus.values())
                .map(status -> statHolder.mapStatusSecondaryCountMap.get(status))
                .map(Object::toString).collect(Collectors.joining("\t"));
        logger.info("Secondary        " + secondary);

        final String supplementary = Arrays.stream(MapStatus.values())
                .map(status -> statHolder.mapStatusSupplementaryCountMap.get(status))
                .map(Object::toString).collect(Collectors.joining("\t"));
        logger.info("Supplementary    " + supplementary);
    }

    private void logStats(final String name, final ReadStatHolder statHolder) {
        logger.info("----");
        logger.info("");
        logger.info("{}", name);
        logger.info("");

        logger.info("Unmapped Reads:              {}", statHolder.getUnmappedReadsCount());
        logger.info("Unplaced Reads:              {}", statHolder.getUnplacedReadsCount());
        logger.info("Secondary Aligned Reads:     {}", statHolder.getSecondaryReadsCount());
        logger.info("Supplementary Aligned Reads: {}", statHolder.getSupplementaryReadsCount());
        logger.info("Primary Aligned Reads:       {}", statHolder.getPrimaryReadsCount());
        logger.info("");
        logger.info("Mapped Reads:                {}", statHolder.getMappedReadsCount());
        logger.info("");
        logger.info("Total Reads:                 {}", statHolder.getTotalReadsCount());
        logger.info("");
        logger.info("Adapter Annotated Reads:     {}", statHolder.getAdapterAnnotatedReadCount());
        logger.info("UMI Annotated Reads:         {}", statHolder.getUmiAnnotatedReadCount());
        logger.info("Raw Barcode Annotated Reads: {}", statHolder.getRawBarcodeAnnotatedReadCount());
        logger.info("Barcode Annotated Reads:     {}", statHolder.getBarcodeAnnotatedReadCount());
        logger.info("");
        logMapAlignmentMatrix(statHolder);
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

