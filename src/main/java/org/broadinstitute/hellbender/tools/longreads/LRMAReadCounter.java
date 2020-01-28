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

    //==================================================================================================================
    // Public Members:
    @Argument(
            fullName  = "output-csv-file",
            optional = true,
            doc = "Output the results to the given CSV file name.")
    private String outputCsvFileName = "";

    @Argument(
            fullName  = "use-file-as-row-header",
            optional = true,
            doc = "Use the file name as the row header for overall counts.")
    private Boolean useFileAsRowHeader = false;

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
            movieNameToStatHolderMap.get(movieName).countReadStats(read);
        }

        overallReadStatHolder.countReadStats(read);
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

        for ( final Transcriptome10xAttribute transAttr : Transcriptome10xAttribute.values() ) {
            for ( final MapStatus status : MapStatus.values() ) {
                for ( final ReadPriorityStatus readPriorityStatus : ReadPriorityStatus.values() ) {
                    sb.append(transAttr.name);
                    sb.append(' ');
                    sb.append(status.name);
                    sb.append(' ');
                    sb.append(readPriorityStatus.name);
                    sb.append(',');
                }
            }
        }

        sb.append("Num Unique ZMWs"); sb.append(',');
        sb.append('\n');

        // Write overall stats first:
        if ( useFileAsRowHeader ){
            sb.append(readArguments.getReadFiles()
                    .stream()
                    .map( x -> x.toURI().toString())
                    .collect(Collectors.joining(", ")));
        }
        else {
            sb.append("OVERALL");
        }
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
                .map(status -> statHolder.mapStatusReadPriorityStatusCountMap.get(status).get(ReadPriorityStatus.PRIMARY))
                .map(Object::toString).collect(Collectors.joining("\t"));
        logger.info("Primary          " + primary);

        final String secondary = Arrays.stream(MapStatus.values())
                .map(status -> statHolder.mapStatusReadPriorityStatusCountMap.get(status).get(ReadPriorityStatus.SECONDARY))
                .map(Object::toString).collect(Collectors.joining("\t"));
        logger.info("Secondary        " + secondary);

        final String supplementary = Arrays.stream(MapStatus.values())
                .map(status -> statHolder.mapStatusReadPriorityStatusCountMap.get(status).get(ReadPriorityStatus.SUPPLEMENTARY))
                .map(Object::toString).collect(Collectors.joining("\t"));
        logger.info("Supplementary    " + supplementary);
    }

    private void logTranscriptTagMapAlignmentMatrices(final ReadStatHolder statHolder) {
        final StringBuilder headerPadding = new StringBuilder();

        final int firstColWidth = 17;
        for (int i = 0; i < firstColWidth; ++i) {
            headerPadding.append(' ');
        }
        final String header = headerPadding.toString() + Arrays.stream(MapStatus.values()).map(Object::toString).collect(Collectors.joining("\t"));

//        final int colPad = 4;
//        final int columnWidth[] = {17,
//                MapStatus.values()[0].name.length() + colPad,
//                MapStatus.values()[1].name.length() + colPad,
//                MapStatus.values()[2].name.length() + colPad};
//
//        final String rowFormatString = "%" + columnWidth[0] + "s%" +
//                "%" + columnWidth[1] + "s%" +
//                "%" + columnWidth[2] + "s%" +
//                "%" + columnWidth[3] + "s%";
//        final String header = String.format(rowFormatString, "", MapStatus.values()[0], MapStatus.values()[1], MapStatus.values()[2]);

        for ( final Transcriptome10xAttribute attr : Transcriptome10xAttribute.values() ) {
            logger.info("--------------------------");
            logger.info("Reads containing " + attr.name() + " (" + attr.value + ")");
            logger.info("");
            logger.info(header);

            for ( final ReadPriorityStatus rps : ReadPriorityStatus.values()) {
                final String row = Arrays.stream(MapStatus.values())
                        .map(status -> statHolder.transcriptAnnotationMapStatusReadPriorityCountMap.get(attr).get(status).get(rps))
                        .map(Object::toString).collect(Collectors.joining("\t"));

                final StringBuilder padding = new StringBuilder();
                for ( int i = 0 ; i < firstColWidth - rps.name.length(); ++i ){
                    padding.append(' ');
                }

                logger.info(rps.name + padding.toString() + row);
            }
            logger.info("");
        }

        logger.info("--------------------------");

        for ( final Transcriptome10xAttribute transAttr : Transcriptome10xAttribute.values() ) {
            for ( final MapStatus status : MapStatus.values() ) {
                for ( final ReadPriorityStatus readPriorityStatus : ReadPriorityStatus.values() ) {
                    logger.info("Reads containing " + transAttr.name() + " (" + transAttr.value + "): " +
                            status.name() + "/" + readPriorityStatus.name() + "\t-\t" +
                            statHolder.transcriptAnnotationMapStatusReadPriorityCountMap.get(transAttr).get(status).get(readPriorityStatus));
                }
            }
        }
        logger.info("--------------------------");


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
        logTranscriptTagMapAlignmentMatrices(statHolder);
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

    private enum MapStatus {
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

    private enum ReadPriorityStatus {
        PRIMARY("PRIMARY"),
        SECONDARY("SECONDARY"),
        SUPPLEMENTARY("SUPPLEMENTARY");

        private String name;

        ReadPriorityStatus(final String name) {
            this.name = name;
        }

        public String getName() { return name; }
        public String toString() { return name; }
    }

    private enum Transcriptome10xAttribute {
        ADAPTER("ADAPTER", "ZA"),
        UMI("UMI", "ZU"),
        RAW_BARCODE("RAW_BARCODE", "CR"),
        BARCODE("BARCODE", "CB");

        private String name;
        private String value;

        Transcriptome10xAttribute(final String name, final String value) {
            this.name = name;
            this.value = value;
        }

        public String getName() { return name; }
        public String getValue() { return value; }
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

        private HashMap<MapStatus, HashMap<ReadPriorityStatus, Integer>> mapStatusReadPriorityStatusCountMap = new HashMap<>();
        private HashMap<Transcriptome10xAttribute, HashMap<MapStatus, HashMap<ReadPriorityStatus, Integer>>> transcriptAnnotationMapStatusReadPriorityCountMap = new HashMap<>();

        ReadStatHolder() {
            Arrays.stream(MapStatus.values()).forEach(
                    mapStatus -> {

                        final HashMap<ReadPriorityStatus, Integer> map = new HashMap<>();
                        Arrays.stream(ReadPriorityStatus.values()).forEach(
                                readPriorityStatus -> map.put(readPriorityStatus, 0)
                        );
                        mapStatusReadPriorityStatusCountMap.put(mapStatus, map);
                    }
            );

            Arrays.stream(Transcriptome10xAttribute.values()).forEach(
                    attr -> {
                        final HashMap<MapStatus, HashMap<ReadPriorityStatus, Integer>> msRpCountMap = new HashMap<>();
                        Arrays.stream(MapStatus.values()).forEach(
                                mapStatus -> {

                                    final HashMap<ReadPriorityStatus, Integer> rpCountMap = new HashMap<>();
                                    Arrays.stream(ReadPriorityStatus.values()).forEach(
                                            readPriorityStatus -> rpCountMap.put(readPriorityStatus, 0)
                                    );
                                    msRpCountMap.put(mapStatus, rpCountMap);
                                }
                        );
                        transcriptAnnotationMapStatusReadPriorityCountMap.put(attr, msRpCountMap);
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
                    createMapStatusAlignmentCsvData() +
                    createTranscriptomeMapStatusAlignmentCsvData();
        }

        String createMapStatusAlignmentCsvData() {
            final StringBuilder sb = new StringBuilder();

            for ( final ReadPriorityStatus readPriorityStatus : ReadPriorityStatus.values() ) {
                for ( final MapStatus status : MapStatus.values() ) {
                    sb.append(mapStatusReadPriorityStatusCountMap.get(status).get(readPriorityStatus));
                    sb.append(',');
                }
            }

            return sb.toString();
        }

        String createTranscriptomeMapStatusAlignmentCsvData() {
            final StringBuilder sb = new StringBuilder();

            for ( final Transcriptome10xAttribute transAttr : Transcriptome10xAttribute.values() ) {
                for ( final MapStatus status : MapStatus.values() ) {
                    for ( final ReadPriorityStatus readPriorityStatus : ReadPriorityStatus.values() ) {
                        sb.append(transcriptAnnotationMapStatusReadPriorityCountMap.get(transAttr).get(status).get(readPriorityStatus));
                        sb.append(',');
                    }
                }
            }

            return sb.toString();
        }

        void incrementCountMap( final MapStatus mapStatus, final ReadPriorityStatus readPriorityStatus ) {
            incrementCountMap(mapStatus, readPriorityStatus, mapStatusReadPriorityStatusCountMap);
        }

        void incrementCountMap( final MapStatus mapStatus, final ReadPriorityStatus readPriorityStatus,
                                final HashMap<MapStatus, HashMap<ReadPriorityStatus, Integer>> map) {
            final Integer count = map.get(mapStatus).get(readPriorityStatus) + 1;
            map.get(mapStatus).put(readPriorityStatus, count);
        }

        void accountForTranscriptome10xFlags(final GATKRead read ) {
            if ( read.getAttributeAsString(Transcriptome10xAttribute.ADAPTER.value) != null ) {
                ++adapterAnnotatedReadCount;
                incrementCountMap(getMapStatus(read), getReadPriorityStatus(read), transcriptAnnotationMapStatusReadPriorityCountMap.get(Transcriptome10xAttribute.ADAPTER));
            }
            if ( read.getAttributeAsString(Transcriptome10xAttribute.BARCODE.value) != null ) {
                ++barcodeAnnotatedReadCount;
                incrementCountMap(getMapStatus(read), getReadPriorityStatus(read), transcriptAnnotationMapStatusReadPriorityCountMap.get(Transcriptome10xAttribute.BARCODE));
            }
            if ( read.getAttributeAsString(Transcriptome10xAttribute.RAW_BARCODE.value) != null ) {
                ++rawBarcodeAnnotatedReadCount;
                incrementCountMap(getMapStatus(read), getReadPriorityStatus(read), transcriptAnnotationMapStatusReadPriorityCountMap.get(Transcriptome10xAttribute.RAW_BARCODE));
            }
            if ( read.getAttributeAsString(Transcriptome10xAttribute.UMI.value) != null ) {
                ++umiAnnotatedReadCount;
                incrementCountMap(getMapStatus(read), getReadPriorityStatus(read), transcriptAnnotationMapStatusReadPriorityCountMap.get(Transcriptome10xAttribute.UMI));
            }
        }

        static MapStatus getMapStatus(final GATKRead read) {
            if ( read.isUnmapped() ) {
                return MapStatus.UNMAPPED;
            }
            else if ( read.isUnplaced() ) {
                return MapStatus.UNPLACED;
            }
            else {
                return MapStatus.MAPPED;
            }
        }

        static ReadPriorityStatus getReadPriorityStatus(final GATKRead read) {
            if (read.isSecondaryAlignment() &&  !read.isSupplementaryAlignment()) {
                return ReadPriorityStatus.SECONDARY;
            }
            else if (!read.isSecondaryAlignment() &&  read.isSupplementaryAlignment()) {
                return ReadPriorityStatus.SUPPLEMENTARY;
            }
            else {
                return ReadPriorityStatus.PRIMARY;
            }
        }

        void countReadStats(final GATKRead read) {
            final MapStatus mapStatus = getMapStatus(read);
            final ReadPriorityStatus readPriorityStatus = getReadPriorityStatus(read);

            ++totalReadsCount;

            switch(mapStatus) {
                case MAPPED: ++mappedReadsCount; break;
                case UNMAPPED: ++unmappedReadsCount; break;
                case UNPLACED: ++unplacedReadsCount; break;
            }

            switch(readPriorityStatus) {
                case PRIMARY: ++primaryReadsCount; break;
                case SECONDARY: ++secondaryReadsCount; break;
                case SUPPLEMENTARY: ++supplementaryReadsCount; break;
            }

            incrementCountMap(mapStatus, readPriorityStatus);
            accountForTranscriptome10xFlags(read);
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

}

