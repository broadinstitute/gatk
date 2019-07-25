package org.broadinstitute.hellbender.utils.codecs.gtf;

import com.sun.tools.javah.Gen;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.LocationAware;
import htsjdk.tribble.AbstractFeatureCodec;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodecHeader;
import htsjdk.tribble.SimpleFeature;
import htsjdk.tribble.annotation.Strand;
import htsjdk.tribble.readers.*;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.codecs.gencode.GencodeGtfFeature;
import scala.Int;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class GtfCodec extends AbstractFeatureCodec<GtfFeature, LineIterator> {

    static final Logger logger = LogManager.getLogger(GtfCodec.class);
    private static final String FIELD_DELIMITER = "\t";

    public static final int NO_FEATURE_ORDER = -1;
    public static final int NO_EXON_NUMBER = -1;

    private static final int NUM_FIELDS = 9;

    private static final int CHROMOSOME_NAME_INDEX = 0;
    private static final int ANNOTATION_SOURCE_INDEX = 1;
    private static final int FEATURE_TYPE_INDEX = 2;
    private static final int START_LOCATION_INDEX = 3;
    private static final int END_LOCATION_INDEX = 4;
    private static final int GENOMIC_STRAND_INDEX = 6;
    private static final int GENOMIC_PHASE_INDEX = 7;
    private static final int EXTRA_FIELDS_INDEX = 8;

    private static final String GTF_FILE_EXTENSION = "gtf";

    private static final String COMMENT_START = "#";

    private int currentLineNum = 0;

    public GtfCodec() {
        super(GtfFeature.class);
    }

    @Override
    public GtfFeature decode(final LineIterator lineIterator) {
        final String line = lineIterator.next();
        currentLineNum++;

        if (line.startsWith(COMMENT_START)) {
            return null;
        }

        final String[] splitLine = line.split(FIELD_DELIMITER, -1);

        if (splitLine.length != NUM_FIELDS) {
            throw new UserException.MalformedFile("Found an invalid number of columns in the given GENCODE file on line "
                    + currentLineNum + " - Given: " + splitLine.length + " Expected: " + NUM_FIELDS + " : " + line);
        }

        try {
            final String contig = splitLine[CHROMOSOME_NAME_INDEX];
            final String source = splitLine[ANNOTATION_SOURCE_INDEX];
            final String type = splitLine[FEATURE_TYPE_INDEX];
            final int start = Integer.valueOf(splitLine[START_LOCATION_INDEX]);
            final int end = Integer.valueOf(splitLine[END_LOCATION_INDEX]);
            final int phase = splitLine[GENOMIC_PHASE_INDEX].equals(".")? -1 : Integer.valueOf(splitLine[GENOMIC_PHASE_INDEX]);
            final Strand strand = Strand.decode(splitLine[GENOMIC_STRAND_INDEX]);
            final Map<String, String> attributes = parseAttributes(splitLine[EXTRA_FIELDS_INDEX]);
            return new GtfFeature(contig, source, type, start, end, strand, phase, attributes);
        } catch (final NumberFormatException ex ) {
            throw new UserException.MalformedFile("Cannot read integer value for start/end position!");
        }
    }

    public static Map<String,String> parseAttributes(final String attributesString) {
        final Map<String, String> attributes = new HashMap<>();

        int currentSubstringStart=0;
        while (currentSubstringStart < attributesString.length()) {
            final int keyStart = currentSubstringStart;
            final int keyEnd = attributesString.indexOf(" ");
            final String key = attributesString.substring(keyStart, keyEnd).trim();

            final int valueStart = attributesString.indexOf("\"", keyEnd);
            final int valueEnd = attributesString.indexOf("\"", valueStart + 1);
            final String value = attributesString.substring(valueStart+1, valueEnd).trim();
            attributes.put(key, value);

            currentSubstringStart = valueEnd + 1;

        }
        return attributes;
    }

    @Override
    public Feature decodeLoc(LineIterator lineIterator) {
        final String line = lineIterator.next();
        if (line.startsWith(COMMENT_START)) {
            return null;
        }

        final String[] splitLine = line.split(FIELD_DELIMITER, -1);

        try {
            return new SimpleFeature(splitLine[CHROMOSOME_NAME_INDEX], Integer.valueOf(splitLine[START_LOCATION_INDEX]), Integer.valueOf(splitLine[END_LOCATION_INDEX]));
        } catch (final NumberFormatException ex ) {
            throw new UserException.MalformedFile("Cannot read integer value for start/end position!");
        }
    }

    @Override
    public boolean canDecode(final String inputFilePath) {
        boolean canDecode;
        try {
            // Simple file and name checks to start with:
            Path p = IOUtil.getPath(inputFilePath);

            canDecode = p.getFileName().toString().toLowerCase().endsWith("." + GTF_FILE_EXTENSION);

            if (canDecode) {

                // Crack open the file and look at the top of it:
                try ( BufferedReader br = new BufferedReader(new InputStreamReader(Files.newInputStream(p))) ) {

                    String line = br.readLine();
                    while (line.startsWith("#")) {
                        line = br.readLine();
                        if ( line == null ) {
                            return false;
                        }
                    }

                    // make sure line conforms to gtf spec
                    final String[] fields = line.split("\t");

                    canDecode &= fields.length == 9;

                    if (canDecode) {
                        // check that start and end fields are integers
                        try {
                            final int start = Integer.parseInt(fields[3]);
                            final int end = Integer.parseInt(fields[4]);
                        } catch (NumberFormatException | NullPointerException nfe) {
                            canDecode = false;
                        }

                        // check for strand

                        canDecode &= fields[6].equals("+") || fields[6].equals("-") || fields[6].equals(".");
                    }
                }

            }
        }
        catch (final FileNotFoundException ex) {
            logger.warn("File does not exist! - " + inputFilePath + " - returning can decode as failure.");
            canDecode = false;
        }
        catch (final IOException ex) {
            logger.warn("Caught IOException on file: " + inputFilePath + " - returning can decode as failure.");
            canDecode = false;
        }

        return canDecode;
    }

    @Override
    public FeatureCodecHeader readHeader(LineIterator lineIterator) throws IOException {

        List<String> header = new ArrayList<>();
        while(lineIterator.hasNext()) {
            String line = lineIterator.peek();
            if (line.startsWith(COMMENT_START)) {
                header.add(line);
                lineIterator.next();
            } else {
                break;
            }
        }

        return new FeatureCodecHeader(header, FeatureCodecHeader.NO_HEADER_END);
    }

    @Override
    public LineIterator makeSourceFromStream(final InputStream bufferedInputStream) {
        return new LineIteratorImpl(new SynchronousLineReader(bufferedInputStream));
    }

    @Override
    public LocationAware makeIndexableSourceFromStream(final InputStream bufferedInputStream) {
        return new AsciiLineReaderIterator(AsciiLineReader.from(bufferedInputStream));
    }

    @Override
    public boolean isDone(final LineIterator lineIterator) {
        return !lineIterator.hasNext();
    }

    @Override
    public void close(final LineIterator lineIterator) {
        CloserUtil.close(lineIterator);
    }
}