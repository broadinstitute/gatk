package org.broadinstitute.hellbender.tools.picard.illumina.parser.readers;

import htsjdk.samtools.util.CloserUtil;
import org.broadinstitute.hellbender.utils.text.parsers.BasicInputParser;

import java.io.File;

import static java.lang.Float.parseFloat;

/**
 * The pos file format is one 3 Illumina formats(pos, locs, and clocs) that stores position data exclusively.
 * pos files store position data for successive clusters in tabbed delimited coordinated pairs, 1 per file row e.g.:
 * <p>
 * xPos1\tyPos1
 * xPos2\tyPos2
 * 102.0\t303.3
 * ...
 * xPosn-1\yPosn-1
 * xPosn\yPosn
 * <p>
 * Where n = the total number of clusters (and therefore lines) in the file.
 */
public class PosFileReader extends AbstractIlluminaPositionFileReader {

    private final BasicInputParser parser;

    public PosFileReader(final File posFile) {
        super(posFile);
        this.parser = new BasicInputParser(true, posFile);
    }

    /**
     * Read a line of text and parse it into two float values, create a PositionInfo and return it
     */
    @Override
    protected PositionInfo unsafeNextInfo() {
        final String[] strVals = this.parser.next();
        if (strVals.length != 2) {
            throw new IlluminaReaderException("Pos file number of values != 2, found (" + strVals.length + ")" + makeExceptionMsg());
        }
        try {
            final float xVal = parseFloat(strVals[0]);
            final float yVal = parseFloat(strVals[1]);

            if (xVal < 0 || yVal < 0) {
                throw new NumberFormatException("X and Y pos values cannot be negative!");
            }

            return new PositionInfo(xVal, yVal, getLane(), getTile());
        } catch (final NumberFormatException nfe) {
            throw new IlluminaReaderException("Bad x or y value in " + makeExceptionMsg(), nfe);
        }
    }

    @Override
    protected String makeExceptionMsg() {
        return "pos file( " + parser.getFileName() +
                " ) on line number( " + parser.getCurrentLineNumber() +
                " ) with current line = " + parser.getCurrentLine();
    }

    public boolean hasNext() {
        return parser.hasNext();
    }

    public void close() {
        CloserUtil.close(parser);
    }
}
