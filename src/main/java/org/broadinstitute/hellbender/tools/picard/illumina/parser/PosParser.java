package org.broadinstitute.hellbender.tools.picard.illumina.parser;

import htsjdk.samtools.util.CloseableIterator;

import org.broadinstitute.hellbender.tools.picard.illumina.parser.readers.AbstractIlluminaPositionFileReader;
import org.broadinstitute.hellbender.tools.picard.illumina.parser.readers.ClocsFileReader;
import org.broadinstitute.hellbender.tools.picard.illumina.parser.readers.LocsFileReader;
import org.broadinstitute.hellbender.tools.picard.illumina.parser.readers.PosFileReader;

import java.io.File;
import java.util.*;

import static htsjdk.samtools.util.CollectionUtil.makeSet;
import static java.util.Collections.unmodifiableSet;
import static org.broadinstitute.hellbender.tools.picard.illumina.parser.IlluminaDataType.Position;

/**
 * PosParser parses multiple files formatted as one of the three file formats that contain position information
 * only (pos, locs, and clocs).  This parser takes a map from tilesToFiles and a FileType enum value indicating
 * whether or not these are POS,LOCS, or CLOCS files.  The only client classes to this class should be IlluminaDataProvider
 * and test classes.  Check out AbstractIlluminaFileReader, PosFileReader, LocsFileReader, and ClocsFileReader for
 * more information on Position related illumina files.
 */
public class PosParser extends PerTileParser<PositionalData> {
    private static Set<IlluminaDataType> supportedTypes = unmodifiableSet(makeSet(Position));

    /**
     * The FileType of the files we are parsing
     */
    private final IlluminaFileUtil.SupportedIlluminaFormat fileType;

    public PosParser(final IlluminaFileMap tilesToFiles, final IlluminaFileUtil.SupportedIlluminaFormat fileType) {
        super(tilesToFiles);
        this.fileType = fileType;
    }

    public PosParser(final IlluminaFileMap tilesToFiles, final int startingTile, final IlluminaFileUtil.SupportedIlluminaFormat fileType) {
        super(tilesToFiles, startingTile);
        this.fileType = fileType;
    }

    /**
     * Make an CloseableIterator<PositionalData> based on the given file and fileType specified at construction.
     * This method wraps a reader in an iterator that converts it's output to the output format expected by
     * IlluminaDataProvider (PositionalData).
     *
     * @param file A file for the current tile being parsed
     * @return An iterator over the PositionalData in that file.
     */
    @Override
    protected CloseableIterator<PositionalData> makeTileIterator(final File file) {

        final AbstractIlluminaPositionFileReader fileReader;
        switch (fileType) {
            case Pos:
                fileReader = new PosFileReader(file);
                break;

            case Locs:
                fileReader = new LocsFileReader(file);
                break;

            case Clocs:
                fileReader = new ClocsFileReader(file);
                break;

            default:
                throw new IlluminaParserException("Unrecognized pos file type " + fileType.name());
        }

        return new CloseableIterator<PositionalData>() {
            private AbstractIlluminaPositionFileReader reader = fileReader;

            public void close() {
                reader.close();
            }

            public boolean hasNext() {
                return reader.hasNext();
            }

            public PositionalData next() {
                final AbstractIlluminaPositionFileReader.PositionInfo nextValue = reader.next();
                return new PositionalData() {
                    public int getXCoordinate() {
                        return nextValue.xQseqCoord;
                    }

                    public int getYCoordinate() {
                        return nextValue.yQseqCoord;
                    }
                };
            }

            public void remove() {
                throw new UnsupportedOperationException();
            }
        };
    }
}
