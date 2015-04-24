package org.broadinstitute.hellbender.tools.picard.illumina.parser.readers;

import java.io.File;
import java.nio.ByteBuffer;

import static java.lang.String.format;
import static org.broadinstitute.hellbender.tools.picard.illumina.parser.readers.MMapBackedIteratorFactory.getLongIterator;

/**
 * Annoyingly, there are two different files with extension .bci in NextSeq output.  This reader handles
 * the file that contains virtual file pointers into a .bcl.bgzf file.  After the header, there is a 64-bit record
 * per tile.
 */
public final class BclIndexReader {
    private static final int BCI_HEADER_SIZE = 8;
    private static final int BCI_VERSION = 0;

    private final BinaryFileIterator<Long> bciIterator;
    private final int numTiles;
    private final File bciFile;
    private int nextRecordNumber = 0;

    public BclIndexReader(final File bclFile) {
        bciFile = new File(bclFile.getAbsolutePath() + ".bci");
        bciIterator = getLongIterator(BCI_HEADER_SIZE, bciFile);
        final ByteBuffer headerBytes = bciIterator.getHeaderBytes();
        final int actualVersion = headerBytes.getInt();
        if (actualVersion != BCI_VERSION) {
            throw new IlluminaReaderException(format("Unexpected version number %d in %s", actualVersion, bciFile.getAbsolutePath()));
        }
        numTiles = headerBytes.getInt();
    }

    public int getNumTiles() {
        return numTiles;
    }

    public long get(final int recordNumber) {
        if (recordNumber < nextRecordNumber) {
            throw new IllegalArgumentException("Can only read forward");
        }
        if (recordNumber > nextRecordNumber) {
            bciIterator.skipElements(recordNumber - nextRecordNumber);
            nextRecordNumber = recordNumber;
        }
        ++nextRecordNumber;
        return bciIterator.getElement();
    }

    public File getBciFile() {
        return bciFile;
    }
}
