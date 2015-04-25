package org.broadinstitute.hellbender.tools.picard.illumina.parser;

import htsjdk.samtools.util.CloseableIterator;
import org.broadinstitute.hellbender.tools.picard.illumina.parser.readers.BarcodeFileReader;

import java.io.File;
import java.util.Set;

import static htsjdk.samtools.util.CollectionUtil.makeSet;
import static java.util.Collections.unmodifiableSet;
import static org.broadinstitute.hellbender.tools.picard.illumina.parser.IlluminaDataType.Barcodes;

/**
 * @author jburke@broadinstitute.org
 */
class BarcodeParser extends PerTileParser<BarcodeData> {

    private static final Set<IlluminaDataType> SUPPORTED_TYPES = unmodifiableSet(makeSet(Barcodes));

    public BarcodeParser(final IlluminaFileMap tilesToFiles) {
        super(tilesToFiles);
    }

    @Override
    protected CloseableIterator<BarcodeData> makeTileIterator(File nextTileFile) {
        return new BarcodeDataIterator(nextTileFile);
    }

    public Set<IlluminaDataType> supportedTypes() {
        return SUPPORTED_TYPES;
    }

    private static class BarcodeDataIterator implements CloseableIterator<BarcodeData> {
        private BarcodeFileReader bfr;

        public BarcodeDataIterator(final File file) {
            bfr = new BarcodeFileReader(file);
        }

        public void close() {
            bfr.close();
        }

        public boolean hasNext() {
            return bfr.hasNext();
        }

        public BarcodeData next() {
            return new BarcodeData() {
                public String getBarcode() {
                    return bfr.next();
                }
            };
        }

        public void remove() {
            throw new UnsupportedOperationException();
        }
    }
}
