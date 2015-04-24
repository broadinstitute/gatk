package org.broadinstitute.hellbender.tools.picard.illumina.parser.readers;

import htsjdk.samtools.util.CloseableIterator;
import org.broadinstitute.hellbender.utils.text.parsers.BasicInputParser;

import java.io.File;

/**
 * Reads a single barcode file line by line and returns the barcode if there was a match or NULL otherwise.
 *
 * Barcode.txt file Format (consists of tab delimited columns, 1 record per row)
 * sequence_read    Matched(Y/N)    BarcodeSequenceMatched
 *
 * sequence read          - the actual bases at barcode position
 * Matched(y/n)           - Y or N indicating if there was a barcode match
 * BarcodeSequenceMatched - matched barcode sequence (empty if read did not match one of the barcodes).
 */
public final class BarcodeFileReader implements CloseableIterator<String> {
    private static final int Y_OR_N_COLUMN = 1;
    private static final int BARCODE_COLUMN = 2;
    private final BasicInputParser textIterator;

    public BarcodeFileReader(final File barcodeFile) {
        this.textIterator = new BasicInputParser(false, barcodeFile);
    }

    @Override
    public String next() {
        final String[] fields = textIterator.next();
        final String barcode;
        if (fields[Y_OR_N_COLUMN].equals("Y")) {
            barcode = fields[BARCODE_COLUMN];
        } else {
            barcode = null;
        }

        return barcode;
    }

    @Override
    public boolean hasNext() {
        return textIterator.hasNext();
    }

    public void remove() {
        throw new UnsupportedOperationException("Remove is not supported by " + BarcodeFileReader.class.getName());
    }

    public void close() {
        textIterator.close();
    }
}
