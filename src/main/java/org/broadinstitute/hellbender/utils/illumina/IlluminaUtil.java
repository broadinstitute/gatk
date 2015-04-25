package org.broadinstitute.hellbender.utils.illumina;

import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.StringUtil;

import java.util.List;

/**
 * Misc utilities for working with Illumina specific files and data
 *
 * @author jburke@broadinstitute.org
 */
public class IlluminaUtil {

    public static final String BARCODE_DELIMITER = "-";

    /**
     * Parse the tile # from the read name.
     * If we find that there are other elements needed from the read name, it might be a good idea to put
     * makeReadName() and various get..() methods into a new class.
     *
     * @param readName As produced by IlluminaUtil.makeReadName()
     * @return tile number, or null if read name is not in correct format.
     */
    public static Integer getTileFromReadName(final String readName) {
        final int first = readName.indexOf(':');
        if (first > 0) {
            final int second = readName.indexOf(':', first+1);
            if (second > 0) {
                final int third = readName.indexOf(':', second+1);
                if (third > 0) {
                    return Integer.parseInt(readName.substring(second + 1, third));
                }
            }
        }

        return null;
    }

    // Strings indented below to make these easier to compare visually.
    /** Describes adapters used on each pair of strands */
    public static enum IlluminaAdapterPair implements AdapterPair {

        PAIRED_END("AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT",  //58 bases
                   "AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG"), // 61 bases

        INDEXED ("AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT",
                 "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG"), // note  8 N's  // 67 bases

        SINGLE_END ("AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT",
                    "AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG"),

        NEXTERA_V1("AATGATACGGCGACCACCGAGATCTACACGCCTCCCTCGCGCCATCAGAGATGTGTATAAGAGACAG",
                   "CTGTCTCTTATACACATCTCTGAGCGGGCTGGCAAGGCAGACCGNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG"),

        NEXTERA_V2("AATGATACGGCGACCACCGAGATCTACACNNNNNNNNTCGTCGGCAGCGTCAGATGTGTATAAGAGACAG",
                   "CTGTCTCTTATACACATCTCCGAGCCCACGAGACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG"),

        DUAL_INDEXED("AATGATACGGCGACCACCGAGATCTNNNNNNNNACACTCTTTCCCTACACGACGCTCTTCCGATCT",
                    "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG"),

        FLUIDIGM("AATGATACGGCGACCACCGAGATCTACACTGACGACATGGTTCTACA",
                 "AGACCAAGTCTCTGCTACCGTANNNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG"),

        TRUSEQ_SMALLRNA("AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGACGATC",
                        "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG"),

        // This one is at the end of the list because its 3' is a subset of several of the 3's above.
        // There are unit tests that try all AdapterPairs, and this one should go at the end os
        // it is checked last.
        ALTERNATIVE_SINGLE_END("AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGACGATC",
                               "TCGTATGCCGTCTTCTGCTTG");


        final String fivePrime, threePrime, fivePrimeReadOrder;
        final byte[]  fivePrimeBytes, threePrimeBytes, fivePrimeReadOrderBytes;


        private IlluminaAdapterPair(final String fivePrime, final String threePrime) {
            this.threePrime = threePrime;
            this.threePrimeBytes = StringUtil.stringToBytes(threePrime);

            this.fivePrime = fivePrime;
            this.fivePrimeReadOrder = SequenceUtil.reverseComplement(fivePrime);
            this.fivePrimeBytes = StringUtil.stringToBytes(fivePrime);
            this.fivePrimeReadOrderBytes = StringUtil.stringToBytes(fivePrimeReadOrder);
        }

        public String get3PrimeAdapter(){ return threePrime; }
        public String get5PrimeAdapter(){ return fivePrime; }
        public String get3PrimeAdapterInReadOrder(){ return threePrime; }
        public String get5PrimeAdapterInReadOrder() { return fivePrimeReadOrder; }
        public byte[] get3PrimeAdapterBytes() { return threePrimeBytes; }
        public byte[] get5PrimeAdapterBytes() { return fivePrimeBytes; }
        public byte[] get3PrimeAdapterBytesInReadOrder() { return threePrimeBytes; }
        public byte[] get5PrimeAdapterBytesInReadOrder()  { return fivePrimeReadOrderBytes; }
        public String getName() { return this.name(); }
    }

    /**
     * Concatenates all the barcode sequences with BARCODE_DELIMITER
     * @param barcodes
     * @return A single string representation of all the barcodes
     */
    public static String barcodeSeqsToString(final List<String> barcodes) {
        return barcodeSeqsToString(barcodes.toArray(new String[barcodes.size()]));
    }

    /**
     * Concatenates all the barcode sequences with BARCODE_DELIMITER
     * @param barcodes
     * @return A single string representation of all the barcodes
     */
    public static String barcodeSeqsToString(final String barcodes[]) {
        final StringBuilder sb = new StringBuilder();
        for (final String bc : barcodes) {
            if (sb.length() > 0) sb.append(BARCODE_DELIMITER);
            sb.append(bc);
        }
        return sb.toString();
    }

    /**
     * Concatenates all the barcode sequences with BARCODE_DELIMITER
     * @param barcodes
     * @return A single string representation of all the barcodes
     */
    public static String barcodeSeqsToString(final byte barcodes[][]) {
        final String bcs[] = new String[barcodes.length];
        for (int i = 0; i < barcodes.length; i++) {
            bcs[i] = StringUtil.bytesToString(barcodes[i]);
        }
        return barcodeSeqsToString(bcs);
    }
}
