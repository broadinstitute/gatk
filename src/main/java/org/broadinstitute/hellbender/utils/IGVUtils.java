package org.broadinstitute.hellbender.utils;

import htsjdk.samtools.util.Locatable;

import java.io.PrintStream;

/**
 * Utilities for interacting with IGV-specific formats. See https://software.broadinstitute.org/software/igv/
 */
public final class IGVUtils {
    private IGVUtils() {}

    /**
     * Prints a header for an IGV-format file. See https://software.broadinstitute.org/software/igv/IGV
     * for details of this format.
     *
     * @param out stream to print the header to
     * @param graphType type of graph (eg., "line")
     * @param columns Column labels for the 5th and subsequent columns, representing the individual "tracks"
     *                in the IGV file. Note that the first 4 column labels in an IGV file are
     *                required to be Chromosome/Start/End/Feature, and should NOT be passed in here.
     */
    public static void printIGVFormatHeader(final PrintStream out, final String graphType, final String ... columns) {
        out.printf("#track graphType=%s%n", graphType);
        out.printf("Chromosome\tStart\tEnd\tFeature\t%s%n", Utils.join("\t", (Object[])columns));
    }

    /**
     * Prints a row in an IGV-format file. See https://software.broadinstitute.org/software/igv/IGV
     * for details of this format.
     *
     * @param out stream to print the row to
     * @param loc genome location associated with this row (in 1-based closed form)
     * @param featureName name of the feature associated with this row
     * @param values values for each track, corresponding to the 5th and subsequent column names declared in the header
     */
    public static void printIGVFormatRow(final PrintStream out, final Locatable loc, final String featureName, final double ... values) {
        // note that start and stop are 0-based in IGV files, but the stop is exclusive so we don't subtract 1 from it
        out.printf("%s\t%d\t%d\t%s", loc.getContig(), loc.getStart() - 1, loc.getEnd(), featureName);
        for ( final double value : values ) {
            out.print(String.format("\t%.5f", value));
        }
        out.println();
    }
}
