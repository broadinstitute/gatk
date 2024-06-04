package org.broadinstitute.hellbender.utils;

import htsjdk.samtools.util.Locatable;
import javassist.Loader;
import org.broadinstitute.hellbender.engine.ReferenceContext;

import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Class with utility methods for working with transcripts.
 */
public class TranscriptUtils {

    //==================================================================================================================
    // Public Static Members:

    /**
     * Extracts the transcript sequence from the reference context given the exons.
     * The given exons are assumed to be all on the same contig.
     * The exons will be sorted in the order they appear in the genome before extracting the sequences.
     * @param refContext the reference context.  Must not be {@code null}.
     * @param exons the exons of the transcript.  Must not be {@code null}.
     * @return the transcript sequence as coded in the given reference context.
     */
    public static final String extractTrascriptFromReference(final ReferenceContext refContext, final List<Locatable> exons) {

        Utils.nonNull(refContext);
        Utils.nonNull(exons);

        final StringBuilder transcript = new StringBuilder();

        // We should iterate through the list of exons in sorted order so we can simply append them together.
        for (final Locatable exon : exons.stream().sorted(Comparator.comparingInt(Locatable::getStart).thenComparing(Locatable::getEnd)).toList() ) {
            transcript.append(new String(refContext.getBases(new SimpleInterval(exon))));
        }
        return transcript.toString();
    }

    //==================================================================================================================
    // Private Static Members:

    //==================================================================================================================
    // Private Members:

    //==================================================================================================================
    // Constructors:

    //==================================================================================================================
    // Override Methods:

    //==================================================================================================================
    // Static Methods:

    //==================================================================================================================
    // Instance Methods:

    //==================================================================================================================
    // Helper Data Types:

}
