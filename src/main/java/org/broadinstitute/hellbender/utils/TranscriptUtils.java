package org.broadinstitute.hellbender.utils;

import javassist.Loader;
import org.broadinstitute.hellbender.engine.ReferenceContext;

import java.util.Arrays;
import java.util.List;

/**
 * Class with utility methods for working with transcripts.
 */
public class TranscriptUtils {

    //==================================================================================================================
    // Public Static Members:

    /**
     * Extracts the transcript sequence from the reference context given the exons.
     * @param refContext the reference context.  Must not be {@code null}.
     * @param exons the exons of the transcript.  Must not be {@code null}.
     * @return the transcript sequence.
     */
    public static final String extractTrascriptFromReference(final ReferenceContext refContext, final List<SimpleInterval> exons) {

        Utils.nonNull(refContext);
        Utils.nonNull(exons);

        final StringBuilder transcript = new StringBuilder();
        for (final SimpleInterval exon : exons) {
            transcript.append(new String(refContext.getBases(exon)));
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
