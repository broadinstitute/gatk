package org.broadinstitute.hellbender.engine.spark.datasources;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.engine.BasicReference;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;

import java.io.IOException;

/**
 * Internal interface to load a reference sequence.
 */
public interface ReferenceSparkSource extends BasicReference {
    ReferenceBases getReferenceBases(SimpleInterval interval) throws IOException;
    SAMSequenceDictionary getReferenceSequenceDictionary(SAMSequenceDictionary optReadSequenceDictionaryToMatch) throws IOException;

    default byte[] getBases( final SimpleInterval window ) {
        try {
            return getReferenceBases(window).getBases();
        } catch ( final IOException ioe ) {
            throw new UserException("Unable to fetch reference bases for " + window, ioe);
        }
    }
}
