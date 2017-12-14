package org.broadinstitute.hellbender.transformers;

import htsjdk.samtools.fastq.FastqConstants;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * Removes /1 or /2 and any whitespace from the end of the read name if present
 */

public class StripMateNumberTransformer implements ReadTransformer {
    public static final long serialVersionUID = 1L;

    @Override
    public GATKRead apply(final GATKRead read) {
        if (read == null) throw new IllegalArgumentException("Read cannot be null");
        final String name = read.getName().trim();
        if (name.endsWith(FastqConstants.FIRST_OF_PAIR) || name.endsWith(FastqConstants.SECOND_OF_PAIR)) {
            read.setName(name.substring(0, name.length() - 2).trim());
        }
        return read;
    }
}
