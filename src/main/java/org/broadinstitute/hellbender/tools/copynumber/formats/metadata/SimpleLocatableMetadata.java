package org.broadinstitute.hellbender.tools.copynumber.formats.metadata;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.utils.Utils;

/**
 * Metadata associated with a collection of locatables.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public class SimpleLocatableMetadata implements LocatableMetadata {
    private final SAMSequenceDictionary sequenceDictionary;

    public SimpleLocatableMetadata(final SAMSequenceDictionary sequenceDictionary) {
        this.sequenceDictionary = Utils.nonNull(sequenceDictionary);
    }

    @Override
    public SAMSequenceDictionary getSequenceDictionary() {
        return sequenceDictionary;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (o == null || getClass() != o.getClass()) {
            return false;
        }

        final SimpleLocatableMetadata that = (SimpleLocatableMetadata) o;
        return sequenceDictionary.isSameDictionary(that.sequenceDictionary);
    }

    @Override
    public int hashCode() {
        return sequenceDictionary.hashCode();
    }

    @Override
    public String toString() {
        return "SimpleLocatableMetadata{" +
                "sequenceDictionary=" + sequenceDictionary +
                '}';
    }

    @Override
    public SAMFileHeader toHeader() {
        return new SAMFileHeader(sequenceDictionary);
    }
}
