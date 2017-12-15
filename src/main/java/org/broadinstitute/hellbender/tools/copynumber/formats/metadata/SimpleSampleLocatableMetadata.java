package org.broadinstitute.hellbender.tools.copynumber.formats.metadata;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.utils.Utils;

/**
 * Metadata associated with a collection of locatables associated with a single sample.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public class SimpleSampleLocatableMetadata implements SampleLocatableMetadata {
    private final SampleMetadata sampleMetadata;
    private final LocatableMetadata locatableMetadata;

    public SimpleSampleLocatableMetadata(final String sampleName,
                                         final SAMSequenceDictionary sequenceDictionary) {
        Utils.nonEmpty(sampleName);
        Utils.nonNull(sequenceDictionary);
        this.sampleMetadata = new SimpleSampleMetadata(sampleName);
        this.locatableMetadata = new SimpleLocatableMetadata(sequenceDictionary);
    }

    @Override
    public String getSampleName() {
        return sampleMetadata.getSampleName();
    }

    @Override
    public SAMSequenceDictionary getSequenceDictionary() {
        return locatableMetadata.getSequenceDictionary();
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (o == null || getClass() != o.getClass()) {
            return false;
        }

        final SimpleSampleLocatableMetadata that = (SimpleSampleLocatableMetadata) o;
        return sampleMetadata.equals(that.sampleMetadata) &&
                locatableMetadata.equals(that.locatableMetadata);
    }

    @Override
    public int hashCode() {
        int result = sampleMetadata.hashCode();
        result = 31 * result + locatableMetadata.hashCode();
        return result;
    }

    @Override
    public String toString() {
        return "SimpleSampleLocatableMetadata{" +
                "sampleMetadata=" + sampleMetadata +
                ", locatableMetadata=" + locatableMetadata +
                '}';
    }

    @Override
    public SAMFileHeader toHeader() {
        final SAMFileHeader header = sampleMetadata.toHeader();
        header.setSequenceDictionary(locatableMetadata.getSequenceDictionary());
        return header;
    }
}
