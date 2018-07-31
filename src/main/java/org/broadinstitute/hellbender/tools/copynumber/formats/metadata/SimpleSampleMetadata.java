package org.broadinstitute.hellbender.tools.copynumber.formats.metadata;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import org.broadinstitute.hellbender.utils.Utils;

/**
 * Metadata associated with a single sample.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public class SimpleSampleMetadata implements SampleMetadata {
    private final String sampleName;

    public SimpleSampleMetadata(final String sampleName) {
        this.sampleName = Utils.nonEmpty(sampleName);
    }

    @Override
    public String getSampleName() {
        return sampleName;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (o == null || getClass() != o.getClass()) {
            return false;
        }

        final SimpleSampleMetadata that = (SimpleSampleMetadata) o;
        return sampleName.equals(that.sampleName);
    }

    @Override
    public int hashCode() {
        return sampleName.hashCode();
    }

    @Override
    public String toString() {
        return "SimpleSampleMetadata{" +
                "sampleName='" + sampleName + '\'' +
                '}';
    }

    @Override
    public SAMFileHeader toHeader() {
        final SAMReadGroupRecord readGroupRecord = new SAMReadGroupRecord(MetadataUtils.GATK_CNV_READ_GROUP_ID);
        readGroupRecord.setAttribute(SAMReadGroupRecord.READ_GROUP_SAMPLE_TAG, sampleName);
        final SAMFileHeader header = new SAMFileHeader();
        header.addReadGroup(readGroupRecord);
        return header;
    }
}
