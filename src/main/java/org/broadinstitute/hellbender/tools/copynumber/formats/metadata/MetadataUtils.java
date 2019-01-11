package org.broadinstitute.hellbender.tools.copynumber.formats.metadata;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.List;
import java.util.stream.Collectors;

/**
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class MetadataUtils {

    static final String GATK_CNV_READ_GROUP_ID = "GATKCopyNumber";

    private MetadataUtils() {}

    public static String readSampleName(final SAMFileHeader header) {
        Utils.nonNull(header);
        Utils.nonEmpty(header.getReadGroups(), "The input header does not contain any read groups.  Cannot determine a sample name.");
        final List<String> sampleNames = header.getReadGroups().stream().map(SAMReadGroupRecord::getSample).distinct().collect(Collectors.toList());
        if (sampleNames.size() > 1) {
            throw new IllegalArgumentException(String.format("The input header contains more than one unique sample name: %s",
                    StringUtils.join(sampleNames, ", ")));
        }
        if (sampleNames.isEmpty()) {
            throw new IllegalArgumentException("The input header does not contain a sample name.");
        }
        return sampleNames.get(0);
    }

    @SuppressWarnings("unchecked")
    public static <T extends Metadata> T fromHeader(final SAMFileHeader header,
                                                    final Metadata.Type metadataType) {
        Utils.nonNull(header);
        Utils.nonNull(metadataType);
        switch (metadataType) {
            case SAMPLE:
                return (T) new SimpleSampleMetadata(readSampleName(header));
            case LOCATABLE:
                return (T) new SimpleLocatableMetadata(header.getSequenceDictionary());
            case SAMPLE_LOCATABLE:
                return (T) new SimpleSampleLocatableMetadata(readSampleName(header), header.getSequenceDictionary());
            default:
                throw new GATKException.ShouldNeverReachHereException("Encountered unknown Metadata.Type.");
        }
    }
}
