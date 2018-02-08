package org.broadinstitute.hellbender.tools.copynumber.formats.metadata;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.copynumber.coverage.readcount.ReadCountFileHeaderKey;
import org.broadinstitute.hellbender.tools.copynumber.coverage.readcount.ReadCountType;
import org.broadinstitute.hellbender.tools.copynumber.coverage.readcount.covariatebin.ReadCountCovariateBinningConfiguration;
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

    /**
     * @return read count type extracted from {@link SAMFileHeader} comments
     */
    public static ReadCountType getReadCountType(final SAMFileHeader header) {
        Utils.nonNull(header);
        final List<String> comments = Utils.nonEmpty(header.getComments());
        return ReadCountType.getReadCountTypeByName(
                ReadCountFileHeaderKey.getHeaderValueForKey(comments, ReadCountFileHeaderKey.READ_COUNT_TYPE));
    }

    /**
     * @return null if read count type is {@link ReadCountType#SIMPLE_COUNT}
     */
    public static List<ReadCountCovariateBinningConfiguration> getBinningConfigurations(
            final SAMFileHeader header,
            final ReadCountType readCountType) {

        Utils.nonNull(header);
        if (readCountType == ReadCountType.SIMPLE_COUNT) {
            return null;
        }

        final List<String> comments = Utils.nonEmpty(header.getComments());
        return ReadCountCovariateBinningConfiguration.parseParameters(
                ReadCountFileHeaderKey.getHeaderValueForKey(comments, ReadCountFileHeaderKey.BINNING_CONFIGURATION));
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
            case BINNING_SAMPLE_LOCATABLE:
                final ReadCountType readCountType = getReadCountType(header);
                return (T) new SimpleBinningSampleLocatableMetadata(readSampleName(header),
                        header.getSequenceDictionary(), readCountType, getBinningConfigurations(header, readCountType));
            default:
                throw new GATKException.ShouldNeverReachHereException("Encountered unknown Metadata.Type.");
        }
    }
}
