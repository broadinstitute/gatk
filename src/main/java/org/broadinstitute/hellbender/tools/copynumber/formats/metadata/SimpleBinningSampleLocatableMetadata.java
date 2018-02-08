package org.broadinstitute.hellbender.tools.copynumber.formats.metadata;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.tools.copynumber.coverage.readcount.ReadCountFileHeaderKey;
import org.broadinstitute.hellbender.tools.copynumber.coverage.readcount.ReadCountType;
import org.broadinstitute.hellbender.tools.copynumber.coverage.readcount.covariatebin.ReadCountCovariateBinningConfiguration;
import org.broadinstitute.hellbender.utils.Utils;

import javax.annotation.Nullable;
import java.util.List;

/**
 * Concrete implementation of {@link BinningSampleLocatableMetadata} that stores metadata associated
 * with a collection of locatables and with a single sample, as well as read count type and list
 * of binning configurations
 *
 * @author Andrey Smirnov &lt;asmirnov@broadinstitute.org&gt;
 */
public class SimpleBinningSampleLocatableMetadata implements BinningSampleLocatableMetadata {

    private final SampleLocatableMetadata sampleLocatableMetadata;
    private final ReadCountType readCountType;
    private final List<ReadCountCovariateBinningConfiguration> binningConfigurations;

    public SimpleBinningSampleLocatableMetadata(
            final String sampleName,
            final SAMSequenceDictionary sequenceDictionary,
            final ReadCountType readCountType,
            @Nullable final List<ReadCountCovariateBinningConfiguration> binningConfigurations) {
        sampleLocatableMetadata = new SimpleSampleLocatableMetadata(sampleName, sequenceDictionary);
        this.readCountType = Utils.nonNull(readCountType);
        this.binningConfigurations = binningConfigurations;
    }

    @Override
    public List<ReadCountCovariateBinningConfiguration> getCovariateBinningConfigurations() {
        return binningConfigurations;
    }

    @Override
    public ReadCountType getReadCountType() {
        return readCountType;
    }

    @Override
    public String getSampleName() {
        return sampleLocatableMetadata.getSampleName();
    }

    @Override
    public SAMSequenceDictionary getSequenceDictionary() {
        return sampleLocatableMetadata.getSequenceDictionary();
    }

    @Override
    public SAMFileHeader toHeader() {
        final SAMFileHeader header = sampleLocatableMetadata.toHeader();
        header.addComment(ReadCountFileHeaderKey.constructReadCountTypeCommentValue(readCountType));
        if (binningConfigurations != null) {
            header.addComment(
                    ReadCountFileHeaderKey.constructCovariateBinningConfigurationComment(binningConfigurations));
        }
        return header;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (o == null || getClass() != o.getClass()) {
            return false;
        }

        final SimpleBinningSampleLocatableMetadata that = (SimpleBinningSampleLocatableMetadata) o;

        if (!sampleLocatableMetadata.equals(that.sampleLocatableMetadata)) {
            return false;
        }
        if (readCountType != that.readCountType) {
            return false;
        }
        return binningConfigurations.equals(that.binningConfigurations);

    }

    @Override
    public int hashCode() {
        int result = sampleLocatableMetadata.hashCode();
        result = 31 * result + readCountType.hashCode();
        result = 31 * result + binningConfigurations.hashCode();
        return result;
    }

    @Override
    public String toString() {
        return "SimpleBinningSampleLocatableMetadata{" +
                "sampleLocatableMetadata=" + sampleLocatableMetadata +
                ", readCountType=" + readCountType +
                ", binningConfigurations=" + binningConfigurations +
                '}';
    }
}
