package org.broadinstitute.hellbender.engine;

import htsjdk.variant.bcf2.BCF2Codec;
import htsjdk.variant.bcf2.BCFVersion;
import org.broadinstitute.hellbender.exceptions.GATKException;

/**
 * GenomicsD uses htslib to generate BCF v2.2 streams that are specially modified for consumption
 * by GATK. GATK uses the htsjdk BCF2Codec to read these streams, but the BCF2Codec only
 * supports BCF 2.1. This class overrides the default, strict version compatibility check used
 * by the BCF2Codec to allow these hybrid streams that are marked a BCF v2.2 to used by the BCF
 * v2.1 codec directly.
 */
public class GenomicsDBBCFCodec extends BCF2Codec {

    @Override
    protected void validateVersionCompatibility(final BCFVersion supportedVersion, final BCFVersion actualVersion) {
        final BCFVersion expectedHTSJDKSupportedVersion = new BCFVersion(2, 1);
        final BCFVersion expectedGenomicsDBSupportedVersion = new BCFVersion(2, 2);

        // make sure the supported htsjdk version matches our expectations
        if (!supportedVersion.equals(expectedHTSJDKSupportedVersion)) {
            // hitting this may mean that this subclass is no longer needed
            throw new GATKException("The htsjdk supported BCF version does not match the expected supported version");
        }

        if (!actualVersion.equals(expectedGenomicsDBSupportedVersion)) {
            throw new IllegalStateException(
                    String.format("A stream supplied by GenomicsDB contains BCF version %s but GATK expects version %s",
                            new BCFVersion(actualVersion.getMajorVersion(), actualVersion.getMinorVersion()).toString(),
                            expectedGenomicsDBSupportedVersion.toString()));
        }
    }
}
