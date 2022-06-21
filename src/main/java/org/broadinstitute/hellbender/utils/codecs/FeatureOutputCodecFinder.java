package org.broadinstitute.hellbender.utils.codecs;

import htsjdk.tribble.Feature;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.util.ArrayList;
import java.util.List;

/**
 * This class knows about all FeatureOutputCodec implementations, and allows you to find an
 * appropriate codec to create a given file type.
 */
public final class FeatureOutputCodecFinder {

    private static final List<FeatureOutputCodec<? extends Feature, ? extends FeatureSink<?>>> outputCodecs =
            new ArrayList<>(10);
    static {
        outputCodecs.add(new BafEvidenceCodec());
        outputCodecs.add(new DepthEvidenceCodec());
        outputCodecs.add(new DiscordantPairEvidenceCodec());
        outputCodecs.add(new SiteDepthCodec());
        outputCodecs.add(new SplitReadEvidenceCodec());
        outputCodecs.add(new BafEvidenceBCICodec());
        outputCodecs.add(new DepthEvidenceBCICodec());
        outputCodecs.add(new DiscordantPairEvidenceBCICodec());
        outputCodecs.add(new SiteDepthBCICodec());
        outputCodecs.add(new SplitReadEvidenceBCICodec());
    }

    public static FeatureOutputCodec<? extends Feature, ? extends FeatureSink<? extends Feature>>
                    find( final GATKPath outputFilePath ) {
        FeatureOutputCodec<?, ?> result = null;
        final String outputFileName = outputFilePath.toString();
        for ( final FeatureOutputCodec<?, ?> codec : outputCodecs ) {
            if ( codec.canDecode(outputFileName) ) {
                if ( result != null ) {
                    throw new GATKException("Found multiple output codecs for " + outputFileName);
                }
                result = codec;
            }
        }
        if ( result == null ) {
            throw new UserException("No feature output codec found for " + outputFileName);
        }
        return result;
    }
}
