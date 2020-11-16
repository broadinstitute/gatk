package org.broadinstitute.hellbender.tools.sv;

import com.google.common.collect.Sets;
import org.broadinstitute.hellbender.utils.codecs.BafEvidenceCodec;
import org.broadinstitute.hellbender.utils.codecs.DepthEvidenceCodec;
import org.broadinstitute.hellbender.utils.codecs.DiscordantPairEvidenceCodec;
import org.broadinstitute.hellbender.utils.codecs.SplitReadEvidenceCodec;

import java.util.Set;

public final class SVIOUtils {

    public static Set<String> EVIDENCE_FILE_EXTENSIONS = Sets.newHashSet(
            DiscordantPairEvidenceCodec.FORMAT_SUFFIX,
            SplitReadEvidenceCodec.FORMAT_SUFFIX,
            BafEvidenceCodec.FORMAT_SUFFIX,
            DepthEvidenceCodec.FORMAT_SUFFIX);

    public static String encodeSVEvidenceFeature(final Object feature) {
        if (feature instanceof BafEvidence) {
            return BafEvidenceCodec.encode((BafEvidence) feature);
        } else if (feature instanceof DepthEvidence) {
            return DepthEvidenceCodec.encode((DepthEvidence) feature);
        } else if (feature instanceof DiscordantPairEvidence) {
            return DiscordantPairEvidenceCodec.encode((DiscordantPairEvidence) feature);
        } else if (feature instanceof SplitReadEvidence) {
            return SplitReadEvidenceCodec.encode((SplitReadEvidence) feature);
        }
        throw new IllegalArgumentException("Unsupported SV evidence class: " + feature.getClass().getSimpleName());
    }
}
