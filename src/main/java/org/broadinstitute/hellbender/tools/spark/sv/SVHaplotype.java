package org.broadinstitute.hellbender.tools.spark.sv;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignmentInterval;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.util.HashSet;
import java.util.List;
import java.util.Objects;

/**
 * Created by valentin on 9/16/17.
 */
public class SVHaplotype extends Haplotype {

    private static final long serialVersionUID = 1L;

    public String getVariantId() {
        return variantId;
    }

    public boolean isContig() {
        return isContig;
    }

    public double getReferenceScore() {
        if (name.equals("ref")) {
            return 0.0;
        } else if (name.equals("alt")) {
            return Double.NEGATIVE_INFINITY;
        } else {
            return referenceAlignmentScore.getValue();
        }
    }

    public double getAlternativeScore() {
        if (name.equals("alt")) {
            return 0.0;
        } else if (name.equals("ref")) {
            return Double.NEGATIVE_INFINITY;
        } else {
            return alternativeAlignmentScore.getValue();
        }
    }

    enum Call {
        REF, ALT, NOCALL;
    }

    private final String name;
    private final String variantId;
    private final List<AlignmentInterval> referenceAlignment;
    private final List<AlignmentInterval> alternativeAlignment;
    private final AlignmentScore referenceAlignmentScore;
    private final AlignmentScore alternativeAlignmentScore;
    private final double callQuality;
    private final Call call;
    private final boolean isContig;

    public static SVHaplotype of(final GATKRead read) {
        final String variantId = getMandatoryAttribute(read, ComposeStructuralVariantHaplotypesSpark.VARIANT_CONTEXT_TAG);
        final List<AlignmentInterval> refAln = getAlignmentIntervalsAttribute(read, ComposeStructuralVariantHaplotypesSpark.REFERENCE_ALIGNMENT_TAG);
        final List<AlignmentInterval> altAln = getAlignmentIntervalsAttribute(read, ComposeStructuralVariantHaplotypesSpark.ALTERNATIVE_ALIGNMENT_TAG);
        final AlignmentScore refScore = getMandatoryAlignmentScore(read, ComposeStructuralVariantHaplotypesSpark.REFERENCE_SCORE_TAG);
        final AlignmentScore altScore = getMandatoryAlignmentScore(read, ComposeStructuralVariantHaplotypesSpark.ALTERNATIVE_SCORE_TAG);
        final boolean isContig = read.getReadGroup().equals("CTG");
        final SimpleInterval location = new SimpleInterval(read.getAssignedContig(), read.getAssignedStart(), read.getAssignedStart());
        return new SVHaplotype(read.getName(), location, variantId, isContig, read.getBases(), refAln, refScore, altAln, altScore);

    }

    private static AlignmentScore getMandatoryAlignmentScore(final GATKRead read, final String tag) {
        return AlignmentScore.valueOf(getMandatoryAttribute(read, tag));
    }


    private static String getMandatoryAttribute(final GATKRead read, final String tag) {
        return ReadUtils.getOptionalStringAttribute(read, tag)
                .orElseThrow(() -> new UserException.BadInput("input read missing '" + tag + "' attribute"));
    }

    private static List<AlignmentInterval> getAlignmentIntervalsAttribute(final GATKRead read, final String tag) {
        return AlignmentInterval.decodeList(getMandatoryAttribute(read, tag));
    }

    private SVHaplotype(final String name, final Locatable loc, final String variantId, final boolean isContig,
                        final byte[] bases, final List<AlignmentInterval> refAln, final AlignmentScore refScore,
                        final List<AlignmentInterval> altAln, final AlignmentScore altScore) {
        super(bases, loc);
        this.name = name;
        this.isContig = isContig;
        this.variantId = variantId;
        this.referenceAlignment = refAln;
        this.alternativeAlignment = altAln;
        this.referenceAlignmentScore = refScore;
        this.alternativeAlignmentScore = altScore;
        final double refScoreValue = refScore.getValue();
        final double altScoreValue = altScore.getValue();
        if (refScoreValue < altScoreValue) {
            this.call = Call.ALT;
            this.callQuality = altScoreValue - refScoreValue;
        } else if (altScoreValue < refScoreValue) {
            this.call = Call.REF;
            this.callQuality = refScoreValue - altScoreValue;
        } else {
            this.call = Call.NOCALL;
            this.callQuality = 0.0;
        }
    }

    public String getName() {
        return name;
    }

}
