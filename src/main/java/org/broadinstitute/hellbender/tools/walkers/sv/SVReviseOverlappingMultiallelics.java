package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;

import java.util.List;
import java.util.ArrayList;
import java.util.Set;
import java.util.HashSet;

@CommandLineProgramProperties(
        summary = "Clean and format SV VCF",
        oneLineSummary = "Clean and format SV VCF",
        programGroup = StructuralVariantDiscoveryProgramGroup.class
)
@BetaFeature
@DocumentedFeature
public class SVReviseOverlappingMultiallelics extends MultiplePassVariantWalker {
    @Argument(
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "Output VCF name"
    )
    private GATKPath outputVcf;

    private VariantContextWriter vcfWriter;

    private final List<VariantContext> overlappingVariantsBuffer = new ArrayList<>();
    private final List<OverlapCandidate> overlapCandidates = new ArrayList<>();
    private final Set<String> filteredVariantIds = new HashSet<>();

    @Override
    protected int numberOfPasses() { return 2; }

    @Override
    protected void afterNthPass(final int n) {
        if (n == 0) {
            resolveRedundantMultiallelics();
        }
    }

    /**
     * Resolves which of each overlapping pair of multiallelic sites is redundant. Candidates must be
     * processed in descending order of the driving (larger) variant's length, so that a larger variant's
     * own redundancy status is always settled before it is used to judge smaller overlapping variants.
     */
    private void resolveRedundantMultiallelics() {
        overlapCandidates.sort((a, b) -> {
            final int cmp = Integer.compare(b.largerLength, a.largerLength);
            if (cmp != 0) {
                return cmp;
            }
            return Integer.compare(b.smallerLength, a.smallerLength);
        });
        for (final OverlapCandidate candidate : overlapCandidates) {
            if (!filteredVariantIds.contains(candidate.largerId)) {
                filteredVariantIds.add(candidate.smallerId);
            }
        }
        overlapCandidates.clear();
        overlappingVariantsBuffer.clear();
    }

    @Override
    public void onTraversalStart() {
        vcfWriter = createVCFWriter(outputVcf);
        vcfWriter.writeHeader(getHeaderForVariants());
    }

    @Override
    public void closeTool() {
        if (vcfWriter != null) {
            vcfWriter.close();
        }
    }

    @Override
    protected void nthPassApply(final VariantContext variant, final ReadsContext readsContext,
                                final ReferenceContext referenceContext, final FeatureContext featureContext, int n) {
        switch (n) {
            case 0:
                firstPassApply(variant);
                break;
            case 1:
                secondPassApply(variant);
                break;
        }
    }

    public void firstPassApply(final VariantContext variant) {
        if (!variant.getFilters().contains(GATKSVVCFConstants.MULTIALLELIC)) {
            return;
        }

        overlappingVariantsBuffer.removeIf(vc -> !vc.getContig().equals(variant.getContig())
                || (vc.getStart() + vc.getAttributeAsInt(GATKSVVCFConstants.SVLEN, 0)) < variant.getStart());
        for (final VariantContext bufferedVariant : overlappingVariantsBuffer) {
            if (overlaps(bufferedVariant, variant)) {
                collectVariantPair(bufferedVariant, variant);
            }
        }
        overlappingVariantsBuffer.add(variant);
    }

    public void secondPassApply(final VariantContext variant) {
        if (filteredVariantIds.contains(variant.getID())) {
            return;
        }

        final VariantContextBuilder builder = new VariantContextBuilder(variant);
        vcfWriter.add(builder.make());
    }

    private void collectVariantPair(final VariantContext v1, final VariantContext v2) {
        // Determine larger variant, swapping if necessary
        VariantContext largerVariant = v1;
        VariantContext smallerVariant = v2;
        final int length1 = v1.getAttributeAsInt(GATKSVVCFConstants.SVLEN, 0);
        final int length2 = v2.getAttributeAsInt(GATKSVVCFConstants.SVLEN, 0);
        if (length2 > length1) {
            largerVariant = v2;
            smallerVariant = v1;
        }

        // Skip if coverage below expected
        final double coverage = getCoverage(largerVariant, smallerVariant);
        if (coverage < 0.5) {
            return;
        }

        overlapCandidates.add(new OverlapCandidate(largerVariant.getID(), smallerVariant.getID(),
                Math.max(length1, length2), Math.min(length1, length2)));
    }

    private static final class OverlapCandidate {
        private final String largerId;
        private final String smallerId;
        private final int largerLength;
        private final int smallerLength;

        private OverlapCandidate(final String largerId, final String smallerId,
                                  final int largerLength, final int smallerLength) {
            this.largerId = largerId;
            this.smallerId = smallerId;
            this.largerLength = largerLength;
            this.smallerLength = smallerLength;
        }
    }

    private boolean overlaps(final VariantContext v1, final VariantContext v2) {
        return v1.getContig().equals(v2.getContig())
                && v1.getStart() <= (v2.getStart() + v2.getAttributeAsInt(GATKSVVCFConstants.SVLEN, 0))
                && v2.getStart() <= (v1.getStart() + v1.getAttributeAsInt(GATKSVVCFConstants.SVLEN, 0));
    }

    private double getCoverage(final VariantContext larger, final VariantContext smaller) {
        final int largerStart = larger.getStart();
        final int smallerStart = smaller.getStart();
        final int largerStop = larger.getEnd();
        final int smallerStop = smaller.getEnd();

        if (largerStart <= smallerStop && smallerStart <= largerStop) {
            final int intersectionSize = Math.min(smallerStop, largerStop) - Math.max(smallerStart, largerStart) + 1;
            return (double) intersectionSize / (smallerStop - smallerStart + 1);
        }
        return 0.0;
    }
}
