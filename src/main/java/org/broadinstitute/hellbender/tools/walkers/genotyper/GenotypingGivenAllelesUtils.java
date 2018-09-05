package org.broadinstitute.hellbender.tools.walkers.genotyper;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.util.List;
import java.util.stream.Collectors;

import static org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils.FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED;
import static org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils.FilteredRecordMergeType.KEEP_UNCONDITIONAL;

/**
 * Compendium of utils to work in GENOTYPE_GIVEN_ALLELES mode.
 */
public final class GenotypingGivenAllelesUtils {

    /**
     * Composes the given allele variant-context providing information about the given variant alleles and reference location.
     * @param tracker the meta data tracker.
     * @param loc the query location.
     * @param keepFiltered whether to include filtered variants
     * @param allelesBinding the target variation context binding containing the given alleles.
     * @return never {@code null}
     */
    public static VariantContext composeGivenAllelesVariantContextFromVariantList(final FeatureContext tracker,
                                                                                  final Locatable loc,
                                                                                  final boolean keepFiltered,
                                                                                  final FeatureInput<VariantContext> allelesBinding) {
        Utils.nonNull(tracker, "tracker may not be null");
        Utils.nonNull(loc, "location may not be null");
        Utils.nonNull(allelesBinding, "alleles binding may not be null");

        final List<VariantContext> variantContextsInFeatureContext = tracker.getValues(allelesBinding, new SimpleInterval(loc));
        return composeGivenAllelesVariantContextFromVariantList(variantContextsInFeatureContext, loc, keepFiltered);
    }

    @VisibleForTesting
    protected static VariantContext composeGivenAllelesVariantContextFromVariantList(final List<VariantContext> variantContextsInFeatureContext,
                                                                                     final Locatable loc,
                                                                                     final boolean keepFiltered) {
        final List<VariantContext> vcsAtLoc = variantContextsInFeatureContext
                .stream()
                .filter(vc -> vc.getStart() == loc.getStart() &&
                        (keepFiltered || vc.isNotFiltered()))
                .collect(Collectors.toList());


        if (vcsAtLoc.isEmpty()) {
            return null;
        }
        final List<String> haplotypeSources = vcsAtLoc.stream().map(VariantContext::getSource).collect(Collectors.toList());
        final VariantContext mergedVc = GATKVariantContextUtils.simpleMerge(vcsAtLoc, haplotypeSources,
                keepFiltered ? KEEP_UNCONDITIONAL : KEEP_IF_ANY_UNFILTERED,
                GATKVariantContextUtils.GenotypeMergeType.PRIORITIZE, false, false, null, false, false);

        return mergedVc;
    }

}
