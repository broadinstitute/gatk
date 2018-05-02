package org.broadinstitute.hellbender.tools.walkers.genotyper;

import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.util.List;
import java.util.stream.Collectors;

/**
 * Compendium of utils to work in GENOTYPE_GIVEN_ALLELES mode.
 */
public final class GenotypingGivenAllelesUtils {

    /**
     * Composes the given allele variant-context providing information about the rods and reference location.
     * @param tracker the meta data tracker.
     * @param loc the query location.
     * @param snpsOnly whether we only should consider SNP variation.
     * @param keepFiltered whether to include filtered variants
     * @param logger where to output warnings.
     * @param allelesBinding the target variation context binding containing the given alleles.
     * @return never {@code null}
     */
    public static VariantContext composeGivenAllelesVariantContextFromRod(final FeatureContext tracker,
                                                                          final Locatable loc,
                                                                          final boolean snpsOnly,
                                                                          final boolean keepFiltered,
                                                                          final Logger logger,
                                                                          final FeatureInput<VariantContext> allelesBinding) {
        Utils.nonNull(tracker, "tracker may not be null");
        Utils.nonNull(loc, "location may not be null");
        Utils.nonNull(allelesBinding, "alleles binding may not be null");

        final List<VariantContext> rodVcsAtLoc = tracker.getValues(allelesBinding, new SimpleInterval(loc))
                .stream()
                .filter(vc -> vc.getStart() == loc.getStart() &&
                        (keepFiltered || vc.isNotFiltered()) &&
                        (!snpsOnly || vc.isSNP()))
                .collect(Collectors.toList());


        if (rodVcsAtLoc.isEmpty()) {
            return null;
        }
        final List<String> haplotypeSources = rodVcsAtLoc.stream().map(VariantContext::getSource).collect(Collectors.toList());
        final VariantContext mergedVc = GATKVariantContextUtils.simpleMerge(rodVcsAtLoc, haplotypeSources,
                GATKVariantContextUtils.FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED,
                GATKVariantContextUtils.GenotypeMergeType.PRIORITIZE, false, false, null, false, false);

        return mergedVc;
    }

}
