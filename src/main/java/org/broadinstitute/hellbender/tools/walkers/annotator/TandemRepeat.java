package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.util.*;

/**
 * Tandem repeat unit composition and counts per allele
 *
 * <p>This annotation tags variants that fall within tandem repeat sets. It also provides the composition of the tandem repeat units and the number of times they are repeated for each allele (including the REF allele).</p>
 *
 * <p>A tandem repeat unit is composed of one or more nucleotides that are repeated multiple times in series. Repetitive sequences are difficult to map to the reference because they are associated with multiple alignment possibilities. Knowing the number of repeat units in a set of tandem repeats tells you the number of different positions the tandem repeat can be placed in. The observation of many tandem repeat units multiplies the number of possible representations that can be made of the region.
 *
 * <h3>Caveat</h3>
 * <ul>
 *     <li>This annotation is currently not compatible with HaplotypeCaller.</li>
 * </ul>
 *
 */
@DocumentedFeature(groupName=HelpConstants.DOC_CAT_ANNOTATORS, groupSummary=HelpConstants.DOC_CAT_ANNOTATORS_SUMMARY, summary="Tandem repeat unit composition and counts per allele (STR, RU, RPA)")
public final class TandemRepeat extends InfoFieldAnnotation implements StandardMutectAnnotation {

    @Override
    public Map<String, Object> annotate(final ReferenceContext ref,
                                        final VariantContext vc,
                                        final ReadLikelihoods<Allele> likelihoods) {
        Utils.nonNull(vc);
        if ( !vc.isIndel()) {
            return Collections.emptyMap();
        }

        final Pair<List<Integer>,byte[]> result = GATKVariantContextUtils.getNumTandemRepeatUnits(vc, getRefBasesStartingAtVariantLocus(ref, vc));
        if (result == null) {
            return Collections.emptyMap();
        }

        final byte[] repeatUnit = result.getRight();
        final List<Integer> numUnits = result.getLeft();

        final Map<String, Object> map = new LinkedHashMap<>();
        map.put(GATKVCFConstants.STR_PRESENT_KEY, true);
        map.put(GATKVCFConstants.REPEAT_UNIT_KEY, new String(repeatUnit));
        map.put(GATKVCFConstants.REPEATS_PER_ALLELE_KEY, numUnits);
        return Collections.unmodifiableMap(map);
    }

    @Override
    public List<String> getKeyNames() {
        return Arrays.asList(
                GATKVCFConstants.STR_PRESENT_KEY,
                GATKVCFConstants.REPEAT_UNIT_KEY,
                GATKVCFConstants.REPEATS_PER_ALLELE_KEY);
    }

    private static byte[] getRefBasesStartingAtVariantLocus(final ReferenceContext ref, final VariantContext vc) {
        final byte[] bases = ref.getBases();
        final int startIndex = vc.getStart() - ref.getWindow().getStart();
        return new String(bases).substring(startIndex).getBytes();
    }

}
