package org.broadinstitute.hellbender.tools.walkers.annotator;

import com.google.common.primitives.Ints;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.List;
import java.util.OptionalInt;

/**
 * Median mapping quality of reads supporting each alt allele.
 *
 * <p>The output is an array containing, for each alt allele, the median mapping quality over all reads that best match that allele.</p>
 */
@DocumentedFeature(groupName=HelpConstants.DOC_CAT_ANNOTATORS, groupSummary=HelpConstants.DOC_CAT_ANNOTATORS_SUMMARY, summary="Median mapping quality of reads supporting each allele (MMQ)")
public class MappingQuality extends PerAlleleAnnotation implements StandardMutectAnnotation {

    // we don't want a GGA mode allele with no reads to prejudice us against a site so we assign a good mapping quality
    private static final int VALUE_FOR_NO_READS = 60;

    @Override
    protected int aggregate(final List<Integer> values) {
        return values.isEmpty() ? VALUE_FOR_NO_READS : MathUtils.median(Ints.toArray(values));
    }

    @Override
    protected boolean includeRefAllele() { return true; }

    @Override
    protected String getVcfKey() { return GATKVCFConstants.MEDIAN_MAPPING_QUALITY_KEY; }

    @Override
    protected String getDescription() { return "median mapping quality"; }

    @Override
    protected OptionalInt getValueForRead(final GATKRead read, final VariantContext vc) {
        Utils.nonNull(read);
        return OptionalInt.of(read.getMappingQuality());
    }
}
