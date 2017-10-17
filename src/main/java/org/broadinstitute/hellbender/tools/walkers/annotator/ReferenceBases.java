package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.logging.OneShotLogger;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;

/**
 * Local reference context at a variant position.
 *
 * </p>The annotation gives ten reference bases each to the left and right of the variant start and the start base for a total of 21 reference bases.
 * Start position is defined as one base before indels.  For example, the reference context AAAAAAAAAACTTTTTTTTTT would apply to a SNV variant
 * context with ref allele C and alt allele G as well as to a deletion variant context with ref allele CT and alt allele C.</p>
 *
 */
@DocumentedFeature(groupName=HelpConstants.DOC_CAT_ANNOTATORS, groupSummary=HelpConstants.DOC_CAT_ANNOTATORS_SUMMARY, summary="Annotate with local reference bases (REF_BASES)")
public class ReferenceBases extends InfoFieldAnnotation {
    public static final String REFERENCE_BASES_KEY = "REF_BASES";

    public static final int NUM_BASES_ON_EITHER_SIDE = 3;

    protected final OneShotLogger warning = new OneShotLogger(this.getClass());

    @Override
    public List<String> getKeyNames() { return Collections.singletonList(REFERENCE_BASES_KEY); }

    @Override
    public Map<String, Object> annotate(final ReferenceContext ref,
                                        final VariantContext vc,
                                        final ReadLikelihoods<Allele> likelihoods) {
        if (ref==null)  {
            warning.warn("REF_BASES requires the reference to annotate, none was provided");
            return Collections.emptyMap();
        }
        final int basesToDiscardInFront = Math.max(vc.getStart() - ref.getWindow().getStart() - NUM_BASES_ON_EITHER_SIDE, 0);
        final String allBases = new String(ref.getBases());
        final String localBases = allBases.substring(basesToDiscardInFront, basesToDiscardInFront + 2 * NUM_BASES_ON_EITHER_SIDE);
        return Collections.singletonMap(REFERENCE_BASES_KEY, localBases );
    }

    @Override
    public List<VCFInfoHeaderLine> getDescriptions() {
        return Arrays.asList(new VCFInfoHeaderLine(ReferenceBases.REFERENCE_BASES_KEY, 1, VCFHeaderLineType.String, "local reference bases."));
    }
}
