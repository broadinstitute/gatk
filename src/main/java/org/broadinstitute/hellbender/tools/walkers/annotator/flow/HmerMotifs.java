package org.broadinstitute.hellbender.tools.walkers.annotator.flow;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.Arrays;
import java.util.List;
import java.util.Map;

/**
 * Flow Annotation: motifs to the left and right of the indel
 */
@DocumentedFeature(groupName=HelpConstants.DOC_CAT_FLOW_ANNOTATORS, groupSummary=HelpConstants.DOC_CAT_FLOW_ANNOTATORS_SUMMARY, summary="Right Motif Flow Annotation")
public class HmerMotifs extends FlowAnnotatorBase implements StandardFlowBasedAnnotation {

    @Override
    public Map<String, Object> annotate(ReferenceContext ref,
                                        VariantContext vc,
                                        AlleleLikelihoods<GATKRead, Allele> likelihoods) {

        final LocalContext        localContext = new LocalContext(ref, vc, likelihoods, true);

        if ( localContext.generateAnnotation ) {
            getLeftMotif(vc, localContext);
            indelClassify(vc, localContext);
            isHmerIndel(vc, localContext);
            getRightMotif(vc, localContext);
        }

        return localContext.asAttributes();
    }

    @Override
    public List<String> getKeyNames() {

        return Arrays.asList(GATKVCFConstants.FLOW_LEFT_MOTIF, GATKVCFConstants.FLOW_RIGHT_MOTIF);
    }


}

