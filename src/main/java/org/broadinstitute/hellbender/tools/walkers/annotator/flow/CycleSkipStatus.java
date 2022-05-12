package org.broadinstitute.hellbender.tools.walkers.annotator.flow;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.annotator.StandardMutectAnnotation;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.logging.OneShotLogger;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.Collections;
import java.util.List;
import java.util.Map;

@DocumentedFeature(groupName=HelpConstants.DOC_CAT_FLOW_ANNOTATORS, groupSummary=HelpConstants.DOC_CAT_FLOW_ANNOTATORS_SUMMARY, summary="Cycle Skip Status Flow Annotation")
public class CycleSkipStatus extends FlowAnnotatorBase implements StandardFlowBasedAnnotation {
    private final Logger logger = LogManager.getLogger(CycleSkipStatus.class);

    @Override
    public Map<String, Object> annotate(ReferenceContext ref,
                                        VariantContext vc,
                                        AlleleLikelihoods<GATKRead, Allele> likelihoods) {

        final LocalContext        localContext = new LocalContext(ref, vc, likelihoods, true);

        if ( localContext.generateAnnotation ) {
            indelClassify(vc, localContext);
            isHmerIndel(vc, localContext);
            getLeftMotif(vc, localContext);
            getRightMotif(vc, localContext);
            cycleSkip(vc, localContext);
        }

        return localContext.asAttributes();
    }

    @Override
    public List<String> getKeyNames() {

        return Collections.singletonList(GATKVCFConstants.FLOW_CYCLESKIP_STATUS);
    }

    protected boolean isActualFlowOrderRequired() {
        return true;
    }


}

