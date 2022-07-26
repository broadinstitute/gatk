package org.broadinstitute.hellbender.tools.walkers.annotator.flow;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.Collections;
import java.util.List;

/**
 * Flow Annotation: length of indel
 */
@DocumentedFeature(groupName=HelpConstants.DOC_CAT_FLOW_ANNOTATORS, groupSummary=HelpConstants.DOC_CAT_FLOW_ANNOTATORS_SUMMARY, summary="Indel Length Flow Annotation")
public class IndelLength extends IndelClassify implements StandardFlowBasedAnnotation {
    private final Logger logger = LogManager.getLogger(IndelLength.class);

    @Override
    public List<String> getKeyNames() {

        return Collections.singletonList(GATKVCFConstants.FLOW_INDEL_LENGTH);
    }


}

