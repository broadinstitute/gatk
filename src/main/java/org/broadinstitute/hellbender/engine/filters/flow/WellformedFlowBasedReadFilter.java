package org.broadinstitute.hellbender.engine.filters.flow;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.engine.filters.AlignmentAgreesWithHeaderReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.filters.WellformedReadFilter;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * Tests whether a flow based read is &quot;well-formed&quot; -- that is, is free of major internal inconsistencies and issues that could lead
 * to errors downstream. If a read passes this filter, the rest of the engine should be able to process it without
 * blowing up. Note that checks already present in WellformedReadFilter are not duplicated here.
 *
 * <p><b>Well-formed flow based reads definition</b></p>
 * <ul>
 *     <li><b>Flow order: read group must have flow order</b></li>
 *     <li><b>Quality:</b> should be symtrical within each hmer.</li>
 *     <li><b>tp attribute:</b> should be symtrical within each hmer.</li>
 *     <li><b>tp attribute:</b> tp+hmer_length should be within [0, maxhmer],</li>
 *     <li><b>Hardclipped hmer:</B> is excempted from above checks.</li>
 *
 *      @see ReadGroupHasFlowOrderReadFilter
 *      @see HmerQualitySymetricReadFilter
 *      @see FlowBasedTPAttributeSymetricReadFilter
 *      @see FlowBasedTPAttributeValidReadFilter
 *
 * </ul>
 */
@DocumentedFeature(groupName=HelpConstants.DOC_CAT_READFILTERS,
        groupSummary=HelpConstants.DOC_CAT_READFILTERS_SUMMARY,
        summary = "Keep only flow based reads that are well-formed",
        extraDocs = {
                ReadGroupHasFlowOrderReadFilter.class,
                HmerQualitySymetricReadFilter.class,
                FlowBasedTPAttributeSymetricReadFilter.class,
                FlowBasedTPAttributeValidReadFilter.class
        }
)
public final class WellformedFlowBasedReadFilter extends ReadFilter {
    private static final long serialVersionUID = 1l;

    private ReadFilter wellFormedFilter = null;

    public WellformedFlowBasedReadFilter() {

    }

    @Override
    public void setHeader(SAMFileHeader header) {
        super.setHeader(header);
        createFilter();
    }

    private void createFilter() {

        wellFormedFilter = (new WellformedReadFilter(samHeader))
                .and(new ReadGroupHasFlowOrderReadFilter(samHeader))
                .and(ReadFilterLibrary.HMER_QUALITY_SYMETRIC_READ_FILTER)
                .and(ReadFilterLibrary.FLOW_BASED_TP_ATTRIBUTE_VALID_READ_FILTER)
                .and(ReadFilterLibrary.FLOW_BASED_TP_ATTRIBUTE_SYMETRIC_READ_FILTER);

    }

    @Override
    public boolean test(final GATKRead read ) {
        return wellFormedFilter.test(read);
    }
}
