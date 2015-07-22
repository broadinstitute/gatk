package org.broadinstitute.hellbender.tools.dataflow.pipelines;


import com.google.common.collect.ImmutableList;
import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.DataFlowProgramGroup;
import org.broadinstitute.hellbender.engine.dataflow.PTransformSAM;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.WellformedReadFilter;
import org.broadinstitute.hellbender.tools.FlagStat;
import org.broadinstitute.hellbender.tools.dataflow.transforms.FlagStatusDataflowTransform;


@CommandLineProgramProperties(summary ="runs FlagStat on dataflow", oneLineSummary = "FlagStat", programGroup = DataFlowProgramGroup.class)
public final class FlagStatDataflow extends DataflowReadsPipeline {
    private static final long serialVersionUID = 1l;

    @Override
    protected ImmutableList<ReadFilter> getReadFilters( SAMFileHeader header ){
        return ImmutableList.of(new WellformedReadFilter(header));
    }

    @Override
    protected PTransformSAM<FlagStat.FlagStatus> getTool() {
        return new FlagStatusDataflowTransform();
    }


}