package org.broadinstitute.hellbender.tools.dataflow.pipelines;


import com.google.common.collect.ImmutableList;
import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.DataFlowProgramGroup;
import org.broadinstitute.hellbender.engine.dataflow.PTransformSAM;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.WellformedReadFilter;
import org.broadinstitute.hellbender.tools.dataflow.transforms.CountBasesDataflowTransform;


@CommandLineProgramProperties(summary = "Count bases in dataflow with filters and transformers", oneLineSummary = "count bases", programGroup = DataFlowProgramGroup.class)
public final class CountBasesDataflow extends DataflowReadsPipeline {
    private static final long serialVersionUID = 1l;

    @Override
    protected ImmutableList<ReadFilter> getReadFilters( final SAMFileHeader header ){
        return ImmutableList.of(new WellformedReadFilter(header));
    }

    @Override
    protected PTransformSAM<Long> getTool() {
        return new CountBasesDataflowTransform();
    }

}



