package org.broadinstitute.hellbender.tools.spark.sv;

import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceRecord;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariationSparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import scala.Tuple2;

import java.util.List;
import java.util.stream.Collectors;

/**
 * Tool to run the sv pipeline up to and including variant discovery
 */
@CommandLineProgramProperties(summary="Master tool to run the structural variation discovery pipeline",
        oneLineSummary="Master tool to run the structural variation discovery pipeline",
        programGroup = StructuralVariationSparkProgramGroup.class)
public class StructuralVariationDiscoveryPipelineSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;
    private final Logger localLogger = LogManager.getLogger(StructuralVariationDiscoveryPipelineSpark.class);


    @Argument(doc = "sam file for aligned contigs", shortName = "contigSAMFile",
            fullName = "contigSAMFile")
    private String outputSAM;

    @Argument(doc = "filename for output vcf", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    private String vcfOutputFileName;

    @ArgumentCollection
    private StructuralVariationDiscoveryArgumentCollection.FindBreakpointEvidenceSparkArgumentCollection reAlignmentStageArgs
            = new StructuralVariationDiscoveryArgumentCollection.FindBreakpointEvidenceSparkArgumentCollection();

    @ArgumentCollection
    private StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection discoverStageArgs
            = new StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection();


    @Override
    public boolean requiresReads()
    {
        return true;
    }

    @Override
    public boolean requiresReference() {return true;}

    @Override
    protected void runTool( final JavaSparkContext ctx ) {

        final SAMFileHeader header = getHeaderForReads();
        final PipelineOptions pipelineOptions = getAuthenticatedGCSOptions();
        final List<String> refNames = new SAMFileHeader(header.getSequenceDictionary())
                .getSequenceDictionary().getSequences().stream()
                .map(SAMSequenceRecord::getSequenceName).collect(Collectors.toList());

        final List<AlignedAssemblyOrExcuse> alignedAssemblyOrExcuseList =
                FindBreakpointEvidenceSpark
                        .gatherEvidenceAndWriteContigSamFile(ctx, pipelineOptions, header, getUnfilteredReads(),
                                reAlignmentStageArgs, localLogger, outputSAM);
        if (alignedAssemblyOrExcuseList.isEmpty()) return;

        final List<Tuple2<Iterable<AlignmentRegion>, byte[]>> parsedAlignments
                = AssemblyAlignmentParser.formatToAlignmentRegions(alignedAssemblyOrExcuseList, refNames);
        if (parsedAlignments.isEmpty()) return;

        DiscoverVariantsFromContigAlignmentsSpark
                .makeSenseAndWrite(ctx.parallelizePairs(parsedAlignments),
                        discoverStageArgs.fastaReference, ctx.broadcast(getReference()),
                        pipelineOptions, vcfOutputFileName, localLogger);

    }

}
