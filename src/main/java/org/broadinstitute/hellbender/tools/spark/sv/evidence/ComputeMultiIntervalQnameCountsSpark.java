package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.utils.HopscotchUniqueMultiMap;

import java.util.List;

@BetaFeature
@CommandLineProgramProperties(
        oneLineSummary = "(Internal) Extracts evidence of structural variations from reads",
        summary =
                "This tool is used in development and should not be of interest to most researchers.  It packages one step" +
                        " of the structural variation workflow as a separate tool for the convenience of developers." +
                        " This tool examines a SAM/BAM/CRAM for reads, or groups of reads, that demonstrate evidence of a structural" +
                        " variation in the vicinity.  It records this evidence as a group of text files in a specified output directory" +
                        " on Spark's HDFS file system.",
        programGroup = StructuralVariantDiscoveryProgramGroup.class)
public class ComputeMultiIntervalQnameCountsSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    protected final Logger logger = LogManager.getLogger(this.getClass());
    @ArgumentCollection
    private final StructuralVariationDiscoveryArgumentCollection.FindBreakpointEvidenceSparkArgumentCollection params =
            new StructuralVariationDiscoveryArgumentCollection.FindBreakpointEvidenceSparkArgumentCollection();

    @Override
    protected void runTool(final JavaSparkContext ctx) {
        final SVReadFilter filter = new SVReadFilter(params);
        final ReadMetadata readMetadata = FindBreakpointEvidenceSpark.buildMetadata(params, getHeaderForReads(), getUnfilteredReads(), filter, logger);
        logger.info("Metadata retrieved.");

        // develop evidence, intervals, and, finally, a set of template names for each interval
        final FindBreakpointEvidenceSpark.EvidenceScanResults
                evidenceScanResults = FindBreakpointEvidenceSpark.getMappedQNamesSet(params, readMetadata, ctx,
                getHeaderForReads(), getUnfilteredReads(), filter, logger);
        final List<SVInterval> intervals = evidenceScanResults.intervals;

        final HopscotchUniqueMultiMap<String, Integer, QNameAndInterval> qNamesMultiMap = evidenceScanResults.qNamesForAssemblyMultiMap;

        final HopscotchUniqueMultiMap<String, Integer, QNameAndInterval> qNamesMultiMapMerged = FindBreakpointEvidenceSpark.mergeAssemblies(qNamesMultiMap, evidenceScanResults.intervals.size(), logger);

    }
}
