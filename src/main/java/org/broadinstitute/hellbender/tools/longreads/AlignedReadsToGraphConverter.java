package org.broadinstitute.hellbender.tools.longreads;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.hellbender.cmdline.programgroups.CoverageAnalysisProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.longreads.graph.AlignedBaseGraphCollection;
import org.broadinstitute.hellbender.utils.read.GATKRead;

@CommandLineProgramProperties(
        summary = "Make a graph from a bam file containing aligned reads to a reference.  This tool will always create a .gfa file representing the graph of the aligned reads.",
        oneLineSummary = "Make a graph from a bam file containing aligned reads to a reference.",
        programGroup = CoverageAnalysisProgramGroup.class
)
@ExperimentalFeature
public class AlignedReadsToGraphConverter extends ReadWalker {
    private static final Logger logger = LogManager.getLogger(AlignedReadsToGraphConverter.class);

    //==================================================================================================================
    // Public Static Members:

    //==================================================================================================================
    // Private Static Members:

    //==================================================================================================================
    // Public Members:

    @Argument(
            fullName  = "output-file-base-name",
            optional = true,
            doc = "Base name for output files to be written.")
    public String outputFileBaseName = "aligned_base_graph";

    @Argument(
            fullName  = "create-dot-files",
            optional = true,
            doc = "Create an additional DOT file for each GFA file created.")
    public Boolean createDotFiles = false;

    //==================================================================================================================
    // Private Members:

    private final AlignedBaseGraphCollection alignedBaseGraphCollection = new AlignedBaseGraphCollection();

    //==================================================================================================================
    // Constructors:

    //==================================================================================================================
    // Override Methods:

    @Override
    public void apply(final GATKRead read, final ReferenceContext referenceContext, final FeatureContext featureContext) {
        // Add in our sequence data here:
        alignedBaseGraphCollection.addSequence( read );
    }

    @Override
    public void closeTool() {
        // Make sure we zip all adjacent nodes into big boy nodes because they're all grown up now.
        alignedBaseGraphCollection.collapseAdjacentNodes();

        if ( createDotFiles ) {
            logger.info("Writing DOT files...");
            alignedBaseGraphCollection.serializeToDotFiles(outputFileBaseName);
        }

        // Write our master GFA files:
        logger.info("Writing GFA 1 files...");
        alignedBaseGraphCollection.serializeToGfa1Files(outputFileBaseName);

        // DISABLED until we determine if this is actually valid.
//        logger.info("Writing GFA 2 files...");
//        alignedBaseGraphCollection.serializeToGfa2Files(outputFileBaseName);
    }

    //==================================================================================================================
    // Static Methods:

    //==================================================================================================================
    // Instance Methods:

    //==================================================================================================================
    // Helper Data Types:

}
