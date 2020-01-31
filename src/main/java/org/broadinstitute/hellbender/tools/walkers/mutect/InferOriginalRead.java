package org.broadinstitute.hellbender.tools.walkers.mutect;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.tools.walkers.mutect.consensus.DuplicateSet;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMFileGATKReadWriter;
import picard.cmdline.argumentcollections.ReferenceArgumentCollection;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;

import java.io.File;
import java.util.*;

import org.broadinstitute.hellbender.engine.filters.UMIReadFilter;

@CommandLineProgramProperties(
        summary = "Prints reads from the input SAM/BAM/CRAM file to the SAM/BAM/CRAM file.",
        oneLineSummary = "Print reads in the SAM/BAM/CRAM file",
        programGroup = ReadDataManipulationProgramGroup.class
)
public class InferOriginalRead extends DuplicateSetWalker {
    private SAMFileGATKReadWriter outputWriter;

    private static final int INIT_DUPLICATE_SIZE = 20;
    ArrayList<GATKRead> currentReads = new ArrayList<>(INIT_DUPLICATE_SIZE);
    String currentUMI = "";
    public static final String FGBIO_MOLECULAR_IDENTIFIER_TAG = "MI";
    // This is the value of the "MI" tag e.g. Z:20/A, where Z is __, 20 means indicates that the read belongs to
    // the 20th molecule in the bam, and A/B indicates Top/Bottom strand
    String currentMolecularIdentifier = "";

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "")
    public File outputBam;
    ReferenceContext lastSeenReferenceContext = null;

    DuplexConsensusCaller engine;

    @Override
    public boolean requiresReference(){ return true; }


    @Override
    public void onTraversalStart(){
        // check for UMI tag
        createOutputBamIndex = true;
        outputWriter = createSAMWriter(IOUtils.getPath(outputBam.getAbsolutePath()), false);
        engine = new InferOriginalReadEngine(getHeaderForReads(), referenceArguments, outputWriter);
    }

    @Override
    public void apply(DuplicateSet duplicateSet, ReferenceContext referenceContext, FeatureContext featureContext) {
        engine.letsDoIt(duplicateSet, referenceContext);
    }

    @Override
    public Object onTraversalSuccess(){
        outputWriter.close();
        return "SUCCESS";
    }
}