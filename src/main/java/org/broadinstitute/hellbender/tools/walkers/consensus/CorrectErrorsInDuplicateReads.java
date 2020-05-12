package org.broadinstitute.hellbender.tools.walkers.consensus;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.DuplicateSetWalker;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.mutect.consensus.DuplicateSet;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMFileGATKReadWriter;
import picard.cmdline.programgroups.OtherProgramGroup;

import java.io.File;
import java.util.Collection;
import java.util.List;

@CommandLineProgramProperties(
        summary = "",
        oneLineSummary = "",
        programGroup = OtherProgramGroup.class
)
public class CorrectErrorsInDuplicateReads extends DuplicateSetWalker {

    private SAMFileGATKReadWriter outputWriter;

    // This is the value of the "MI" tag e.g. Z:20/A, where Z is __, 20 means indicates that the read belongs to
    // the 20th molecule in the bam, and A/B indicates Top/Bottom strand
    String currentMolecularIdentifier = "";

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "")
    public File outputBam;
    ReferenceContext lastSeenReferenceContext = null;

    ErrorCorrector engine;

    @Override
    public boolean requiresReference(){ return true; }


    @Override
    public void onTraversalStart(){
        // check for UMI tag
        createOutputBamIndex = true;
        outputWriter = createSAMWriter(IOUtils.getPath(outputBam.getAbsolutePath()), false);
        engine = new PileupBasedErrorCorrectionEngine(getHeaderForReads(), referenceArguments, outputWriter);
    }

    @Override
    public void apply(DuplicateSet duplicateSet, ReferenceContext referenceContext, FeatureContext featureContext) {
        if (engine.filterDuplicateSet(duplicateSet)){
            return;
        }

        final Collection<GATKRead> errorCorrectedRead = engine.correctErrorsAndUpdateQualities(duplicateSet, referenceContext);
        errorCorrectedRead.forEach(r -> outputWriter.addRead(r));
    }

    @Override
    public Object onTraversalSuccess(){
        outputWriter.close();
        return "SUCCESS";
    }

}
