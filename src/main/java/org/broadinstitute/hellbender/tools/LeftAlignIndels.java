package org.broadinstitute.hellbender.tools;


import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.WorkflowProperties;
import org.broadinstitute.barclay.argparser.WorkflowOutput;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.GATKPath;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.read.*;

/**
 * Left-aligns indels in read data
 *
 * <p>
 * This tool left-aligns any indels in the read data contained in a BAM or CRAM file. The same indel can often be
 * placed at multiple positions and still represent the same haplotype.  While it is a commonly used convention to place
 * an indel at the left-most position, this doesn't always happen (either because upstream tools broke ties between
 * equivalent representations randomly or used different aligning conventions), so this tool can be used to left-align
 * them according to convention. </p>
 *
 *<p><b>This tool will left-align reads with one and only one indel</b>.</p> 
 *
 * <h3>Input</h3>
 * <p>
 * A BAM or CRAM file to left-align.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * A left-aligned BAM or CRAM file.
 * </p>
 *
 * <h3>Usage example</h3>
 * <pre>
 * gatk LeftAlignIndels \
 *   -R reference.fasta \
 *   -I input.bam \
 *   -O output.bam
 * </pre>
 *
 */
@CommandLineProgramProperties(
        summary = "Left-aligns indels from reads in a SAM/BAM/CRAM file.",
        oneLineSummary = "Left-aligns indels from reads in a SAM/BAM/CRAM file",
        programGroup = ReadDataManipulationProgramGroup.class
)
@DocumentedFeature
@WorkflowProperties
public final class LeftAlignIndels extends ReadWalker {

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,doc="Output BAM")
    @WorkflowOutput(optionalCompanions={StandardArgumentDefinitions.OUTPUT_INDEX_COMPANION})
    private GATKPath output;

    private SAMFileGATKReadWriter outputWriter = null;

    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    public void onTraversalStart() {
        outputWriter = createSAMWriter(output, true);
    }

    @Override
    public void apply( GATKRead read, ReferenceContext ref, FeatureContext featureContext ) {
        // we can not deal with screwy records, and reads with a single cigar element are a trivial case
        if ( read.isUnmapped() || read.numCigarElements() <= 1 ) {
            outputWriter.addRead(read);
            return;
        }

        final CigarBuilder.Result result = AlignmentUtils.leftAlignIndels(read.getCigar(), ref.getBases(), read.getBases(), 0);

        read.setCigar(result.getCigar());
        if (result.getLeadingDeletionBasesRemoved() > 0) {
            read.setPosition(read.getContig(), read.getStart() + result.getLeadingDeletionBasesRemoved());
        }

        outputWriter.addRead(read);
    }

    @Override
    public void closeTool() {
        if ( outputWriter != null ) {
            outputWriter.close();
        }
    }
}
