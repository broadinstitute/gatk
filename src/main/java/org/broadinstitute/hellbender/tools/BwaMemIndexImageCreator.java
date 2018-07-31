package org.broadinstitute.hellbender.tools;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.programgroups.ReferenceProgramGroup;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndex;
/**
 * Create a BWA-MEM index image file for use with GATK BWA tools
 *
 * <p>Tools that utilize BWA-MEM (e.g. BwaSpark, PathSeqBwaSpark) require an index image file of the
 * reference sequences. This tool generates the image file from a reference FASTA file.</p>
 *
 * <h3>Input</h3>
 *
 * <ul>
 *     <li>Reference FASTA file</li>
 * </ul>
 *
 * <h4>Output</h4>
 *
 * <ul>
 *     <li>BWA-MEM index image file of the reference</li>
 * </ul>
 *
 * <h3>Usage example</h3>
 *
 * <pre>
 * gatk BwaMemIndexImageCreator \
 *     -I reference.fasta \
 *     -O reference.fasta.img
 * </pre>
 *
 */
@DocumentedFeature
@CommandLineProgramProperties(
        summary = "Create a BWA-MEM index image file for use with GATK BWA tools",
        oneLineSummary = "Create a BWA-MEM index image file for use with GATK BWA tools",
        programGroup = ReferenceProgramGroup.class
)
public final class BwaMemIndexImageCreator extends CommandLineProgram {

    @Argument(fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME,
            doc = "Input reference FASTA file location.")
    private String referenceFastaLoc = null;

    /**
     * If not provided, the default image file path will be the same as the reference FASTA with the extension ".img".
     */
    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "Output reference index image file (ending in \".img\").",
            optional = true)
    private String referenceIndexImageOutputLoc = null;

    @Override
    protected final Object doWork() {
        if (referenceIndexImageOutputLoc == null) {
            referenceIndexImageOutputLoc = referenceFastaLoc + ".img";
        }
        BwaMemIndex.createIndexImageFromFastaFile(referenceFastaLoc, referenceIndexImageOutputLoc);
        return null;
    }
}
