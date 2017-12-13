package org.broadinstitute.hellbender.tools;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMFileGATKReadWriter;

/**
 * Copy reads in a SAM/BAM/CRAM file verbatim to a SAM file
 *
 * <p>
 *     This tool copies reads from a SAM/BAM/CRAM file to a SAM file. If you are interested in printing reads in
 *     certain intervals only, use the -L argument.
 *
 * </p>
 *
 * <h3> Input </h3>
 * <ul>
 *     <li> SAM/BAM/CRAM file to be copied to a SAM file </li>
 *     <li> (optional) interval file with the -L argument </li>
 * </ul>
 *
 * <h3> Output </h3>
 * SAM file with reads copied from the source
 *
 *
 * <h3>Examples</h3>
 * <pre>
 *  gatk PrintReads \
 *  -R reference.fasta \
 *  -I input.bam \
 *  -L 20 \
 *  -O chr20.sam
 * </pre>
 *
 *
 * {@GATK.walkertype ReadWalker}
 */
@CommandLineProgramProperties(
	summary = "Prints reads from the input SAM/BAM/CRAM file to the SAM/BAM/CRAM file.",
    oneLineSummary = "Print reads in the SAM/BAM/CRAM file",
    programGroup = ReadProgramGroup.class
)
@DocumentedFeature
public final class PrintReads extends ReadWalker {

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="Write output to this file")
    public String output;

    private SAMFileGATKReadWriter outputWriter;

    @Override
    public void onTraversalStart() {
        outputWriter = createSAMWriter(IOUtils.getPath(output), true);
    }

    @Override
    public void apply( GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext ) {
        outputWriter.addRead(read);
    }

    @Override
    public void closeTool() {
        if ( outputWriter != null ) {
            outputWriter.close();
        }
    }
}
