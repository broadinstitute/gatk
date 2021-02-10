package org.broadinstitute.hellbender.tools;

import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.WorkflowProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.argumentcollections.OptionalTextOutputArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.CoverageAnalysisProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * Calculate and print to the standard output (and optionally a file) the overall number of bases in a SAM/BAM/CRAM file
 *
 * <h3>Input</h3>
 * <ul>
 *     <li> A single BAM file</li>
 * </ul>
 *
 * <h3>Example</h3>
 *
 * <pre>
 *   gatk CountBases \
 *     -I input_reads.bam
 * </pre>
 */
@DocumentedFeature
@CommandLineProgramProperties(
	summary = " Calculate and print to the standard output (and optionally a file) the overall number of bases in a SAM/BAM/CRAM file\n" +
            "\n" +
            " <h3>Input</h3>\n" +
            " <ul>\n" +
            "     <li> A single BAM file</li>\n" +
            " </ul>\n" +
            "\n" +
            " <h3>Example</h3>\n" +
            "\n" +
            " <pre>\n" +
            "   gatk CountBases \\\n" +
            "     -I input_reads.bam\n" +
            " </pre>",
	oneLineSummary = "Count bases in a SAM/BAM/CRAM file",
    programGroup = CoverageAnalysisProgramGroup.class
)
@WorkflowProperties
public final class CountBases extends ReadWalker {

    private long count = 0;

    @ArgumentCollection
    final public OptionalTextOutputArgumentCollection out = new OptionalTextOutputArgumentCollection();

    @Override
    public void apply( GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext ) {
        count += read.getLength();
    }

    @Override
    public Object onTraversalSuccess() {
        out.print(count);
        return count;
    }
}
