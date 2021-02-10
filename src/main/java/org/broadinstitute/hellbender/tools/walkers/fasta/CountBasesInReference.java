package org.broadinstitute.hellbender.tools.walkers.fasta;

import com.google.common.annotations.VisibleForTesting;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.argumentcollections.OptionalTextOutputArgumentCollection;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceWalker;
import picard.cmdline.programgroups.ReferenceProgramGroup;

/**
 * Counts the number of times each base occurs in a reference, and prints the counts to standard output
 * (and optionally to a file).
 *
 * <h3> Input </h3>
 * <ul>
 *     <li> A single reference file. </li>
 * </ul>
 *
 * <h3> Usage example: </h3>
 * <pre>
 *  gatk CountBasesInReference \
 *       -R reference.fasta
 * </pre>
 */
@DocumentedFeature
@CommandLineProgramProperties(
        oneLineSummary = "Count the numbers of each base in a reference file",
        summary = "Count the number of times each individual base occurs in a reference file and print to standard output " +
                "(and optionally to a file).",
        programGroup = ReferenceProgramGroup.class
)
public class CountBasesInReference extends ReferenceWalker {

    @ArgumentCollection
    final public OptionalTextOutputArgumentCollection out = new OptionalTextOutputArgumentCollection();

    @VisibleForTesting
    final long[] baseCounts = new long[256];

    @Override
    public void apply(ReferenceContext referenceContext, ReadsContext read, FeatureContext featureContext) {
        baseCounts[referenceContext.getBase()]++;
    }

    @Override
    public Object onTraversalSuccess() {
        final StringBuilder sb = new StringBuilder();
        for (int i = 0; i < baseCounts.length; i++) {
            final long count = baseCounts[i];
            if (count > 0) {
                sb.append((char) i).append(" : ").append(count).append("\n");
            }
        }
        out.print(sb);
        System.out.print(sb.toString());

        return 0;
    }
}
