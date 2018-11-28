package org.broadinstitute.hellbender.tools.walkers.fasta;

import com.google.common.annotations.VisibleForTesting;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceWalker;
import picard.cmdline.programgroups.ReferenceProgramGroup;

@DocumentedFeature
@CommandLineProgramProperties(
  summary = "Count the numbers of each base in a reference file",
  oneLineSummary = "Count the number of each individual base which occurs in a reference file.",
  programGroup = ReferenceProgramGroup.class
)
public class CountBasesInReference extends ReferenceWalker {
    @VisibleForTesting
    final long[] baseCounts = new long[256];

    @Override
    public void apply(ReferenceContext referenceContext, ReadsContext read, FeatureContext featureContext) {
      baseCounts[referenceContext.getBase()]++;
    }

    @Override
    public Object onTraversalSuccess(){
        for (int i = 0; i < 256; i++) {
            final long count = baseCounts[i];
            if (count > 0){
                System.out.println((char)i + " : " + count );
            }
        }
        return 0;
    }
}
