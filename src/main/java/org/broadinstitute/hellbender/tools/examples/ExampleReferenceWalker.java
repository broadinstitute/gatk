package org.broadinstitute.hellbender.tools.examples;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ExampleProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.Comparator;
import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;

/**
 * Counts the number of times each reference context is seen as well as how many times it's overlapped by reads and variants.
 * Why you would want to do this is a mystery.
 *
 * This is a toy example to illustrate how to implement a ReferenceWalker.
 */
@CommandLineProgramProperties(
        summary = "Example of how to implement a ReferenceWalker that uses reads and features as well as a custom window",
        oneLineSummary = "Example of how to implement a ReferenceWalker",
        programGroup = ExampleProgramGroup.class,
        omitFromCommandLine = true)
public class ExampleReferenceWalker extends ReferenceWalker {

    @Argument(fullName = StandardArgumentDefinitions.VARIANT_LONG_NAME,
            shortName = StandardArgumentDefinitions.VARIANT_SHORT_NAME,
            doc="variants to count overlaps of")
    private FeatureInput<VariantContext> variants;

    @VisibleForTesting
    final Map<String, OverlapCounts> contextCounts = new HashMap<>();

    @Override
    protected SimpleInterval getReferenceWindow(SimpleInterval locus) {
        // add one base of padding to each side of each reference locus.
        return new SimpleInterval(locus.getContig(), Math.max(locus.getStart() -1, 1), locus.getStart()+1);
    }

    @Override
    public void apply(ReferenceContext referenceContext, ReadsContext readsContext, FeatureContext featureContext) {
        final byte[] bases = referenceContext.getBases();
        final String baseString = new String(bases);
        final OverlapCounts counts = contextCounts.getOrDefault(baseString, new OverlapCounts());
        if(readsContext.iterator().hasNext()){
            counts.overlappedByReads++;
        }

        if(!featureContext.getValues(variants).isEmpty()){
            counts.overlappedByVariants++;
        }

        counts.timesSeen++;

        contextCounts.put(baseString, counts);
    }

    @Override
    public Object onTraversalSuccess(){
        contextCounts.entrySet().stream()
                .sorted(Comparator.comparing(Entry::getKey))
                .forEachOrdered(entry -> System.out.println(entry.getKey() + " " + entry.getValue().toString()));
        return 0;
    }

    @VisibleForTesting
    static class OverlapCounts {
        long timesSeen = 0;
        long overlappedByReads = 0;
        long overlappedByVariants = 0;

        @Override
        public String toString(){
            return String.format("Seen: %d Reads %d Variants %d", timesSeen, overlappedByReads, overlappedByVariants);
        }
    }
}
