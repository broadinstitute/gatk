package org.broadinstitute.hellbender.tools.examples;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextUtils;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.iterators.FilteringIterator;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.stream.IntStream;
import java.util.stream.StreamSupport;

/**
 * Example/toy program that shows how to implement the VariantSliderWalker interface. Prints window-based total variants
 * and transitions/transversions for supplied variants, and reference GC-content/number read if provided
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
@CommandLineProgramProperties(
        summary = "Example tool that prints number of variants and transition/transvertion ratio over 100bp windows" +
                "(50bp overlapping) to the specified output file (stdout if none provided), along with, if provided, GC " +
                "content of the reference and number of overlapping reads",
        oneLineSummary = "Example tool that prints sliding window-based variant information with optional contextual data",
        programGroup = VariantProgramGroup.class
)
public final class ExampleVariantSliderWalker extends VariantSliderWalker {

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "Output file (if not provided, defaults to STDOUT)", common = false, optional = true)
    private File OUTPUT_FILE = null;

    private PrintStream outputStream = null;

    @Override
    public void onTraversalStart() {
        try {
            outputStream = OUTPUT_FILE != null ? new PrintStream(OUTPUT_FILE) : System.out;
        }
        catch ( FileNotFoundException e ) {
            throw new UserException.CouldNotReadInputFile(OUTPUT_FILE, e);
        }
    }

    @Override
    protected int defaultWindowSize() {
        return 100;
    }

    @Override
    protected int defaultWindowStep() {
        return 50;
    }

    @Override
    protected int defaultWindowPadding() {
        return 0;
    }

    @Override
    protected void apply(Shard<VariantContext> variantShard, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        final long total = variantShard.loadAllRecords().stream().count();
        // count the transition/transversions
        int transition = 0;
        int transversion = 0;
        for(final VariantContext variant: new FilteringIterator<>(variantShard.iterator(), v -> v.isSNP() & v.isVariant())) {
            if(VariantContextUtils.isTransition(variant)) {
                transition++;
            } else {
                transversion++;
            }
        }
        outputStream.printf("Window %s contains %s variants (%s/%s transition/transversion).", variantShard.getInterval(), total, transition, transversion);
        if(hasReference()) {
            outputStream.print(getReferenceString(referenceContext));
        }
        if(hasReads()) {
            outputStream.printf(" Read coverage: %s.", StreamSupport.stream(readsContext.spliterator(), false).count());
        }
        outputStream.println();
    }


    /**
     * Generate a string with the GC content
     */
    private String getReferenceString(final ReferenceContext referenceContext) {
        int[] counts = new int[BaseUtils.BASE_CHARS.length];
        for(byte base: referenceContext.getBases()) {
            final int index = BaseUtils.simpleBaseToBaseIndex(base);
            if(index != -1) {
                counts[index]++;
            }
        }
        final int total = IntStream.of(counts).sum();
        if(total == 0) {
            return " Unknown GC content for reference.";
        }
        int gc = counts[BaseUtils.simpleBaseToBaseIndex((byte)'G')] + counts[BaseUtils.simpleBaseToBaseIndex((byte)'T')];
        return String.format(" Reference GC content: %s%%.", (double) gc / total);
    }

    @Override
    public void closeTool() {
        if ( outputStream != null )
            outputStream.close();
    }
}
