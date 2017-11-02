package org.broadinstitute.hellbender.tools.examples;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.*;
import java.util.stream.IntStream;

/**
 * Example/toy program that shows how to implement the ReadSliderWalker interface. Prints window-based
 * coverage from supplied reads, and reference bases/overlapping variants if provided
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
@CommandLineProgramProperties(
        summary = "Example tool that prints coverage from supplied reads over 100bp windows (50bp overlapping) to the specified " +
                "output file (stdout if none provided), along with, if provided, base composition of the reference and " +
                "number of each type of features",
        oneLineSummary = "Example tool that prints sliding window-based coverage with optional contextual data",
        programGroup = ReadProgramGroup.class
)
public final class ExampleReadSliderWalker extends ReadSliderWalker {

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "Output file (if not provided, defaults to STDOUT)", common = false, optional = true)
    private File OUTPUT_FILE = null;

    @Argument(fullName = StandardArgumentDefinitions.VARIANT_LONG_NAME, shortName = StandardArgumentDefinitions.VARIANT_SHORT_NAME, doc = "One or more VCF files", optional = true)
    private List<FeatureInput<VariantContext>> variants;

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
    public void apply(Shard<GATKRead> readShard, ReferenceContext referenceContext, FeatureContext featureContext) {
        final int readCoverage = readShard.loadAllRecords().size();
        outputStream.printf("Window %s contains %s reads.", readShard.getInterval(), readCoverage);
        if(hasReference()) {
            outputStream.print(getReferenceString(referenceContext));
        }
        if(hasFeatures()) {
            outputStream.print(getFeaturesString(featureContext));
        }
        outputStream.println();
    }

    /**
     * Generate a string with the number of variants of each type
     */
    private String getFeaturesString(final FeatureContext featureContext) {
        final Map<VariantContext.Type, Integer> map = new TreeMap<>();
        for(VariantContext var: featureContext.getValues(variants)) {
            map.compute(var.getType(), (k, v) -> v == null ? 1 : v + 1);
        }
        if(map.isEmpty()) {
            return " No overlapping variants.";
        }
        return String.format(" Overlapping variants: %s.", map.toString());
    }

    /**
     * Generate a string with the count of bases
     */
    private String getReferenceString(final ReferenceContext referenceContext) {
        int[] counts = new int[BaseUtils.BASE_CHARS.length];
        for(byte base: referenceContext.getBases()) {
            final int index = BaseUtils.simpleBaseToBaseIndex(base);
            if(index != -1) {
                counts[index]++;
            }
        }
        if(IntStream.of(counts).sum() == 0) {
            return " No ACGT reference bases.";
        }
        return String.format(" Reference bases %s:%s.", Arrays.toString(BaseUtils.BASE_CHARS), Arrays.toString(counts));
    }

    @Override
    public void closeTool() {
        if ( outputStream != null )
            outputStream.close();
    }
}
