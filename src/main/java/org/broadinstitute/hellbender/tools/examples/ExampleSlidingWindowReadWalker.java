package org.broadinstitute.hellbender.tools.examples;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.OptionalShardingArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.ShardingArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.ExampleProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.engine.SlidingWindowReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.Shard;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.IOException;
import java.io.PrintStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

/**
 * Example/toy program that shows how to implement the {@link SlidingWindowReadWalker} interface. Prints window-based
 * coverage from supplied reads, and reference bases/overlapping variants if provided.
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
@CommandLineProgramProperties(
        summary = "Example tool that prints coverage from supplied reads over 100bp windows by default (50bp overlapping) to the specified " +
                "output file (stdout if none provided), along with, if provided, base composition of the reference and " +
                "number of each type of features",
        oneLineSummary = "Example tool that prints sliding window-based coverage with optional contextual data",
        programGroup = ExampleProgramGroup.class
)
public final class ExampleSlidingWindowReadWalker extends SlidingWindowReadWalker {

    public static final int DEFAULT_WINDOW_SIZE = 100;
    public static final int DEFAULT_WINDOW_STEP = 50;
    public static final int DEFAULT_WINDOW_PAD = 0;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "Output file (if not provided, defaults to STDOUT)", common = false, optional = true)
    private String outputFile = null;

    @Argument(fullName = StandardArgumentDefinitions.VARIANT_LONG_NAME, shortName = StandardArgumentDefinitions.VARIANT_SHORT_NAME, doc = "One or more VCF files", optional = true)
    private List<FeatureInput<VariantContext>> variants;

    private PrintStream outputStream = System.out;

    /** Initialize the output if provided. */
    @Override
    public void onTraversalStart() {
        // initialize the output file if present
        if (outputFile != null) {
            try {
                final Path path = IOUtils.getPath(outputFile);
                outputStream = new PrintStream(Files.newOutputStream(path));
            } catch (final IOException e) {
                throw new UserException.CouldNotCreateOutputFile(outputFile, e);
            }
        }
    }

    /** Returns an optional argument collection. */
    @Override
    public ShardingArgumentCollection getShardingArgs() {
        // fixed sliding window parameters
        return new OptionalShardingArgumentCollection(DEFAULT_WINDOW_SIZE, DEFAULT_WINDOW_STEP, DEFAULT_WINDOW_PAD);
    }

    @Override
    public void apply(final Shard<GATKRead> reads, final ReferenceContext reference, final FeatureContext features) {
        // compute read coverage over window
        final long readCoverage = Utils.stream(reads.iterator()).count();

        // output coverage  over window
        outputStream.printf("%s (%s padded): %d reads.",
                reads.getInterval(),
                reads.getPaddedInterval(),
                readCoverage);

        // add reference if there is data source for it
        if(reference.hasBackingDataSource()) {
            outputStream.print(getReferenceString(reference));
        }
        // add features if there is data source for it
        if(features.hasBackingDataSource()) {
            outputStream.print(getFeaturesString(features));
        }
        // print a new line
        outputStream.println();
    }

    /**
     * Generate a string with the number of variants of each type
     */
    private String getFeaturesString(final FeatureContext featureContext) {
        // initialize map
        final Map<VariantContext.Type, Integer> map = new TreeMap<>();
        // add the variants
        for(VariantContext var: featureContext.getValues(variants)) {
            map.compute(var.getType(), (k, v) -> v == null ? 1 : v + 1);
        }
        // return result
        return String.format(" Variants: %s.", map.toString());
    }

    /**
     * Generate a string with the count of bases
     */
    private String getReferenceString(final ReferenceContext referenceContext) {
        // initialize the map
        final Map<Character, Integer> baseCounts = new LinkedHashMap<>(BaseUtils.BASE_CHARS.length);
        for (final char c: BaseUtils.BASE_CHARS) {
            baseCounts.put(c, 0);
        }
        // count bases
        for(byte base: referenceContext.getBases()) {
            if (BaseUtils.isRegularBase(base)) {
                baseCounts.compute((char) base, (c, i) -> i + 1);
            }
        }
        // return result
        return String.format(" ACGT reference bases %s.", baseCounts);
    }

    @Override
    public void closeTool() {
        // close the output stream if present
        if ( outputStream != null ) {
            outputStream.close();
        }
    }
}
