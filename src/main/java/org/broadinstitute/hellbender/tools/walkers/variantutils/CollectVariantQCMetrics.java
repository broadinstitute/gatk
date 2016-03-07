package org.broadinstitute.hellbender.tools.walkers.variantutils;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.function.Function;
import java.util.stream.Collectors;

@CommandLineProgramProperties(
        summary = "VariantQC",
        oneLineSummary = "VariantQC",
        programGroup = VariantProgramGroup.class
)
public final class CollectVariantQCMetrics extends VariantWalker{

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="File to which variants should be written")
    public File out = null;

    private PrintStream outputStream = null;

    private Map<String, Function<VariantContext, String>> metricComputers;//Note: it's a list because the order is defined

    @Override
    public void onTraversalStart() {
        outputStream = initializeOutputStream();
        metricComputers = initializeMetrics();
        printHeader();
    }

    private void printHeader() {
        final String header = String.join("\t", metricComputers.keySet());
        outputStream.println(
                "contig" + "\t" +
                "start" + "\t" +
                "end" + "\t" +
                "ref" + "\t" +
                "alt" + "\t" +

                header);
    }

    private Map<String, Function<VariantContext, String>> initializeMetrics() {
        final Map<String, Function<VariantContext, String>> results = new LinkedHashMap<>();
        results.put("callRate", vc -> formatDouble(callRate(vc)));
        return results;
    }

    private double callRate(final VariantContext vc) {
        final int nGenotypes = vc.getGenotypes().size();
        return (double) vc.getGenotypes().stream().filter(gt -> gt.isCalled()).count() / nGenotypes;
    }

    private static String formatDouble(final double d){
        return String.format("%.3f", d);
    }

    private PrintStream initializeOutputStream() {
        try {
            return new PrintStream(out);
        } catch ( final FileNotFoundException e ) {
            throw new UserException.CouldNotCreateOutputFile(out, e);
        }
    }

    @Override
    public void apply(final VariantContext variant, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext) {
        final String line = metricComputers.values().stream().map(computer -> computer.apply(variant)).collect(Collectors.joining("\t"));
        outputStream.println(variant.getContig() + "\t" +
                             variant.getStart() + "\t" +
                             variant.getEnd()  + "\t" +
                             variant.getReference()  + "\t" +
                             variant.getAlternateAlleles()  + "\t" +
                             line);
    }

    /**
     * Close out the new variants file.
     */
    @Override
    public Object onTraversalDone() {
        try {
            return null;
        } finally {
            outputStream.close();
        }
    }
}
