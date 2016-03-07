package org.broadinstitute.hellbender.tools.walkers.variantutils;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;
import org.broadinstitute.hellbender.engine.filters.VariantFilter;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.function.Function;
import java.util.function.Predicate;
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

    protected VariantFilter makeVariantFilter() {
        return vc -> vc.isBiallelic(); //Note: this restriction should be removed in the future
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
        //Note: for now, this is very slow because the same work is redone over and over.
        // We'll optimize once all functionality is in place.
        final Map<String, Function<VariantContext, String>> results = new LinkedHashMap<>();
        results.put("callRate",           vc -> formatDouble(callRate(vc)));
        results.put("nCalled",            vc -> formatInt(countGenotypes(vc, Genotype::isCalled)));
        results.put("nNotCalled",         vc -> formatInt(countGenotypes(vc, gt -> !gt.isCalled())));
        results.put("nHomRef",            vc -> formatInt(countGenotypes(vc, Genotype::isHomRef)));
        results.put("nHet",               vc -> formatInt(countGenotypes(vc, Genotype::isHet)));
        results.put("nHomVar",            vc -> formatInt(countGenotypes(vc, Genotype::isHomVar)));
        results.put("nNonRef",            vc -> formatInt(countGenotypes(vc, Genotype::isHet) + countGenotypes(vc, Genotype::isHomVar)));
        results.put("rHeterozygosity",    vc -> formatDouble(safeDivide(countGenotypes(vc, Genotype::isHet), countGenotypes(vc, Genotype::isCalled))));
        results.put("rHetHomVar",         vc -> formatDouble(safeDivide(countGenotypes(vc, Genotype::isHet),countGenotypes(vc, Genotype::isHomVar))));
        results.put("MAC",        vc -> formatInt(countGenotypes(vc, Genotype::isHet) +
                                                2*countGenotypes(vc, Genotype::isHomVar)));
        results.put("MAF",        vc -> {
                                    final int hets = countGenotypes(vc, Genotype::isHet);
                                    final int refs = 2*countGenotypes(vc, Genotype::isHomRef) + hets;
                                    final int alts = 2*countGenotypes(vc, Genotype::isHomVar) + hets;
                                    return formatDouble(safeDivide(alts, refs+alts));
                                    }
                );
        return results;
    }

    private static double safeDivide(final int d1, final int d2){
        if (d2 == 0){
            return Double.NaN;
        } else {
            return (double) d1 / d2;
        }
    }

    private static int countGenotypes(final VariantContext vc, final Predicate<Genotype> f) {
        return (int) vc.getGenotypes().stream().filter(f).count();
    }

    private static double callRate(final VariantContext vc) {
        return safeDivide(countGenotypes(vc, Genotype::isCalled), vc.getGenotypes().size());
    }

    private static String formatInt(final int n){
        return String.format("%d", n);
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
