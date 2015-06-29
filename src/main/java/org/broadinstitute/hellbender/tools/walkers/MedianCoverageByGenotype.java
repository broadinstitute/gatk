package org.broadinstitute.hellbender.tools.walkers;

import htsjdk.samtools.util.Histogram;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeType;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
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
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

@CommandLineProgramProperties(
        usage = "Compute the median coverage for each genotype.",
        usageShort = "Compute the median coverage for each genotype.",
        programGroup = VariantProgramGroup.class
)
public final class MedianCoverageByGenotype extends VariantWalker {

    @Argument(fullName = "output", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc="Write output to this file")
    public File OUTPUT;

    private PrintWriter writer;

    @Override
    public void onTraversalStart() {
        try {
            writer = new PrintWriter(OUTPUT);
        } catch (final FileNotFoundException e) {
            throw new UserException.CouldNotCreateOutputFile(OUTPUT, e);
        }
        writer.println("CHR\tPOS\tID\tHOMREF\tHET\tHOMVAR");
    }

    @Override
    public void apply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        final Map<GenotypeType, Double> medianCoverageByGenotype = new HashMap<>();
         variant.getGenotypes().stream()
                .collect(Collectors.groupingBy(Genotype::getType))
                .forEach((gt, v) -> {
                    double[] dpFields = v.stream().map(g -> (g.getDP())).mapToDouble(i -> (double)i).toArray();
                    if (dpFields.length > 1) {
                        medianCoverageByGenotype.put(gt, new DescriptiveStatistics(dpFields).getPercentile(50));
                    }
                    else if (dpFields.length == 1){ //getPercentile returns NaN if length = 1
                        medianCoverageByGenotype.put(gt, dpFields[0]);
                    }
                    else {
                        medianCoverageByGenotype.put(gt, Double.NaN);
                    }
                });
        writer.print(variant.getChr());
        writer.print('\t');
        writer.print(variant.getStart());
        writer.print('\t');
        writer.print(variant.getID());
        writer.print('\t');
        writer.print(medianCoverageByGenotype.getOrDefault(GenotypeType.HOM_REF, Double.NaN));
        writer.print('\t');
        writer.print(medianCoverageByGenotype.getOrDefault(GenotypeType.HET, Double.NaN));
        writer.print('\t');
        writer.print(medianCoverageByGenotype.getOrDefault(GenotypeType.HOM_VAR, Double.NaN));
        writer.println();
    }

    @Override
    public Object onTraversalDone() {
        writer.close();
        return "WORKED!";
    }

}
