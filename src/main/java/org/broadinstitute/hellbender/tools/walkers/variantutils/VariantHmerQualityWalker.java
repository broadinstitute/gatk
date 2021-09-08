package org.broadinstitute.hellbender.tools.walkers.variantutils;

import htsjdk.tribble.Feature;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.codecs.table.TableFeature;
import picard.cmdline.programgroups.ReferenceProgramGroup;

import java.io.PrintStream;
import java.util.List;
import java.util.Objects;


@DocumentedFeature
@CommandLineProgramProperties(
        oneLineSummary = "Walk gvcf and track GQ at hmer sites.",
        summary = "Walker used to collect quality data from variant file along hmer sites.",
        programGroup = ReferenceProgramGroup.class
)
public class VariantHmerQualityWalker extends VariantWalker {

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "Output file (if not provided, defaults to STDOUT)", common = false, optional = true)
    private GATKPath outputFile = null;

    @Argument(fullName="auxiliaryVariants", shortName="av", doc="Auxiliary set of variants", optional=true)
    private FeatureInput<VariantContext> auxiliaryVariants;

    private PrintStream outputStream = null;

    @Argument(shortName = "F", fullName = "feature_file", doc = "Feature file (eg., VCF or BED file)")
    public FeatureInput<TableFeature> featuresFile = null;


    public String findGT(Genotype g) {
        if (g.isHomRef()) {
            return "HomRef";
        } else if (g.isHet()) {
            return "Het";
        } else if (g.isHomVar()) {
            return "HomVar";
        } else {
            return "Unknown";
        }
    }


    @Override
    public void onTraversalStart() {
        outputStream = outputFile != null ? new PrintStream(outputFile.getOutputStream()) : System.out;
        outputStream.println("HEADER GT GQ GCHmerContent LinHmerContent QuadHmerContent MaxHmerLength");
    }

    @Override
    public void apply( VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext ) {
        /** Create string for current variant interval */
        final String var_loc = String.format("%s:%s-%s", variant.getContig(), variant.getStart(), variant.getEnd());

        /** Get list of overlapping hmer intervals from table */
        final List<TableFeature> overlappingHmerFeatures = featureContext.getValues(featuresFile);

        /** Initialize hmer metrics to collect */
        int linearGCCount = 0;
        int linearHmerCount = 0;
        int quadHmerCount = 0;
        int maxHmerLength = 0;

        /** Loop over overlapping hmer intervals and update metrics */
        for (TableFeature hmerIntervalLine : overlappingHmerFeatures) {

            /** Get base and length from hmer table */
            String base = hmerIntervalLine.get("base");
            int length = Integer.parseInt(hmerIntervalLine.get("length"));

            /** Get max overlapping hmer length */
            if (length > maxHmerLength):
                maxHmerLength = length;

            /** Update counts accordingly */
            linearHmerCount += length;
            quadHmerCount += length * length;
            if ( (Objects.equals(base, "G")) || (Objects.equals(base, "C"))) {
                linearGCCount += length;
            }
        }

        /** Normalize metrics to get hmer content */
        final int intervalWidth = variant.getEnd() - variant.getStart()+1;
        final float linearHmerContent = (float)linearHmerCount/intervalWidth;
        final float quadHmerContent = (float)quadHmerCount/intervalWidth;
        final float linearGCContent = (float)linearGCCount/linearHmerContent;


        outputStream.printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
                var_loc, findGT(variant.getGenotype(0)), variant.getGenotype(0).getGQ(),
                linearGCContent, linearHmerContent, quadHmerContent, maxHmerLength);
    }


    @Override
    public void closeTool() {
        if ( outputStream != null ) {
            outputStream.close();
        }
    }
}
