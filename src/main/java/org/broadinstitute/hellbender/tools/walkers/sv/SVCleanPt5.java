package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.samtools.util.OverlapDetector;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;

import java.util.Arrays;
import java.util.List;
import java.util.ArrayList;
import java.util.Set;
import java.util.Map;
import java.util.HashSet;
import java.util.HashMap;
import java.util.Collections;

/**
 * Completes an initial series of cleaning steps for a VCF produced by the GATK-SV pipeline.
 *
 * <h3>Inputs</h3>
 * <ul>
 *     <li>
 *         VCF containing structural variant (SV) records from the GATK-SV pipeline.
 *     </li>
 *     <li>
 *         TODO
 *     </li>
 * </ul>
 *
 * <h3>Output</h3>
 * <ul>
 *     <li>
 *         Cleansed VCF.
 *     </li>
 * </ul>
 *
 * <h3>Usage Example</h3>
 * <pre>
 *     gatk SVCleanPt5 \
 *       -V input.vcf.gz \
 *       --sample-list samples.txt \
 * 	     --multi-cnv-list multi.cnvs.txt
 * 	     --output-prefix result
 * </pre>
 *
 * <h3>Processing Steps</h3>
 * <ol>
 *     <li>
 *         TODO
 *     </li>
 * </ol>
 */
@CommandLineProgramProperties(
        summary = "Clean and format SV VCF",
        oneLineSummary = "Clean and format SV VCF",
        programGroup = StructuralVariantDiscoveryProgramGroup.class
)
@BetaFeature
@DocumentedFeature
public class SVCleanPt5 extends MultiplePassVariantWalker { // MultiVariantWalker?
    @Override
    protected int numberOfPasses() {
        return 2;
    }

    @Override
    public void onTraversalStart() {
        /*
        try {
            revisedCnWriter = Files.newBufferedWriter(Paths.get(outputPrefix + ".txt"));

            sampleWhitelist = new HashSet<>(Files.readAllLines(sampleListPath.toPath()));
            multiallelicCnvs = new HashSet<>(Files.readAllLines(multiCnvPath.toPath()));
        } catch (IOException e) {
            throw new RuntimeException("Error reading input file", e);
        }
         */
        return;
    }

    @Override
    public Object onTraversalSuccess() {
        /*
        try {
            List<String> variantIDs = new ArrayList<>(revisedCopyNumbers.keySet());
            Collections.sort(variantIDs);

            for (String variantID : variantIDs) {
                Map<String, Integer> sampleMap = revisedCopyNumbers.get(variantID);

                List<String> samples = new ArrayList<>(sampleMap.keySet());
                Collections.sort(samples);

                for (String sample : samples) {
                    int rdCn = sampleMap.get(sample);
                    revisedCnWriter.write(variantID + "\t" + sample + "\t" + rdCn);
                    revisedCnWriter.newLine();
                }
            }

            if (revisedCnWriter != null) {
                revisedCnWriter.close();
            }

            return null;
        } catch (IOException e) {
            throw new RuntimeException("Error writing multiallelic CNVs", e);
        }
         */
        return null;
    }

    @Override
    protected void nthPassApply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext, int n) {
        /*
        switch (n) {
            case 0:
                firstPassApply(variant);
                break;
            case 1:
                secondPassApply(variant);
                break;
            default:
                throw new IllegalArgumentException("Invalid pass number: " + n);
        }
         */
        return;
    }

    @Override
    protected void afterNthPass(int n) {
        return;
    }
}
