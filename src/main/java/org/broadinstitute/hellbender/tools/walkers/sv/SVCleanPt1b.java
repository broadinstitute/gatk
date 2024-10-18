package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

import java.util.Arrays;
import java.util.List;
import java.util.ArrayList;
import java.util.Set;
import java.util.Map;
import java.util.HashSet;
import java.util.HashMap;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;

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
 *     gatk SVCleanPt1b \
 *       -V input.vcf.gz \
 *       -O output.vcf.gz
 * </pre>
 *
 * <h3>Cleaning Steps</h3>
 * <ol>
 *     <li>
 *         TODO
 *     </li>
 * </ol>
 */
@CommandLineProgramProperties(
        summary = "Clean and format structural variant VCFs per Step 1b",
        oneLineSummary = "Clean and format structural variant VCFs per Step 1b",
        programGroup = StructuralVariantDiscoveryProgramGroup.class
)
@BetaFeature
@DocumentedFeature
public class SVCleanPt1b extends MultiplePassVariantWalker {
    public static final String CNV_FILE_LONG_NAME = "cnv-file";

    @Argument(
            fullName = CNV_FILE_LONG_NAME,
            doc = "Output CNVs file name",
            optional = true
    )
    private final GATKPath outputCnvs = new GATKPath(GATKSVVCFConstants.CNVS_DEFAULT_FILE);

    @Argument(
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "Output VCF name"
    )
    private GATKPath outputVcf;

    private VariantContextWriter vcfWriter;
    private BufferedWriter cnvsWriter;

    private final List<VariantContext> overlappingVariantsBuffer = new ArrayList<>();
    final private Set<String> multiCnvs = new HashSet<>();
    final private Map<String, Map<String, Pair<String, String>>> revisedEventsAll = new HashMap<>();
    final private Map<String, Set<String>> revisedEventsFiltered = new HashMap<>();
    final private Map<String, Map<String, Integer>> revisedRdCn = new HashMap<>();

    @Override
    protected int numberOfPasses() {
        return 3;
    }

    @Override
    public Object onTraversalSuccess() {
        try {
            cnvsWriter = new BufferedWriter(new FileWriter(outputCnvs.toPath().toFile()));
            for (String variantId : multiCnvs) {
                cnvsWriter.write(variantId);
                cnvsWriter.newLine();
            }
            return null;
        } catch (IOException e) {
            throw new RuntimeException("Error creating CNVs file", e);
        }
    }

    @Override
    public void closeTool() {
        try {
            if (vcfWriter != null) {
                vcfWriter.close();
            }
            if (cnvsWriter != null) {
                cnvsWriter.close();
            }
        } catch (IOException e) {
            throw new RuntimeException("Error closing output file", e);
        }
    }

    @Override
    protected void nthPassApply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext, int n) {
        switch (n) {
            case 0:
                firstPassApply(variant);
                break;
            case 1:
                secondPassApply(variant);
                break;
            case 2:
                thirdPassApply(variant);
                break;
        }
    }

    @Override
    protected void afterNthPass(int n) {
        switch (n) {
            case 0:
                processCollectedVariants();
                break;
            case 1:
                vcfWriter = createVCFWriter(outputVcf);
                vcfWriter.writeHeader(getHeaderForVariants());
                break;
        }
    }

    public void firstPassApply(final VariantContext variant) {
        final String svType = variant.getAttributeAsString(GATKSVVCFConstants.SVTYPE, "");
        final boolean isDelDup = svType.equals(GATKSVVCFConstants.SYMB_ALT_STRING_DUP) || svType.equals(GATKSVVCFConstants.SYMB_ALT_STRING_DEL);
        final boolean isLarge = Math.abs(variant.getAttributeAsInt(GATKSVVCFConstants.SVLEN, 0)) >= 5000;
        if (isDelDup && isLarge) {
            overlappingVariantsBuffer.removeIf(vc -> vc.getEnd() < variant.getStart());
            for (VariantContext bufferedVariant : overlappingVariantsBuffer) {
                if (overlaps(bufferedVariant, variant)) {
                    processOverlap(bufferedVariant, variant);
                }
            }
            overlappingVariantsBuffer.add(variant);
        }
    }

    public void secondPassApply(final VariantContext variant) {
        if (revisedEventsFiltered.containsKey(variant.getID())) {
            initializeRdCn(variant);
        }
    }

    public void thirdPassApply(final VariantContext variant) {
        VariantContextBuilder builder = new VariantContextBuilder(variant);
        if (revisedEventsAll.containsKey(variant.getID())) {
            processVariant(builder, variant);
        }

        final String svType = variant.getAttributeAsString(GATKSVVCFConstants.SVTYPE, "");
        final boolean isDelDup = svType.equals(GATKSVVCFConstants.SYMB_ALT_STRING_DUP) || svType.equals(GATKSVVCFConstants.SYMB_ALT_STRING_DEL);
        final boolean isLarge = variant.getEnd() - variant.getStart() >= 1000;
        if (isDelDup && isLarge) {
            processCnvs(variant);
        }
        vcfWriter.add(builder.make());
    }

    private void processOverlap(VariantContext v1, VariantContext v2) {
        // Get overlap data
        VariantContext wider;
        VariantContext narrower;
        if (v1.getLengthOnReference() > v2.getLengthOnReference()) {
            wider = v1;
            narrower = v2;
        } else if (v2.getLengthOnReference() > v1.getLengthOnReference()) {
            wider = v2;
            narrower = v1;
        } else {
            return;
        }
        String widerID = wider.getID();
        String narrowerID = narrower.getID();

        // Skip processing if same variant ID, SV type or samples
        String widerSvType = wider.getAttributeAsString(GATKSVVCFConstants.SVTYPE, "");
        String narrowerSvType = narrower.getAttributeAsString(GATKSVVCFConstants.SVTYPE, "");
        Set<String> widerSamples = getNonReferenceSamples(wider);
        Set<String> narrowerSamples = getNonReferenceSamples(narrower);
        if (widerID.equals(narrowerID) || widerSvType.equals(narrowerSvType) || widerSamples.equals(narrowerSamples)) {
            return;
        }

        // Get samples present in wider but not in narrower
        Set<String> nonCommonSamples = new HashSet<>(widerSamples);
        nonCommonSamples.removeAll(narrowerSamples);
        if (nonCommonSamples.isEmpty()) {
            return;
        }

        // Revise variant if coverage exceeds threshold
        double coverage = getCoverage(wider, narrower);
        if (coverage >= 0.5) {
            for (String sample : nonCommonSamples) {
                revisedEventsAll.computeIfAbsent(narrowerID, k -> new HashMap<>())
                        .put(sample, new ImmutablePair<>(widerID, widerSvType));
            }
        }
    }

    private void processCollectedVariants() {
        for (final Map.Entry<String, Map<String, Pair<String, String>>> entry : revisedEventsAll.entrySet()) {
            for (final Map.Entry<String, Pair<String, String>> innerEntry : entry.getValue().entrySet()) {
                // Identifies variant-sample pairs we need RD_CN values for to improve speed
                final String sampleName = innerEntry.getKey();
                final String variantId = entry.getKey();
                final String widerVariantId = innerEntry.getValue().getLeft();
                final String svType = innerEntry.getValue().getRight();
                if (svType.equals(GATKSVVCFConstants.SYMB_ALT_STRING_DUP) || svType.equals(GATKSVVCFConstants.SYMB_ALT_STRING_DEL)) {
                    revisedEventsFiltered.computeIfAbsent(variantId, k -> new HashSet<>()).add(sampleName);
                    revisedEventsFiltered.computeIfAbsent(widerVariantId, k -> new HashSet<>()).add(sampleName);
                }
            }
        }
    }

    private void initializeRdCn(final VariantContext variant) {
        // Initialize data structures
        final String variantId = variant.getID();
        final Set<String> samples = revisedEventsFiltered.get(variantId);
        final Map<String, Integer> variantRdCn = new HashMap<>();

        // Initialize revisedRdCn value for each variant
        for (final String sampleName : samples) {
            final Genotype genotype = variant.getGenotype(sampleName);
            final String rdCn = (String) genotype.getExtendedAttribute(GATKSVVCFConstants.RD_CN);
            variantRdCn.put(sampleName, Integer.parseInt(rdCn));
        }
        revisedRdCn.put(variantId, variantRdCn);
    }

    private void processVariant(final VariantContextBuilder builder, final VariantContext variant) {
        // Initialize data structures
        final String variantId = variant.getID();
        final Map<String, Pair<String, String>> variantEvents = revisedEventsAll.get(variantId);
        final List<Genotype> newGenotypes = new ArrayList<>();

        // Create updated genotypes
        for (String sample : variant.getSampleNamesOrderedByName()) {
            final Genotype oldGenotype = variant.getGenotype(sample);
            final Pair<String, String> event = variantEvents.get(sample);

            if (event != null) {
                final String widerVariantId = event.getLeft();
                final String widerSvType = event.getRight();
                final int currentRdCn = revisedRdCn.get(variantId).getOrDefault(sample, 0);
                final int widerRdCn = revisedRdCn.getOrDefault(widerVariantId, new HashMap<>()).getOrDefault(sample, 0);
                if (!revisedEventsFiltered.getOrDefault(widerVariantId, new HashSet<>()).contains(sample)) {
                    System.err.println(sample + " " + widerVariantId);
                }

                int newVal = -1;
                if (widerSvType.equals(GATKSVVCFConstants.SYMB_ALT_STRING_DUP) && currentRdCn == 2 && widerRdCn == 3) {
                    newVal = 1;
                } else if (widerSvType.equals(GATKSVVCFConstants.SYMB_ALT_STRING_DEL) && currentRdCn == 2 && widerRdCn == 1) {
                    newVal = 3;
                }

                if (newVal != -1) {
                    final GenotypeBuilder gb = new GenotypeBuilder(oldGenotype);
                    gb.alleles(Arrays.asList(variant.getReference(), variant.getAlternateAllele(0)));
                    gb.GQ(Integer.parseInt((String) oldGenotype.getExtendedAttribute(GATKSVVCFConstants.RD_GQ)));
                    newGenotypes.add(gb.make());
                } else {
                    newGenotypes.add(oldGenotype);
                }
            } else {
                newGenotypes.add(oldGenotype);
            }
        }
        builder.genotypes(newGenotypes);
    }

    private void processCnvs(VariantContext variant) {
        final boolean isDel = variant.getAttributeAsString(GATKSVVCFConstants.SVTYPE, "").equals(GATKSVVCFConstants.SYMB_ALT_STRING_DEL);
        for (String sample : variant.getSampleNamesOrderedByName()) {
            final Genotype genotype = variant.getGenotype(sample);
            final String rdCnString = (String) genotype.getExtendedAttribute(GATKSVVCFConstants.RD_CN);
            final int rdCn = Integer.parseInt(rdCnString);
            if ((isDel && rdCn > 3) || (!isDel && (rdCn < 1 || rdCn > 4))) {
                multiCnvs.add(variant.getID());
                break;
            }
        }
    }

    private boolean overlaps(VariantContext v1, VariantContext v2) {
        return v1.getContig().equals(v2.getContig()) && v1.getStart() <= v2.getEnd() && v2.getStart() <= v1.getEnd();
    }

    private Set<String> getNonReferenceSamples(VariantContext variant) {
        Set<String> samples = new HashSet<>();
        for (String sampleName : variant.getSampleNames()) {
            Genotype genotype = variant.getGenotype(sampleName);
            if (genotype.isCalled() && !genotype.isHomRef()) {
                samples.add(sampleName);
            }
        }
        return samples;
    }

    private double getCoverage(VariantContext wider, VariantContext narrower) {
        int nStart = narrower.getStart();
        int nStop = narrower.getEnd();
        int wStart = wider.getStart();
        int wStop = wider.getEnd();

        if (wStart <= nStop && nStart <= wStop) {
            int intersectionSize = Math.min(nStop, wStop) - Math.max(nStart, wStart) + 1;
            return (double) intersectionSize / (nStop - nStart + 1);
        }
        return 0.0;
    }
}
