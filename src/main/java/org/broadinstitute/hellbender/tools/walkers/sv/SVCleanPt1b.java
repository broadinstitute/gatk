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

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.util.zip.GZIPInputStream;
import java.nio.file.Files;
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
 *         Annotated VCF.
 *     </li>
 * </ul>
 *
 * <h3>Usage Example</h3>
 * <pre>
 *     gatk SVCleanPt1b \
 *       -V structural.vcf.gz \
 *       -O cleansed.vcf.gz
 *       --bed-file overlap.bed
 * </pre>
 *
 * <h3>Cleaning Steps</h3>
 * <ol>
 *     <li>
 *         Calculates new copy numbers for variant genotypes that match an overlapping variant.
 *     </li>
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
    public static final String BED_FILE_LONG_NAME = "bed-file";
    public static final String CNV_FILE_LONG_NAME = "cnv-file";

    @Argument(
            fullName = BED_FILE_LONG_NAME,
            doc = "BED file"
    )
    private GATKPath bedFile;

    @Argument(
            fullName = CNV_FILE_LONG_NAME,
            doc = "Output CNVs file name",
            optional = true
    )
    private GATKPath outputCnvs = new GATKPath(GATKSVVCFConstants.CNVS_DEFAULT_FILE);

    @Argument(
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "Output VCF name"
    )
    private GATKPath outputVcf;

    private VariantContextWriter vcfWriter;
    private BufferedWriter cnvsWriter;

    final private Set<String> multiCnvs = new HashSet<>();
    final private Map<String, Map<String, Pair<String, String>>> revisedEventsAll = new HashMap<>();
    final private Map<String, Set<String>> revisedEventsFiltered = new HashMap<>();
    final private Map<String, Map<String, Integer>> revisedRdCn = new HashMap<>();

    @Override
    protected int numberOfPasses() {
        return 3;
    }

    @Override
    protected void afterNthPass(int n) {
        return;
    }

    @Override
    public void onTraversalStart() {
        // Pre-process BED file
        processBedFile();

        // Write header
        vcfWriter = createVCFWriter(outputVcf);
        vcfWriter.writeHeader(getHeaderForVariants());
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
        switch(n) {
            case 1:
                firstPassApply(variant, readsContext, referenceContext, featureContext);
                break;
            case 2:
                secondPassApply(variant, readsContext, referenceContext, featureContext);
                break;
            case 3:
                thirdPassApply(variant, readsContext, referenceContext, featureContext);
                break;
        }
    }

    public void firstPassApply(final VariantContext variant, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext) {
        return;
    }

    public void secondPassApply(final VariantContext variant, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext) {
        if (shouldInitializeRdCn(variant)) {
            initializeRdCn(variant);
        }
    }

    public void thirdPassApply(final VariantContext variant, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext) {
        VariantContextBuilder builder = new VariantContextBuilder(variant);
        if (shouldProcessVariant(variant)) {
            processVariant(builder, variant);
        }
        if (shouldProcessCnvs(variant)) {
            processCnvs(variant);
        }
        vcfWriter.add(builder.make());
    }

    private boolean shouldInitializeRdCn(final VariantContext variant) {
        return revisedEventsFiltered.containsKey(variant.getID());
    }

    private boolean shouldProcessVariant(final VariantContext variant) {
        return revisedEventsAll.containsKey(variant.getID());
    }

    private boolean shouldProcessCnvs(final VariantContext variant) {
        final String svType = variant.getAttributeAsString(GATKSVVCFConstants.SVTYPE, "");
        final boolean isDelDup = svType.equals(GATKSVVCFConstants.SYMB_ALT_STRING_DUP) || svType.equals(GATKSVVCFConstants.SYMB_ALT_STRING_DEL);
        final boolean isLarge = variant.getEnd() - variant.getStart() >= 1000;
        return isDelDup && isLarge;
    }

    private void initializeRdCn(final VariantContext variant) {
        // Initialize data structures
        final String variantId = variant.getID();
        final Set<String> samples = revisedEventsFiltered.get(variantId);
        final Map<String, Integer> variantRdCn = new HashMap<>();

        // Initialize revisedRdCn value for each variant
        for (final String sampleName : samples) {
            final Genotype genotype = variant.getGenotype(sampleName);
            final int rdCn = (int) genotype.getExtendedAttribute(GATKSVVCFConstants.RD_CN);
            variantRdCn.put(sampleName, rdCn);
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
            final int rdCn = (int) genotype.getExtendedAttribute(GATKSVVCFConstants.RD_CN);
            if ((isDel && rdCn > 3) || (!isDel && (rdCn < 1 || rdCn > 4))) {
                multiCnvs.add(variant.getID());
                break;
            }
        }
    }

    private void processBedFile() {
        try {
            String line;
            final BufferedReader reader = new BufferedReader(new InputStreamReader(new GZIPInputStream(Files.newInputStream(bedFile.toPath()))));
            while ((line = reader.readLine()) != null) {
                final String[] fields = line.split("\t");
                if (fields.length < 12) continue;

                final String[] wider = Integer.parseInt(fields[2]) - Integer.parseInt(fields[1]) >= Integer.parseInt(fields[8]) - Integer.parseInt(fields[7])
                        ? Arrays.copyOfRange(fields, 0, 6)
                        : Arrays.copyOfRange(fields, 6, 12);
                final String[] narrower = Integer.parseInt(fields[2]) - Integer.parseInt(fields[1]) >= Integer.parseInt(fields[8]) - Integer.parseInt(fields[7])
                        ? Arrays.copyOfRange(fields, 6, 12)
                        : Arrays.copyOfRange(fields, 0, 6);
                if (wider[5].equals(GATKSVVCFConstants.BLANK_SAMPLES)) continue;

                double coverage = getCoverage(wider, narrower);
                if (coverage >= 0.5) {
                    Set<String> widerSamples = new HashSet<>(Arrays.asList(wider[5].split(",")));
                    Set<String> narrowerSamples = new HashSet<>(Arrays.asList(narrower[5].split(",")));
                    Set<String> uniqueSamples = new HashSet<>(widerSamples);
                    uniqueSamples.removeAll(narrowerSamples);

                    for (String sample : uniqueSamples) {
                        revisedEventsAll.computeIfAbsent(narrower[3], k -> new HashMap<>())
                                .put(sample, new ImmutablePair<>(wider[3], wider[4]));
                    }
                }
            }
            reader.close();

            for (final Map.Entry<String, Map<String, Pair<String, String>>> entry : revisedEventsAll.entrySet()) {
                for (final Map.Entry<String, Pair<String, String>> innerEntry : entry.getValue().entrySet()) {
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
        } catch (IOException e) {
            throw new RuntimeException("Error reading bed file", e);
        }
    }

    private double getCoverage(String[] wider, String[] narrower) {
        final int nStart = Integer.parseInt(narrower[1]);
        final int nStop = Integer.parseInt(narrower[2]);
        final int wStart = Integer.parseInt(wider[1]);
        final int wStop = Integer.parseInt(wider[2]);

        if (wStart <= nStop && nStart <= wStop) {
            final int intersectionSize = Math.min(nStop, wStop) - Math.max(nStart, wStart);
            return (double) intersectionSize / (nStop - nStart);
        }
        return 0.0;
    }
}
