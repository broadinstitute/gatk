package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

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
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.Map;
import java.util.HashSet;
import java.util.HashMap;

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
 *     gatk SVCleanPt2 \
 *       -V input.vcf.gz \
 *       --sample-list samples.txt \
 * 	     --multi-cnv-list multi.cnvs.txt
 * 	     --output-prefix result
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
        summary = "Clean and format SV VCF",
        oneLineSummary = "Clean and format SV VCF",
        programGroup = StructuralVariantDiscoveryProgramGroup.class
)
@BetaFeature
@DocumentedFeature
public class SVCleanPt2 extends MultiplePassVariantWalker {
    public static final String SAMPLE_LIST_LONG_NAME = "sample-list";
    public static final String MULTI_CNV_LONG_NAME = "multi-cnv-list";
    public static final String OUTPUT_PREFIX_LONG_NAME = "output-prefix";

    @Argument(
            fullName = SAMPLE_LIST_LONG_NAME,
            doc = "Samples to include"
    )
    private GATKPath sampleListPath;

    @Argument(
            fullName = MULTI_CNV_LONG_NAME,
            doc = "List of multiallelic CNVs"
    )
    private GATKPath multiCnvPath;

    @Argument(
            fullName = OUTPUT_PREFIX_LONG_NAME,
            doc = "Prefix for output files"
    )
    private String outputPrefix;

    private BufferedWriter revisedCnWriter;

    private Set<String> sampleWhitelist;
    private Set<String> multiallelicCnvs;

    private Set<String> wasRevisedToNormal = new HashSet<>();
    private Map<String, Map<String, Integer>> revisedCopyNumbers = new HashMap<>();
    private final Map<String, Set<String>> variantToSamplesWithAbnormalCN = new HashMap<>();
    private final List<VariantContext> variantBuffer = new ArrayList<>();
    private final Map<String, Integer> variantLengths = new HashMap<>();

    @Override
    protected int numberOfPasses() {
        return 3;
    }

    @Override
    public void onTraversalStart() {
        try {
            sampleWhitelist = new HashSet<>(Files.readAllLines(sampleListPath.toPath()));
            multiallelicCnvs = new HashSet<>(Files.readAllLines(multiCnvPath.toPath()));
        } catch (IOException e) {
            throw new RuntimeException("Error reading input file", e);
        }
    }

    @Override
    public Object onTraversalSuccess() {
        try {
            for (Map.Entry<String, Map<String, Integer>> entry : revisedCopyNumbers.entrySet()) {
                String variantID = entry.getKey();
                for (Map.Entry<String, Integer> sampleEntry : entry.getValue().entrySet()) {
                    String sample = sampleEntry.getKey();
                    int rdCn = sampleEntry.getValue();
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
            default:
                throw new IllegalArgumentException("Invalid pass number: " + n);
        }
    }

    @Override
    protected void afterNthPass(int n) {
        if (n == 2) {
            try {
                revisedCnWriter = Files.newBufferedWriter(Paths.get(outputPrefix + ".txt"));
            } catch (IOException e) {
                throw new RuntimeException("Error creating output files", e);
            }
        }
    }

    private void firstPassApply(VariantContext variant) {
        // Skip variants not in DEL or DUP
        String svType = variant.getAttributeAsString(GATKSVVCFConstants.SVTYPE, "");
        if (!svType.equals("<DEL>") && !svType.equals("<DUP>")) {
            return;
        }

        // Process each sample
        for (String sample : variant.getSampleNames()) {
            if (!sampleWhitelist.contains(sample)) {
                continue;
            }
            Genotype genotype = variant.getGenotype(sample);
            if (!genotype.isCalled()) {
                continue;
            }
            Integer rdCn = genotype.hasExtendedAttribute("RD_CN") ?
                    Integer.parseInt(genotype.getExtendedAttribute("RD_CN").toString()) : null;
            if (rdCn == null || rdCn == 2) {
                continue;
            }
            if ((svType.equals("<DEL>") && rdCn < 2) || (svType.equals("<DUP>") && rdCn > 2)) {
                variantToSamplesWithAbnormalCN.computeIfAbsent(variant.getID(), k -> new HashSet<>()).add(sample);
            }
        }

        // Store variant length
        int variantLength = Math.abs(variant.getAttributeAsInt("SVLEN", 0));
        variantLengths.put(variant.getID(), variantLength);

        // Add to variant buffer for overlap detection in the next pass
        variantBuffer.add(variant);
    }

    private void secondPassApply(VariantContext variant) {
        String variantID = variant.getID();
        VariantContext currentVariant = variantBuffer.stream()
                .filter(vc -> vc.getID().equals(variantID))
                .findFirst()
                .orElse(null);
        if (currentVariant == null) {
            return;
        }

        // Find overlapping variants
        for (VariantContext otherVariant : variantBuffer) {
            if (variantID.equals(otherVariant.getID())) {
                continue;
            }
            if (variantsOverlap(currentVariant, otherVariant)) {
                // Apply the logic from the script to adjust RD_CN values
                adjustCopyNumbers(currentVariant, otherVariant);
            }
        }
    }

    private void thirdPassApply(VariantContext variant) {
        VariantContextBuilder builder = new VariantContextBuilder(variant);
        String variantID = variant.getID();
        Map<String, Integer> revisedRdCnForVariant = revisedCopyNumbers.getOrDefault(variantID, Collections.emptyMap());
        List<Genotype> newGenotypes = new ArrayList<>();

        // Build the set of alleles for the variant
        List<Allele> variantAlleles = new ArrayList<>(variant.getAlleles());
        boolean variantAllelesModified = false;

        for (Genotype genotype : variant.getGenotypes()) {
            String sample = genotype.getSampleName();
            Integer revisedRdCn = revisedRdCnForVariant.get(sample);
            if (revisedRdCn != null) {
                // Create a new genotype with the revised RD_CN
                GenotypeBuilder gb = new GenotypeBuilder(genotype);
                gb.attribute("RD_CN", revisedRdCn);

                // Adjust GT and alleles if necessary
                if (revisedRdCn == 2) {
                    // Homozygous reference
                    gb.alleles(Arrays.asList(variant.getReference(), variant.getReference()));
                    gb.GQ(99); // Example GQ value for homozygous reference
                } else {
                    // Heterozygous or other genotype
                    Allele altAllele;
                    if (variant.getAlternateAlleles().isEmpty()) {
                        // Need to create ALT allele
                        String svType = variant.getAttributeAsString("SVTYPE", null);
                        if (svType == null) {
                            throw new IllegalArgumentException("SVTYPE is missing for variant " + variantID);
                        }
                        altAllele = Allele.create("<" + svType + ">", false);
                        variantAlleles.add(altAllele);
                        variantAllelesModified = true;
                    } else {
                        altAllele = variant.getAlternateAllele(0);
                    }
                    gb.alleles(Arrays.asList(variant.getReference(), altAllele));
                }

                newGenotypes.add(gb.make());
            } else {
                newGenotypes.add(genotype);
            }
        }

        // Update the variant's alleles if modified
        if (variantAllelesModified) {
            builder.alleles(variantAlleles);
        }

        builder.genotypes(newGenotypes);
        VariantContext updatedVariant = builder.make();
        identifyMultiallelicCnvs(updatedVariant);
    }

    private boolean variantsOverlap(VariantContext v1, VariantContext v2) {
        return v1.getContig().equals(v2.getContig()) &&
                v1.getStart() <= v2.getEnd() &&
                v2.getStart() <= v1.getEnd();
    }

    private void adjustCopyNumbers(VariantContext v1, VariantContext v2) {
        // Determine larger and smaller variants
        VariantContext largerVariant = (variantLengths.get(v1.getID()) >= variantLengths.get(v2.getID())) ? v1 : v2;
        VariantContext smallerVariant = (largerVariant == v1) ? v2 : v1;

        // Calculate overlap
        int overlapStart = Math.max(largerVariant.getStart(), smallerVariant.getStart());
        int overlapEnd = Math.min(largerVariant.getEnd(), smallerVariant.getEnd());
        int overlapLength = overlapEnd - overlapStart + 1;
        double overlapPercentageSmaller = (double) overlapLength / (smallerVariant.getEnd() - smallerVariant.getStart() + 1);
        double overlapPercentageLarger = (double) overlapLength / (largerVariant.getEnd() - largerVariant.getStart() + 1);

        // Apply logic based on support type and other conditions
        // (Implementation of specific conditions from the script)
        // For brevity, let's assume we have a method that applies these conditions
        applyAdjustmentLogic(largerVariant, smallerVariant, overlapPercentageSmaller, overlapPercentageLarger);
    }

    private void applyAdjustmentLogic(VariantContext largerVariant, VariantContext smallerVariant,
                                      double overlapSmaller, double overlapLarger) {

        String smallerVariantID = smallerVariant.getID();
        String largerVariantID = largerVariant.getID();
        Map<String, Integer> smallerVariantRdCn = getRdCnForVariant(smallerVariant);
        Map<String, Integer> largerVariantRdCn = getRdCnForVariant(largerVariant);
        Map<String, String> smallerVariantSupport = getSupportForVariant(smallerVariant);
        Map<String, String> largerVariantSupport = getSupportForVariant(largerVariant);
        Map<String, String> smallerVariantGT = getGTForVariant(smallerVariant);
        Map<String, String> largerVariantGT = getGTForVariant(largerVariant);
        String svtype1 = smallerVariant.getAttributeAsString("SVTYPE", "");
        String svtype2 = largerVariant.getAttributeAsString("SVTYPE", "");

        // Lengths of the variants
        int length1 = smallerVariant.getEnd() - smallerVariant.getStart();
        int length2 = largerVariant.getEnd() - largerVariant.getStart();

        // Iterate over samples present in both variants
        Set<String> samples = new HashSet<>(smallerVariant.getSampleNames());
        samples.retainAll(largerVariant.getSampleNames());

        for (String sample : samples) {
            String id1 = smallerVariantID + "@" + sample;
            String id2 = largerVariantID + "@" + sample;

            // Check if id1 has already been revised to normal
            if (wasRevisedToNormal.contains(id1)) {
                continue;
            }

            // Retrieve or update RD_CN values if they have been revised already
            Integer RD_CN1 = revisedCopyNumbers.getOrDefault(smallerVariantID, Collections.emptyMap()).getOrDefault(sample, smallerVariantRdCn.get(sample));
            Integer RD_CN2 = revisedCopyNumbers.getOrDefault(largerVariantID, Collections.emptyMap()).getOrDefault(sample, largerVariantRdCn.get(sample));

            String support1 = smallerVariantSupport.get(sample);
            String support2 = largerVariantSupport.get(sample);
            String GT1 = smallerVariantGT.get(sample);
            String GT2 = largerVariantGT.get(sample);

            // Ensure RD_CN values are not null
            if (RD_CN1 == null || RD_CN2 == null) {
                continue;
            }

            // Calculate overlaps
            boolean smallOverlap50 = overlapSmaller > 0.5;
            boolean largeOverlap50 = overlapLarger > 0.5;

            // Apply the conditions from the shell script

            // Condition 1: Smaller depth call is being driven by larger
            if (support1.contains("RD") && !support1.equals("RD") && support2.equals("RD") &&
                    smallOverlap50 && !multiallelicCnvs.contains(smallerVariantID)) {

                if (RD_CN1 == 0) {
                    makeRevision(id2, RD_CN2 + 2);
                } else if (RD_CN1 == 1) {
                    makeRevision(id2, RD_CN2 + RD_CN1);
                } else if (RD_CN1 > 1) {
                    int newCN = RD_CN2 - RD_CN1 + 2;
                    newCN = Math.max(newCN, 0);
                    makeRevision(id2, newCN);
                }
            }

            // Condition 2: Smaller CNV driving larger CNV genotype
            else if (support1.equals("RD") && support2.contains("RD") && !support2.equals("RD") &&
                    smallOverlap50 && !multiallelicCnvs.contains(largerVariantID) &&
                    !GT2.equals("0/0") && largeOverlap50) {

                if (RD_CN2 == 0) {
                    makeRevision(id1, RD_CN1 + 2);
                } else if (RD_CN2 == 1) {
                    makeRevision(id1, RD_CN1 + RD_CN2);
                } else if (RD_CN2 > 1) {
                    int newCN = RD_CN1 - RD_CN2 + 2;
                    newCN = Math.max(newCN, 0);
                    makeRevision(id1, newCN);
                }
            }

            // Condition 3: Depth-only calls where smaller call is being driven by larger
            else if (support1.equals("RD") && support2.equals("RD") && smallOverlap50 &&
                    svtype1.equals(svtype2) && !multiallelicCnvs.contains(smallerVariantID)) {

                if (RD_CN1 == 0 && !RD_CN1.equals(RD_CN2)) {
                    makeRevision(id2, RD_CN2 + 2);
                } else if (RD_CN1 == 1 && RD_CN1 > RD_CN2) {
                    makeRevision(id2, 1);
                } else if (RD_CN1 > 1 && RD_CN1 < RD_CN2) {
                    int newCN = RD_CN2 - RD_CN1 + 2;
                    newCN = Math.max(newCN, 0);
                    makeRevision(id2, newCN);
                } else {
                    makeRevision(id2, 2);
                }
            }

            // Condition 4: Any other time a larger call is driving a smaller call
            else if (support1.contains("RD") && smallOverlap50 && length2 > 5000 &&
                    !multiallelicCnvs.contains(smallerVariantID)) {

                if (RD_CN1 == 0) {
                    makeRevision(id2, RD_CN2 + 2);
                } else if (RD_CN1 == 1) {
                    makeRevision(id2, RD_CN2 + RD_CN1);
                } else if (RD_CN1 > 1) {
                    int newCN = RD_CN2 - RD_CN1 + 2;
                    newCN = Math.max(newCN, 0);
                    makeRevision(id2, newCN);
                }
            }
        }
    }

    private void makeRevision(String id, int val) {
        // id is in the format variantID@sample
        String[] tokens = id.split("@");
        String variantID = tokens[0];
        String sample = tokens[1];
        revisedCopyNumbers.computeIfAbsent(variantID, k -> new HashMap<>()).put(sample, val);
        if (val == 2) {
            wasRevisedToNormal.add(id);
        }
    }

    private Map<String, String> getSupportForVariant(VariantContext variant) {
        Map<String, String> supportMap = new HashMap<>();
        for (String sample : variant.getSampleNames()) {
            Genotype genotype = variant.getGenotype(sample);
            String support = genotype.hasExtendedAttribute("EV") ?
                    genotype.getExtendedAttribute("EV").toString() : "";
            supportMap.put(sample, support);
        }
        return supportMap;
    }

    private Map<String, String> getGTForVariant(VariantContext variant) {
        Map<String, String> gtMap = new HashMap<>();
        for (String sample : variant.getSampleNames()) {
            Genotype genotype = variant.getGenotype(sample);
            String gt = genotype.isCalled() ? genotype.getGenotypeString() : "./.";
            gtMap.put(sample, gt);
        }
        return gtMap;
    }

    private Map<String, Integer> getRdCnForVariant(VariantContext variant) {
        Map<String, Integer> rdCnMap = new HashMap<>();
        for (String sample : variant.getSampleNames()) {
            Genotype genotype = variant.getGenotype(sample);
            if (genotype.hasExtendedAttribute("RD_CN")) {
                rdCnMap.put(sample, Integer.parseInt(genotype.getExtendedAttribute("RD_CN").toString()));
            }
        }
        return rdCnMap;
    }

    private void identifyMultiallelicCnvs(VariantContext variant) {
        String svType = variant.getAttributeAsString(GATKSVVCFConstants.SVTYPE, "");
        boolean isDel = svType.equals("<DEL>");
        boolean isDup = svType.equals("<DUP>");
        int variantLength = variantLengths.getOrDefault(variant.getID(), 0);
        if ((isDel || isDup) && variantLength >= 5000) {
            for (Genotype genotype : variant.getGenotypes()) {
                Integer rdCn = genotype.hasExtendedAttribute("RD_CN") ?
                        Integer.parseInt(genotype.getExtendedAttribute("RD_CN").toString()) : null;
                if (rdCn != null) {
                    if (isDel && rdCn > 3) {
                        // Multiallelic deletion
                        multiallelicCnvs.add(variant.getID());
                        break;
                    } else if (isDup && (rdCn < 1 || rdCn > 4)) {
                        // Multiallelic duplication
                        multiallelicCnvs.add(variant.getID());
                        break;
                    }
                }
            }
        }
    }
}
