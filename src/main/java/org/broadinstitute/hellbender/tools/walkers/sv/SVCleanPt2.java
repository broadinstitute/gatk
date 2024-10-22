package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
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

    private final Map<String, Set<String>> abnormalRdCn = new HashMap<>();
    private OverlapDetector<VariantContext> overlapDetector = new OverlapDetector<>(0, 0);
    private final Map<String, Map<String, Integer>> revisedCopyNumbers = new HashMap<>(); // STATUS: To Be Verified
    private final Set<String> revisedComplete = new HashSet<>(); // STATUS:  To Be Verified

    private static final int MIN_VARIANT_SIZE = 5000;

    @Override
    protected int numberOfPasses() {
        return 3;
    }

    @Override
    public void onTraversalStart() {
        try {
            revisedCnWriter = Files.newBufferedWriter(Paths.get(outputPrefix + ".txt"));

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
        return;
    }

    private void firstPassApply(VariantContext variant) {
        // Skip if not expected SVTYPE or SVLEN
        if (!isDelDup(variant) || !isLargeVariant(variant, MIN_VARIANT_SIZE)) {
            return;
        }

        // Flag sample as having abnormal copy number if it passes various conditions
        for (String sample : variant.getSampleNames()) {
            Genotype genotype = variant.getGenotype(sample);
            Integer rdCn = genotype.hasExtendedAttribute("RD_CN") ? Integer.parseInt(genotype.getExtendedAttribute("RD_CN").toString()) : null;
            if (!sampleWhitelist.contains(sample) || !genotype.isCalled() || rdCn == null || rdCn == 2) {
                continue;
            }

            String svType = variant.getAttributeAsString(GATKSVVCFConstants.SVTYPE, "");
            if ((svType.equals("DEL") && rdCn < 2) || (svType.equals("DUP") && rdCn > 2)) {
                abnormalRdCn.computeIfAbsent(variant.getID(), k -> new HashSet<>()).add(sample);
            }
        }

        // Add variant to overlap detector
        overlapDetector.addLhs(variant, variant);
    }

    private void secondPassApply(VariantContext variant) {
        // Skip if not expected SVTYPE or SVLEN
        if (!isDelDup(variant) || !isLargeVariant(variant, MIN_VARIANT_SIZE)) {
            return;
        }

        // Adjust copy numbers for overlapping variants
        Set<VariantContext> overlappingVariants = overlapDetector.getOverlaps(variant);
        for (VariantContext otherVariant : overlappingVariants) {
            if (!variant.getID().equals(otherVariant.getID())) {
                adjustCopyNumbers(variant, otherVariant);
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
                    gb.alleles(Arrays.asList(variant.getReference(), variant.getReference()));
                    gb.GQ(99);
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

    private void adjustCopyNumbers(VariantContext v1, VariantContext v2) {
        // Define data structures to store metadata
        String variantId1 = v1.getID();
        String variantId2 = v2.getID();
        Map<String, Integer> variantRdCn1 = getRdCnForVariant(v1);
        Map<String, Integer> variantRdCn2 = getRdCnForVariant(v2);
        Map<String, Set<String>> variantSupport1 = getSupportForVariant(v1);
        Map<String, Set<String>> variantSupport2 = getSupportForVariant(v2);
        Map<String, Genotype> variantGt1 = getGTForVariant(v1);
        Map<String, Genotype> variantGt2 = getGTForVariant(v2);
        String svtype1 = v1.getAttributeAsString("SVTYPE", "");
        String svtype2 = v2.getAttributeAsString("SVTYPE", "");

        // Calculate overlap
        int length1 = v1.getEnd() - v1.getStart();
        int length2 = v2.getEnd() - v2.getStart();
        int lengthOverlap = Math.min(v2.getEnd(), v1.getEnd()) - Math.max(v1.getStart(), v2.getStart());
        double overlap1 = (double) lengthOverlap / (double) length1;
        double overlap2 = (double) lengthOverlap / (double) length2;

        // Get samples with abnormal CN across both variants
        Set<String> samples = new HashSet<>(abnormalRdCn.getOrDefault(variantId1, Collections.emptySet()));
        samples.retainAll(abnormalRdCn.getOrDefault(variantId2, Collections.emptySet()));

        // Iterate through samples to test against conditions
        for (String sample : samples) {
            // Validate baseline filters
            String id1 = variantId1 + "@" + sample;
            String id2 = variantId2 + "@" + sample;
            Integer rdCn1 = revisedCopyNumbers.getOrDefault(variantId1, Collections.emptyMap()).getOrDefault(sample, variantRdCn1.get(sample));
            Integer rdCn2 = revisedCopyNumbers.getOrDefault(variantId2, Collections.emptyMap()).getOrDefault(sample, variantRdCn2.get(sample));
            if (revisedComplete.contains(id1) || rdCn1 == null || rdCn2 == null) {
                continue;
            }

            // Initialize fields for evaluation
            Set<String> support1 = variantSupport1.get(sample);
            Set<String> support2 = variantSupport2.get(sample);
            Genotype genotype1 = variantGt1.get(sample);
            Genotype genotype2 = variantGt2.get(sample);

            // Condition 1: Smaller depth call is being driven by a larger call
            if (support1.contains("RD") && support1.size() > 1 && support2.equals(Collections.singleton("RD"))
                    && overlap2 > 0.5 && !multiallelicCnvs.contains(variantId1)) {
                if (rdCn1 == 0) {
                    makeRevision(id2, rdCn2 + 2);
                } else if (rdCn1 == 1) {
                    makeRevision(id2, rdCn2 + rdCn1);
                } else if (rdCn1 > 1) {
                    int newCN = rdCn2 - rdCn1 + 2;
                    newCN = Math.max(newCN, 0);
                    makeRevision(id2, newCN);
                }
            }

            // Condition 2: Smaller CNV is driven by a larger CNV genotype
            else if (support1.equals(Collections.singleton("RD")) && support2.contains("RD") && support2.size() > 1
                    && overlap1 > 0.5 && overlap2 > 0.5 && !multiallelicCnvs.contains(variantId2) && !genotype2.isHomRef()) {
                if (rdCn2 == 0) {
                    makeRevision(id1, rdCn1 + 2);
                } else if (rdCn2 == 1) {
                    makeRevision(id1, rdCn1 + rdCn2);
                } else if (rdCn2 > 1) {
                    int newCN = rdCn1 - rdCn2 + 2;
                    newCN = Math.max(newCN, 0);
                    makeRevision(id1, newCN);
                }
            }

            // Condition 3: Depth-only calls where smaller call is driven by a larger call
            else if (support1.equals(Collections.singleton("RD")) && support2.equals(Collections.singleton("RD"))
                    && overlap2 > 0.5 && !multiallelicCnvs.contains(variantId1) && svtype1.equals(svtype2)) {
                if (rdCn1 == 0 && !rdCn1.equals(rdCn2)) {
                    makeRevision(id2, rdCn2 + 2);
                } else if (rdCn1 == 1 && rdCn1 > rdCn2) {
                    makeRevision(id2, 1);
                } else if (rdCn1 > 1 && rdCn1 < rdCn2) {
                    int newCN = rdCn2 - rdCn1 + 2;
                    newCN = Math.max(newCN, 0);
                    makeRevision(id2, newCN);
                } else {
                    makeRevision(id2, 2);
                }
            }

            // Condition 4: Any other time a larger call drives a smaller call
            else if (support1.contains("RD") && overlap2 > 0.5 && !multiallelicCnvs.contains(variantId1)
                    && length2 > MIN_VARIANT_SIZE) {
                if (rdCn1 == 0) {
                    makeRevision(id2, rdCn2 + 2);
                } else if (rdCn1 == 1) {
                    makeRevision(id2, rdCn2 + rdCn1);
                } else if (rdCn1 > 1) {
                    int newCN = rdCn2 - rdCn1 + 2;
                    newCN = Math.max(newCN, 0);
                    makeRevision(id2, newCN);
                }
            }
        }
    }

    private boolean isDelDup(VariantContext variant) {
        String svType = variant.getAttributeAsString(GATKSVVCFConstants.SVTYPE, "");
        return svType.equals("DEL") || svType.equals("DUP");
    }

    private boolean isLargeVariant(VariantContext variant, int minSize) {
        int variantLength = Math.abs(variant.getAttributeAsInt("SVLEN", 0));
        return variantLength >= minSize;
    }

    private Map<String, Set<String>> getSupportForVariant(VariantContext variant) {
        Map<String, Set<String>> supportMap = new HashMap<>();
        for (String sample : variant.getSampleNames()) {
            Genotype genotype = variant.getGenotype(sample);
            String supportStr = genotype.hasExtendedAttribute("EV") ? genotype.getExtendedAttribute("EV").toString() : "";
            Set<String> supportSet = new HashSet<>();
            if (!supportStr.isEmpty()) {
                supportSet.addAll(Arrays.asList(supportStr.split(",")));
            }
            supportMap.put(sample, supportSet);
        }
        return supportMap;
    }

    private Map<String, Genotype> getGTForVariant(VariantContext variant) {
        Map<String, Genotype> gtMap = new HashMap<>();
        for (String sample : variant.getSampleNames()) {
            Genotype genotype = variant.getGenotype(sample);
            gtMap.put(sample, genotype);
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

    private void makeRevision(String id, int val) {
        String[] tokens = id.split("@");
        String variantId = tokens[0];
        String sample = tokens[1];
        revisedCopyNumbers.computeIfAbsent(variantId, k -> new HashMap<>()).put(sample, val);
        if (val == 2) {
            revisedComplete.add(id);
        }
    }

    private void identifyMultiallelicCnvs(VariantContext variant) {
        if (isDelDup(variant) &&  isLargeVariant(variant, MIN_VARIANT_SIZE)) {
            for (Genotype genotype : variant.getGenotypes()) {
                Integer rdCn = genotype.hasExtendedAttribute("RD_CN") ? Integer.parseInt(genotype.getExtendedAttribute("RD_CN").toString()) : null;
                if (rdCn != null) {
                    String svType = variant.getAttributeAsString(GATKSVVCFConstants.SVTYPE, "");
                    if (svType.equals("DEL") && rdCn > 3) {
                        multiallelicCnvs.add(variant.getID());
                        break;
                    } else if (svType.equals("DUP") && (rdCn < 1 || rdCn > 4)) {
                        multiallelicCnvs.add(variant.getID());
                        break;
                    }
                }
            }
        }
    }
}
