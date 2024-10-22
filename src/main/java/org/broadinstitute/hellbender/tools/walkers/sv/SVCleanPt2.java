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
    private final Map<String, Map<String, Integer>> revisedCopyNumbers = new HashMap<>();
    private final Set<String> revisedComplete = new HashSet<>();

    private static final int MIN_VARIANT_SIZE = 5000;

    @Override
    protected int numberOfPasses() {
        return 2;
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
    }

    @Override
    protected void nthPassApply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext, int n) {
        // Skip if not expected SVTYPE or below SVLEN threshold
        if (!isDelDup(variant) || !isLargeVariant(variant, MIN_VARIANT_SIZE)) {
            return;
        }

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
    }

    @Override
    protected void afterNthPass(int n) {
        return;
    }

    private void firstPassApply(VariantContext variant) {
        // Flag sample as having an abnormal copy number if it passes certain conditions
        for (String sample : variant.getSampleNames()) {
            Genotype genotype = variant.getGenotype(sample);
            Integer rdCn = genotype.hasExtendedAttribute("RD_CN") ? Integer.parseInt(genotype.getExtendedAttribute("RD_CN").toString()) : null;
            if (!sampleWhitelist.contains(sample) || !genotype.isCalled() || rdCn == null || rdCn == 2) {
                continue;
            }

            String svType = variant.getAttributeAsString(GATKSVVCFConstants.SVTYPE, "");
            if ((svType.equals(GATKSVVCFConstants.SYMB_ALT_STRING_DEL) && rdCn < 2) || (svType.equals(GATKSVVCFConstants.SYMB_ALT_STRING_DUP) && rdCn > 2)) {
                abnormalRdCn.computeIfAbsent(variant.getID(), k -> new HashSet<>()).add(sample);
            }
        }

        // Add variant to overlap detector
        overlapDetector.addLhs(variant, variant);
    }

    private void secondPassApply(VariantContext variant) {
        // Check if copy number needs to be adjusted for samples within overlapping variants
        Set<VariantContext> overlappingVariants = overlapDetector.getOverlaps(variant);
        for (VariantContext otherVariant : overlappingVariants) {
            if (!variant.getID().equals(otherVariant.getID())) {
                adjustCopyNumber(variant, otherVariant);
            }
        }
    }

    private void adjustCopyNumber(VariantContext v1, VariantContext v2) {
        // Track metadata through data structures
        String variantId1 = v1.getID();
        String variantId2 = v2.getID();
        Map<String, Integer> variantRdCn1 = getRdCnForVariant(v1);
        Map<String, Integer> variantRdCn2 = getRdCnForVariant(v2);
        Map<String, Set<String>> variantSupport1 = getSupportForVariant(v1);
        Map<String, Set<String>> variantSupport2 = getSupportForVariant(v2);
        String svtype1 = v1.getAttributeAsString("SVTYPE", "");
        String svtype2 = v2.getAttributeAsString("SVTYPE", "");

        // Calculate overlap metadata
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
            Genotype genotype2 = v2.getGenotype(sample);

            // Condition 1: Smaller depth call is being driven by a larger call
            if (support1.contains(GATKSVVCFConstants.EV_VALUES.get(1)) && support1.size() > 1
                    && support2.equals(Collections.singleton(GATKSVVCFConstants.EV_VALUES.get(1)))
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
            else if (support1.equals(Collections.singleton(GATKSVVCFConstants.EV_VALUES.get(1)))
                    && support2.contains(GATKSVVCFConstants.EV_VALUES.get(1)) && support2.size() > 1
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
            else if (support1.equals(Collections.singleton(GATKSVVCFConstants.EV_VALUES.get(1)))
                    && support2.equals(Collections.singleton(GATKSVVCFConstants.EV_VALUES.get(1)))
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
            else if (support1.contains(GATKSVVCFConstants.EV_VALUES.get(1)) && overlap2 > 0.5 && !multiallelicCnvs.contains(variantId1)
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
        return svType.equals(GATKSVVCFConstants.SYMB_ALT_STRING_DEL) || svType.equals(GATKSVVCFConstants.SYMB_ALT_STRING_DUP);
    }

    private boolean isLargeVariant(VariantContext variant, int minSize) {
        int variantLength = Math.abs(variant.getAttributeAsInt(GATKSVVCFConstants.SVLEN, 0));
        return variantLength >= minSize;
    }

    private Map<String, Set<String>> getSupportForVariant(VariantContext variant) {
        Map<String, Set<String>> supportMap = new HashMap<>();
        for (String sample : variant.getSampleNames()) {
            Genotype genotype = variant.getGenotype(sample);
            String supportStr = genotype.hasExtendedAttribute(GATKSVVCFConstants.EV) ? genotype.getExtendedAttribute(GATKSVVCFConstants.EV).toString() : "";
            Set<String> supportSet = new HashSet<>();
            if (!supportStr.isEmpty()) {
                supportSet.addAll(Arrays.asList(supportStr.split(",")));
            }
            supportMap.put(sample, supportSet);
        }
        return supportMap;
    }

    private Map<String, Integer> getRdCnForVariant(VariantContext variant) {
        Map<String, Integer> rdCnMap = new HashMap<>();
        for (String sample : variant.getSampleNames()) {
            Genotype genotype = variant.getGenotype(sample);
            if (genotype.hasExtendedAttribute(GATKSVVCFConstants.RD_CN)) {
                rdCnMap.put(sample, Integer.parseInt(genotype.getExtendedAttribute(GATKSVVCFConstants.RD_CN).toString()));
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
}
