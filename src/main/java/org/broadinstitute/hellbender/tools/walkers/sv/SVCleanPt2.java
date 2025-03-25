package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;

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
 *     TODO
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
public class SVCleanPt2 extends VariantWalker {
    public static final String SAMPLE_LIST_LONG_NAME = "sample-list";
    public static final String OUTPUT_REVISED_LIST_LONG_NAME = "output-revised-list";

    @Argument(
            fullName = SAMPLE_LIST_LONG_NAME,
            doc = "File with samples to include"
    )
    private GATKPath sampleListPath;

    @Argument(
            fullName = OUTPUT_REVISED_LIST_LONG_NAME,
            doc = "Prefix for output files"
    )
    private GATKPath outputRevisedList;

    private BufferedWriter revisedCnWriter;

    private Set<String> sampleWhitelist;

    private final Map<String, Set<String>> abnormalRdCn = new HashMap<>();
    private final List<VariantContext> overlappingVariantsBuffer = new ArrayList<>();
    private final Map<String, Map<String, Integer>> revisedCopyNumbers = new HashMap<>();
    private final Set<String> revisedComplete = new HashSet<>();

    private static final int MIN_VARIANT_SIZE = 5000;

    @Override
    public void onTraversalStart() {
        try {
            revisedCnWriter = Files.newBufferedWriter(Paths.get(outputRevisedList.toString()));

            sampleWhitelist = new HashSet<>(Files.readAllLines(sampleListPath.toPath()));
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
            throw new RuntimeException("Error writing output file", e);
        }
    }

    @Override
    public void apply(final VariantContext variant, final ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        // Skip if not expected SVTYPE or below SVLEN threshold
        if (!isDelDup(variant) || !isLarge(variant, MIN_VARIANT_SIZE)) {
            return;
        }

        // Flag sample as having an abnormal copy number if it passes certain conditions
        for (String sample : variant.getSampleNames()) {
            Genotype genotype = variant.getGenotype(sample);
            int rdCn = Integer.parseInt(genotype.getExtendedAttribute("RD_CN").toString());
            if (!sampleWhitelist.contains(sample) || !genotype.isCalled() || rdCn == 2) {
                continue;
            }

            String svType = variant.getAttributeAsString(GATKSVVCFConstants.SVTYPE, "");
            if ((svType.equals(GATKSVVCFConstants.SYMB_ALT_STRING_DEL) && rdCn < 2) || (svType.equals(GATKSVVCFConstants.SYMB_ALT_STRING_DUP) && rdCn > 2)) {
                abnormalRdCn.computeIfAbsent(variant.getID(), k -> new HashSet<>()).add(sample);
            }
        }

        // Process overlaps with variants in the buffer
        overlappingVariantsBuffer.removeIf(vc -> !vc.getContig().equals(variant.getContig()) || vc.getEnd() < variant.getStart());
        for (VariantContext bufferedVariant : overlappingVariantsBuffer) {
            if (overlaps(variant, bufferedVariant)) {
                adjustCopyNumber(bufferedVariant, variant);
                adjustCopyNumber(variant, bufferedVariant);
            }
        }
        overlappingVariantsBuffer.add(variant);
    }

    private void adjustCopyNumber(final VariantContext v1, final VariantContext v2) {
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
            int rdCn1 = revisedCopyNumbers.getOrDefault(variantId1, Collections.emptyMap()).getOrDefault(sample, variantRdCn1.get(sample));
            int rdCn2 = revisedCopyNumbers.getOrDefault(variantId2, Collections.emptyMap()).getOrDefault(sample, variantRdCn2.get(sample));
            if (revisedComplete.contains(id1)) {
                continue;
            }

            // Initialize fields for evaluation
            Set<String> support1 = variantSupport1.get(sample);
            Set<String> support2 = variantSupport2.get(sample);
            Genotype genotype2 = v2.getGenotype(sample);

            // Condition 1: Smaller depth call is being driven by a larger call
            if (support1.contains(GATKSVVCFConstants.EV_VALUES.get(1)) && support1.size() > 1
                    && support2.equals(Collections.singleton(GATKSVVCFConstants.EV_VALUES.get(1)))
                    && overlap2 > 0.5 && !v1.hasAttribute(GATKSVVCFConstants.MULTI_CNV)) {
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
                    && overlap1 > 0.5 && overlap2 > 0.5 && !v2.hasAttribute(GATKSVVCFConstants.MULTI_CNV)
                    && !genotype2.isHomRef()) {
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
                    && overlap2 > 0.5 && !v1.hasAttribute(GATKSVVCFConstants.MULTI_CNV) && svtype1.equals(svtype2)) {
                if (rdCn1 == 0 && rdCn1 != rdCn2) {
                    makeRevision(id2, rdCn2 + 2);
                } else if (rdCn1 == 1 && rdCn1 > rdCn2) {
                    makeRevision(id2, 1);
                } else if (rdCn1 > 1 && rdCn1 < rdCn2) {
                    makeRevision(id2, Math.max(rdCn2 - rdCn1 + 2, 0));
                } else {
                    makeRevision(id2, 2);
                }
            }

            // Condition 4: Any other time a larger call drives a smaller call
            else if (support1.contains(GATKSVVCFConstants.EV_VALUES.get(1))
                    && overlap2 > 0.5 && !v1.hasAttribute(GATKSVVCFConstants.MULTI_CNV) && length2 > MIN_VARIANT_SIZE) {
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

    private boolean overlaps(final VariantContext v1, final VariantContext v2) {
        return v1.getContig().equals(v2.getContig()) && v1.getStart() <= v2.getEnd() && v2.getStart() <= v1.getEnd();
    }

    private boolean isDelDup(final VariantContext variant) {
        String svType = variant.getAttributeAsString(GATKSVVCFConstants.SVTYPE, "");
        return svType.equals(GATKSVVCFConstants.SYMB_ALT_STRING_DEL) || svType.equals(GATKSVVCFConstants.SYMB_ALT_STRING_DUP);
    }

    private boolean isLarge(final VariantContext variant, final int minSize) {
        int variantLength = Math.abs(variant.getAttributeAsInt(GATKSVVCFConstants.SVLEN, 0));
        return variantLength >= minSize;
    }

    private Map<String, Set<String>> getSupportForVariant(final VariantContext variant) {
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

    private Map<String, Integer> getRdCnForVariant(final VariantContext variant) {
        Map<String, Integer> rdCnMap = new HashMap<>();
        for (String sample : variant.getSampleNames()) {
            Genotype genotype = variant.getGenotype(sample);
            if (genotype.hasExtendedAttribute(GATKSVVCFConstants.RD_CN)) {
                rdCnMap.put(sample, Integer.parseInt(genotype.getExtendedAttribute(GATKSVVCFConstants.RD_CN).toString()));
            }
        }
        return rdCnMap;
    }

    private void makeRevision(final String id, final int val) {
        String[] tokens = id.split("@");
        String variantId = tokens[0];
        String sample = tokens[1];
        revisedCopyNumbers.computeIfAbsent(variantId, k -> new HashMap<>()).put(sample, val);
        if (val == 2) {
            revisedComplete.add(id);
        }
    }
}