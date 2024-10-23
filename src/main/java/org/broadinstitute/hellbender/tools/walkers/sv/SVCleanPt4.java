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
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
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
import java.util.zip.GZIPOutputStream;

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
 *     gatk SVCleanPt4 \
 *       -V input.vcf.gz \
 *       --revised-cn-list revised.txt \
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
        summary = "SClean and format SV VCF",
        oneLineSummary = "Clean and format SV VCF",
        programGroup = StructuralVariantDiscoveryProgramGroup.class
)
@BetaFeature
@DocumentedFeature
public class SVCleanPt4 extends VariantWalker {
    public static final String REVISED_CN_LIST_LONG_NAME = "revised-cn-list";
    public static final String OUTPUT_PREFIX_LONG_NAME = "output-prefix";

    @Argument(
            fullName = REVISED_CN_LIST_LONG_NAME,
            doc = "File with variant IDs, sample IDs, and RD_CN values"
    )
    private GATKPath cnReviseList;

    @Argument(
            fullName = OUTPUT_PREFIX_LONG_NAME,
            doc = "Prefix for output files"
    )
    private String outputPrefix;

    private VariantContextWriter vcfWriter;
    private BufferedWriter multiGenoWriter;

    private Map<String, Map<String, Integer>> revisedCopyNumbers;
    private final Set<String> multiGenoIds = new HashSet<>();

    private double recordStart;
    private double recordEnd;
    private int maxVF;
    private long recordIdx;

    @Override
    public void onTraversalStart() {
        // Read revised copy numbers
        revisedCopyNumbers = readRevisedEvents(cnReviseList);

        // Parse batchNum and totalBatch from file name
        String cnReviseListFileName = cnReviseList.toPath().getFileName().toString();
        String[] regenoFileNameTokens = cnReviseListFileName.split("\\.");
        String[] batchTokens = regenoFileNameTokens[1].split("_");
        int batchNum = Math.max(Integer.parseInt(batchTokens[0]), 1);
        int totalBatch = Math.max(Integer.parseInt(batchTokens[1]), 1);

        // Get variant count
        long totalNumVariants = 0;
        String inputVcfPath = getDrivingVariantsFeatureInput().getFeaturePath();
        try (FeatureDataSource<VariantContext> dataSource = new FeatureDataSource<>(inputVcfPath)) {
            for (VariantContext vc : dataSource) {
                totalNumVariants++;
            }
        }

        // Initialize metadata variables
        double segments = totalNumVariants / (double) totalBatch;
        recordStart = (batchNum - 1) * segments;
        recordEnd = batchNum * segments;
        maxVF = Math.max((int) (getHeaderForVariants().getGenotypeSamples().size() * 0.01), 2);
        recordIdx = 0;

        // Create output writers
        try {
            vcfWriter = createVCFWriter(Paths.get(outputPrefix + ".revised_vcf_lines.txt"));
            vcfWriter.writeHeader(getHeaderForVariants());

            multiGenoWriter = Files.newBufferedWriter(Paths.get(outputPrefix + ".multi_geno_ids.txt"));
        } catch (IOException e) {
            throw new RuntimeException("Error creating output file", e);
        }
    }

    @Override
    public Object onTraversalSuccess() {
        try {
            List<String> variantIDs = new ArrayList<>(multiGenoIds);
            Collections.sort(variantIDs);
            for (String variantID : variantIDs) {
                multiGenoWriter.write(variantID);
                multiGenoWriter.newLine();
            }
            return null;
        } catch (IOException e) {
            throw new RuntimeException("Error writing to output file ", e);
        }
    }

    public void closeTool() {
        try {
            if (vcfWriter != null) {
                vcfWriter.close();
            }
            if (multiGenoWriter != null) {
                multiGenoWriter.close();
            }
        } catch (IOException e) {
            throw new RuntimeException("Error closing output file ", e);
        }
    }

    @Override
    public void apply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        // Initialize data structures
        recordIdx++;
        VariantContextBuilder variantBuilder = new VariantContextBuilder(variant);
        List<Genotype> genotypes = variant.getGenotypes();

        // Modify genotypes if variant appears in revise list
        if (revisedCopyNumbers.containsKey(variant.getID())) {
            Map<String, Integer> sampleCnMap = revisedCopyNumbers.get(variant.getID());
            List<Genotype> newGenotypes = new ArrayList<>();
            for (Genotype genotype : genotypes) {
                String sampleName = genotype.getSampleName();
                if (sampleCnMap.containsKey(sampleName)) {
                    GenotypeBuilder gb = new GenotypeBuilder(genotype);
                    gb.alleles(Arrays.asList(variant.getReference(), variant.getAlternateAllele(0)));
                    gb.attribute(GATKSVVCFConstants.RD_CN, sampleCnMap.get(sampleName));
                    newGenotypes.add(gb.make());
                } else {
                    newGenotypes.add(genotype);
                }
            }
            variantBuilder.genotypes(newGenotypes);
            vcfWriter.add(variantBuilder.make());
        }

        // Identify multiple genotypes if within recordStart and recordEnd
        if (recordIdx >= recordStart && recordIdx < recordEnd) {
            int numGtOver2 = 0;
            for (Genotype genotype : genotypes) {
                Integer peGt = genotype.hasExtendedAttribute(GATKSVVCFConstants.PE_GT) ?
                        Integer.parseInt(genotype.getExtendedAttribute(GATKSVVCFConstants.PE_GT).toString()) : null;
                Integer srGt = genotype.hasExtendedAttribute(GATKSVVCFConstants.SR_GT) ?
                        Integer.parseInt(genotype.getExtendedAttribute(GATKSVVCFConstants.SR_GT).toString()) : null;
                Integer gt = null;
                if (peGt == null) {
                    continue;
                } else if (srGt == null) {
                    gt = peGt;
                } else if (peGt > 0 && srGt == 0) {
                    gt = peGt;
                } else if (peGt == 0) {
                    gt = srGt;
                } else {
                    Integer peGq = genotype.hasExtendedAttribute(GATKSVVCFConstants.PE_GQ) ?
                            Integer.parseInt(genotype.getExtendedAttribute(GATKSVVCFConstants.PE_GQ).toString()) : null;
                    Integer srGq = genotype.hasExtendedAttribute(GATKSVVCFConstants.SR_GQ) ?
                            Integer.parseInt(genotype.getExtendedAttribute(GATKSVVCFConstants.SR_GQ).toString()) : null;
                    if (peGq != null && srGq != null && peGq >= srGq) {
                        gt = peGt;
                    } else {
                        gt = srGt;
                    }
                }
                if (gt > 2) {
                    numGtOver2++;
                }
            }
            if (numGtOver2 > maxVF) {
                multiGenoIds.add(variant.getID());
            }
        }
    }

    private Integer getIntegerAttribute(Genotype genotype, String attributeName) {
        if (genotype.hasExtendedAttribute(attributeName)) {
            Object attr = genotype.getExtendedAttribute(attributeName);
            if (attr instanceof Integer) {
                return (Integer) attr;
            } else if (attr instanceof String) {
                try {
                    return Integer.parseInt((String) attr);
                } catch (NumberFormatException e) {
                    return null;
                }
            }
        }
        return null;
    }

    private Map<String, Map<String, Integer>> readRevisedEvents(final GATKPath filePath) {
        try (BufferedReader reader = new BufferedReader(new FileReader(filePath.toPath().toFile()))) {
            final Map<String, Map<String, Integer>> result = new HashMap<>();
            String line;
            while ((line = reader.readLine()) != null) {
                String[] fields = line.split("\t");
                if (fields.length < 3) continue;

                String variantId = fields[0];
                String sampleId = fields[1];
                int rdCn = Integer.parseInt(fields[2]);

                result.computeIfAbsent(variantId, k -> new HashMap<>()).put(sampleId, rdCn);
            }
            return result;
        } catch (IOException e) {
            throw new RuntimeException("Error reading input file", e);
        }
    }
}
