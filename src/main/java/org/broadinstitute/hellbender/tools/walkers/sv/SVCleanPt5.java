package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;

import htsjdk.variant.vcf.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKSVVariantContextUtils;

import java.util.*;

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
public class SVCleanPt5 extends MultiplePassVariantWalker {
    @Argument(
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "Output VCF name"
    )
    private GATKPath outputVcf;

    private VariantContextWriter vcfWriter;

    private final List<VariantContext> overlappingVariantsBuffer = new ArrayList<>();
    private final Set<String> filteredVariantIds = new HashSet<>();

    @Override
    protected int numberOfPasses() { return 2; }

    @Override
    public void onTraversalStart() {
        final VCFHeader header = getHeaderForVariants();
        final Set<VCFHeaderLine> originalHeaderLines = header.getMetaDataInInputOrder();

        // Add new header lines
        final Set<VCFHeaderLine> newHeaderLines = new HashSet<>();
        for (final VCFHeaderLine line : header.getMetaDataInInputOrder()) {
            if (line instanceof VCFInfoHeaderLine) {
                if (GATKSVVCFConstants.FILTER_VCF_INFO_LINES.contains(((VCFInfoHeaderLine) line).getID())) {
                    continue;
                }
            } else if (line instanceof VCFAltHeaderLine) {
                if (((VCFAltHeaderLine) line).getID().equals(GATKSVVCFConstants.UNR)) {
                    continue;
                }
            }
            if (GATKSVVCFConstants.FILTER_VCF_LINES.stream().anyMatch(line.toString()::contains)) {
                continue;
            }
            newHeaderLines.add(line);
        }

        // Write header
        VCFHeader newHeader = new VCFHeader(newHeaderLines, header.getGenotypeSamples());
        vcfWriter = createVCFWriter(outputVcf);
        vcfWriter.writeHeader(newHeader);
    }

    @Override
    public void closeTool() {
        if (vcfWriter != null) {
            vcfWriter.close();
        }
    }

    @Override
    protected void nthPassApply(final VariantContext variant, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext, int n) {
        switch (n) {
            case 0:
                firstPassApply(variant);
                break;
            case 1:
                secondPassApply(variant);
                break;
        }
    }

    @Override
    protected void afterNthPass(int n) {}

    public void firstPassApply(final VariantContext variant) {
        if (!variant.getFilters().contains(GATKSVVCFConstants.MULTIALLELIC)) {
            return;
        }

        overlappingVariantsBuffer.removeIf(vc -> !vc.getContig().equals(variant.getContig()) || vc.getEnd() < variant.getStart());
        for (VariantContext bufferedVariant : overlappingVariantsBuffer) {
            if (overlaps(bufferedVariant, variant)) {
                processVariantPair(bufferedVariant, variant);
                processVariantPair(variant, bufferedVariant);
            }
        }
        overlappingVariantsBuffer.add(variant);
    }

    public void secondPassApply(final VariantContext variant) {
        if (filteredVariantIds.contains(variant.getID())) {
            return;
        }

        VariantContextBuilder builder = new VariantContextBuilder(variant);
        processSvType(variant, builder);
        vcfWriter.add(builder.make());
    }

    private void processVariantPair(VariantContext largerVariant, VariantContext smallerVariant) {
        int lengthLarger = largerVariant.getEnd() - largerVariant.getStart() + 1;
        int lengthSmaller = smallerVariant.getEnd() - smallerVariant.getStart() + 1;
        if (lengthLarger < lengthSmaller) {
            return;
        }

        int overlapStart = Math.max(largerVariant.getStart(), smallerVariant.getStart());
        int overlapEnd = Math.min(largerVariant.getEnd(), smallerVariant.getEnd());
        int overlapLength = overlapEnd - overlapStart + 1;
        if (overlapLength <= 0) {
            return;
        }

        double smallCoverage = (double) overlapLength / lengthSmaller;
        if (smallCoverage > 0.5) {
            if (!filteredVariantIds.contains(largerVariant.getID())) {
                filteredVariantIds.add(smallerVariant.getID());
            }
        }
    }

    private void processSvType(final VariantContext variant, final VariantContextBuilder builder) {
        final String svType = variant.getAttributeAsString(GATKSVVCFConstants.SVTYPE, null);
        boolean hasMobileElement = variant.getAlleles().stream()
                .map(GATKSVVariantContextUtils::getSymbolicAlleleSymbols)
                .flatMap(Arrays::stream)
                .anyMatch(symbol -> symbol.equals(GATKSVVCFConstants.ME));
        if (svType == null || hasMobileElement) {
            return;
        }

        List<Genotype> genotypes = variant.getGenotypes();
        List<Genotype> updatedGenotypes = new ArrayList<>(genotypes.size());
        for (Genotype genotype : genotypes) {
            GenotypeBuilder gb = new GenotypeBuilder(genotype);
            gb.alleles(Arrays.asList(variant.getReference(), Allele.create("<" + svType + ">", false)));
            updatedGenotypes.add(gb.make());
        }

        final Allele refAllele = variant.getReference();
        final Allele altAllele = Allele.create("<" + svType + ">", false);
        List<Allele> newAlleles = Arrays.asList(refAllele, altAllele);
        builder.alleles(newAlleles);
        builder.genotypes(updatedGenotypes);
    }

    private boolean overlaps(final VariantContext v1, final VariantContext v2) {
        return v1.getContig().equals(v2.getContig()) && v1.getStart() <= v2.getEnd() && v2.getStart() <= v1.getEnd();
    }
}
