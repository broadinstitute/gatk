package org.broadinstitute.hellbender.tools;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.FlowBasedProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.*;

@CommandLineProgramProperties(
        summary = "Divides annotations that were summed across samples by genomicsDB by the number of samples with het or hom var calls. " +
                "This is an approximation of taking the average of these annotations.",
        oneLineSummary = "Divides annotations that were summed by genomicsDB by number of samples to calculate average.",
        programGroup = FlowBasedProgramGroup.class
)
@ExperimentalFeature
public final class CalculateAverageCombinedAnnotations extends VariantWalker {
    public static final String ANNOTATION_LIST_LONG_NAME = "summed-annotation-to-divide";
    public static final String ANNOTATION_LIST_SHORT_NAME = "summed-annotation";

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "Output file (if not provided, defaults to STDOUT)", common = false, optional = false)
    private GATKPath outputFile = null;

    @Argument(fullName = ANNOTATION_LIST_LONG_NAME, shortName = ANNOTATION_LIST_SHORT_NAME, doc = "INFO Annotations in VCF that have been summed by GenomicsDB and need to be divided by the number of het or homvar samples to calculate the average value. Must use annotation string as it's defined in the VCF.", common = false, optional = false)
    private List<String> annotations = new ArrayList<>();

    private VariantContextWriter vcfWriter = null;

    @Override
    public void onTraversalStart() {
        vcfWriter = createVCFWriter(outputFile);
        if (annotations.size() == 0) {
            throw new UserException("--" + ANNOTATION_LIST_LONG_NAME + " must be provided.");
        }
        VCFHeader header = getHeaderForVariants();
        for(String annot : annotations) {
            header.addMetaDataLine(new VCFInfoHeaderLine("AVERAGE_" + annot, 1, VCFHeaderLineType.Float, "Average of "+ annot +" annotation across samples. See "+ annot +" header line for more information."));
        }
        vcfWriter.writeHeader(header);
    }

    @Override
    public void apply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        if (!variant.hasAttribute(GATKVCFConstants.RAW_GENOTYPE_COUNT_KEY)){
            throw new UserException(String.format("Need annotation %s at site %s:%d", GATKVCFConstants.RAW_GENOTYPE_COUNT_KEY, variant.getContig(), variant.getStart()));
        }
        List<String> genotypeCounts = variant.getAttributeAsStringList(GATKVCFConstants.RAW_GENOTYPE_COUNT_KEY, "");
        double counter = Double.parseDouble(genotypeCounts.get(1)) + Double.parseDouble(genotypeCounts.get(2)); //Het and hom var counts, all alleles are lumped together.
        if (counter > 0 ) {
            Map<String, Object> finalAnnotations = new HashMap<>();
            for (String annot : annotations) {
                if (variant.hasAttribute(annot)) {
                    finalAnnotations.put("AVERAGE_" + annot, variant.getAttributeAsDouble(annot, 0) / counter);
                }
            }
            VariantContext vc = new VariantContextBuilder(variant).putAttributes(finalAnnotations).make();
            vcfWriter.add(vc);
        } else {
            vcfWriter.add(variant);
        }
    }

    @Override
    public void closeTool() {
        vcfWriter.close();
    }
}
