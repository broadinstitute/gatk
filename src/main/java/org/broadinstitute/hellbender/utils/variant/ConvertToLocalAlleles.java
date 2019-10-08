package org.broadinstitute.hellbender.utils.variant;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.GATKPathSpecifier;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;
import picard.cmdline.programgroups.VariantManipulationProgramGroup;

import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

@CommandLineProgramProperties(summary="Convert an input VCF with standard GT/AD/PL etc fields into local alleles for with LGT/LAA/LAD/LPL instead.",
        oneLineSummary="Convert VCF to local alleles format",
        programGroup= VariantManipulationProgramGroup.class)
@DocumentedFeature
public class ConvertToLocalAlleles extends VariantWalker {

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            optional = false,
            doc = "Output path to write VCF with locallized alleles to")
    public GATKPathSpecifier output;

    @Argument(fullName="keep-non-local", doc="Keep original non-local form of allele annotations in the output vcf.", optional = true)
    public boolean keepNonLocal= false;

    @Argument(fullName="all-output-local", doc="Convert to local form even if the local alleles are the same size as the original allele list", optional = true)
    public boolean allOutputLocal = false;

    public static final VCFFormatHeaderLine LGT_HEADER_LINE = new VCFFormatHeaderLine(LocalAlleler.LGT, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Like GT but with respect to the local alleles (LAA)");
    public static final VCFFormatHeaderLine LAA_HEADER_LINE = new VCFFormatHeaderLine(LocalAlleler.LAA, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, "Local Alleles, a list of the alleles that this genotype contains information about");
    public static final VCFFormatHeaderLine LAD_HEADER_LINE = new VCFFormatHeaderLine(LocalAlleler.LAD, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "AD but with respect to the local alleles");
    public static final VCFFormatHeaderLine LPL_HEADER_LINE = new VCFFormatHeaderLine(LocalAlleler.LPL, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, "PL but with respect to the local alleles");


    private VariantContextWriter writer;

    @Override
    public void onTraversalStart() {
        writer = createVCFWriter(output.toPath());
        writer.writeHeader(createHeader());
    }

    private VCFHeader createHeader() {
        final VCFHeader inputHeader = getHeaderForVariants();
        final Set<VCFHeaderLine> newMetaData = new LinkedHashSet<>(inputHeader.getMetaDataInInputOrder());
        if(inputHeader.getFormatHeaderLine(VCFConstants.GENOTYPE_KEY) != null) {
            newMetaData.add(LAA_HEADER_LINE);
            newMetaData.add(LGT_HEADER_LINE);

        }
        if(inputHeader.getFormatHeaderLine(VCFConstants.GENOTYPE_ALLELE_DEPTHS) != null) {
            newMetaData.add(LAD_HEADER_LINE);
        }
        if(inputHeader.getFormatHeaderLine(VCFConstants.GENOTYPE_PL_KEY) != null) {
            newMetaData.add(LPL_HEADER_LINE);
        }

        if(!keepNonLocal && allOutputLocal){
            newMetaData.removeIf( line -> {
                        if (line.getKey().equals("FORMAT")) {
                            final String id = ((VCFFormatHeaderLine) line).getID();
                            return (id.equals(VCFConstants.GENOTYPE_KEY)) || id.equals(VCFConstants.GENOTYPE_ALLELE_DEPTHS) || id.equals(VCFConstants.GENOTYPE_PL_KEY);
                        }
                        return false;
                    }
            );
        }
        return new VCFHeader(newMetaData, inputHeader.getSampleNamesInOrder());
    }

    @Override
    public void apply(final VariantContext variant, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext) {
        writer.add(localize(variant));
    }

    private VariantContext localize(final VariantContext vc) {
        final VariantContextBuilder vcb = new VariantContextBuilder(vc);
        final List<Genotype> localizedGenotypes = vc.getGenotypes()
                .stream()
                .map(genotype -> LocalAlleler.addLocalFields(genotype, vc, !keepNonLocal, !allOutputLocal))
                .collect(Collectors.toList());
        return vcb.genotypes(localizedGenotypes).make();
    }

    @Override
    public void closeTool() {
        if(writer != null){
            writer.close();
        }
    }
}
