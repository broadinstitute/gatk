package org.broadinstitute.hellbender.tools.walkers.validation;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.apache.commons.math3.util.MathArrays;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;
import org.broadinstitute.hellbender.utils.MathUtils;

import java.io.File;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * Given mixing weights of different samples in a pooled bam, annotate a corresponding
 * vcf containing individual sample genotypes.
 *
 * <p>The formula is:</p>
 *
 * <p>
 *     Expected allele fraction = SUM_samples {mixing_fraction(sample) * [0 if hom ref, 0.5 is het, 1.0 if hom var]}
 * </p>
 *
 * <h3>Example</h3>
 *
 * <pre>
 * java -jar gatk.jar AnnotateVcfWithExpectedAlleleFraction \
 *   -V input.vcf \
 *   -O output.vcf \
 *   -mixingFractions mixingFractions.table
 * </pre>
 *
 * Created by David Benjamin on 1/31/17.
 */
@CommandLineProgramProperties(
        summary = "Annotate a multi-sample vcf with expected allele fractions in pooled sequencing given mixing" +
                " fractions of the different samples in the pool via the formula" +
                " Expected allele fraction = SUM_samples {mixing_fraction(sample) * [0 if hom ref, 0.5 is het, 1.0 if hom var]}",
        oneLineSummary = "(Internal) Annotate a vcf with expected allele fractions in pooled sequencing",
        programGroup = VariantProgramGroup.class
)
@DocumentedFeature
public class AnnotateVcfWithExpectedAlleleFraction extends VariantWalker {

    public static final String MIXING_FRACTIONS_TABLE_NAME = "mixingFractions";

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "The output annotated VCF file",
            optional = false)
    private final File outputVcf = null;

    @Argument(fullName = MIXING_FRACTIONS_TABLE_NAME,
            shortName = MIXING_FRACTIONS_TABLE_NAME,
            doc = "The input mixing fractions table",
            optional = false)
    private final File inputMixingFractions = null;

    private VariantContextWriter vcfWriter;

    public static final String EXPECTED_ALLELE_FRACTION_NAME = "AF_EXP";

    double[] mixingFractionsInSampleOrder;

    @Override
    public void onTraversalStart() {
        final VCFHeader inputHeader = getHeaderForVariants();
        final Set<VCFHeaderLine> headerLines = new HashSet<>(inputHeader.getMetaDataInSortedOrder());
        headerLines.add(new VCFInfoHeaderLine(EXPECTED_ALLELE_FRACTION_NAME, 1, VCFHeaderLineType.Float, "expected allele fraction in pooled bam"));
        final VCFHeader vcfHeader = new VCFHeader(headerLines, inputHeader.getGenotypeSamples());
        vcfWriter = createVCFWriter(outputVcf);
        vcfWriter.writeHeader(vcfHeader);

        final List<MixingFraction> mixingFractionsList = MixingFraction.readMixingFractions(inputMixingFractions);
        final Map<String, Double> mixingfractionsMap = mixingFractionsList.stream()
                .collect(Collectors.toMap(MixingFraction::getSample, MixingFraction::getMixingFraction));
        mixingFractionsInSampleOrder = inputHeader.getSampleNamesInOrder().stream()
                .mapToDouble(mixingfractionsMap::get).toArray();
    }

    @Override
    public void apply(final VariantContext vc, final ReadsContext readsContext, final ReferenceContext refContext, final FeatureContext fc) {
        final double[] weights = vc.getGenotypes().stream().mapToDouble(g -> weight(g)).toArray();
        final double expectedAlleleFraction = MathUtils.sum(MathArrays.ebeMultiply(weights, mixingFractionsInSampleOrder));
        vcfWriter.add(new VariantContextBuilder(vc).attribute(EXPECTED_ALLELE_FRACTION_NAME, expectedAlleleFraction).make());
    }

    @Override
    public void closeTool() {
        if ( vcfWriter != null ) {
            vcfWriter.close();
        }
    }

    private static double weight(final Genotype genotype) {
        if (genotype.isHomVar()) {
            return 1.0;
        } else if (genotype.isHet()) {
            return 0.5;
        } else {
            return 0.0;
        }
    }
}
