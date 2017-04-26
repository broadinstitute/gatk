package org.broadinstitute.hellbender.tools.spark.sv;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.RequiredVariantInputArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariationSparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.engine.spark.datasources.VariantsSparkSource;
import org.broadinstitute.hellbender.utils.variant.GATKVariant;

import java.io.File;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Created by valentin on 4/20/17.
 */
@CommandLineProgramProperties(summary = "genotype SV variant call files",
        oneLineSummary = "genotype SV variant call files",
        programGroup = StructuralVariationSparkProgramGroup.class)
public class GenotypeStructuralVariantsSpark extends GATKSparkTool {

    private static final long serialVersionUID = 1L;

    @ArgumentCollection
    private RequiredVariantInputArgumentCollection variantArguments = new RequiredVariantInputArgumentCollection();

    @Argument(doc = "output VCF file",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            optional = true)
    private File outputFile = null;

    private VariantsSparkSource variantsSource;

    @Override
    public boolean requiresReads() {
        return true;
    }

    @Override
    public boolean requiresReference() {
        return true;
    }

    private void setUp(final JavaSparkContext ctx) {
        if (!outputFile.getParentFile().isDirectory()) {
            throw new CommandLineException.BadArgumentValue(StandardArgumentDefinitions.OUTPUT_SHORT_NAME, "the output file location is not a directory:");
        } else if (outputFile.exists() && !outputFile.isFile()) {
            throw new CommandLineException.BadArgumentValue(StandardArgumentDefinitions.OUTPUT_SHORT_NAME, "the output file makes reference to something that is not a file");
        }
        variantsSource = new VariantsSparkSource(ctx);
    }

    @Override
    protected void runTool(final JavaSparkContext ctx) {
        setUp(ctx);
        final JavaRDD<VariantContext> variants = variantsSource.getParallelVariantContexts(variantArguments.variantFiles.get(0).getFeaturePath(), getIntervals());
        final JavaRDD<VariantContext> outputVariants = processVariants(variants, ctx);
        final VCFHeader header = composeOutputHeader();
        SVVCFWriter.writeVCF(getAuthenticatedGCSOptions(), outputFile.getParent(), outputFile.getName(), referenceArguments.getReferenceFile().getAbsolutePath(), outputVariants, header, logger);
        tearDown(ctx);
    }

    private VCFHeader composeOutputHeader() {
        final SAMFileHeader readHeader = getHeaderForReads();
        final List<String> samples = readHeader.getReadGroups().stream()
                .map(SAMReadGroupRecord::getSample)
                .distinct()
                .sorted()
                .collect(Collectors.toList());
        final VCFHeader result = new VCFHeader(Collections.emptySet(), samples);
        result.setSequenceDictionary(getReferenceSequenceDictionary());
        result.addMetaDataLine(new VCFInfoHeaderLine(VCFConstants.END_KEY, 1, VCFHeaderLineType.Integer, "last base position of the variant"));
        return result;
    }

    private JavaRDD<VariantContext> processVariants(final JavaRDD<VariantContext> variants, final JavaSparkContext ctx) {
        return variants;
    }

    private void tearDown(final JavaSparkContext ctx) {
        variantsSource = null;
    }
}
