package org.broadinstitute.hellbender.tools.walkers.variantutils;

import com.google.common.collect.Lists;
import htsjdk.tribble.Feature;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberStandardArgument;
import org.broadinstitute.hellbender.tools.copynumber.utils.annotatedinterval.AnnotatedInterval;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.DataSourceUtils;
import org.broadinstitute.hellbender.tools.funcotator.metadata.FuncotationMetadata;
import org.broadinstitute.hellbender.tools.funcotator.metadata.VcfFuncotationMetadata;
import org.broadinstitute.hellbender.transformers.VariantTransformer;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.codecs.table.TableFeature;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.tsv.SimpleXSVWriter;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;
import picard.cmdline.programgroups.VariantEvaluationProgramGroup;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.util.*;
import java.util.stream.Collectors;


@CommandLineProgramProperties(
        summary = "Validates a VCF file with an extra strict set of criteria.",
        oneLineSummary = "Validate VCF",
        programGroup = VariantEvaluationProgramGroup.class
)
@BetaFeature
public class SummarizeActiveRegionOutAgainstVCF extends FeatureWalker<TableFeature> {
    private static final Logger logger = LogManager.getLogger(SummarizeActiveRegionOutAgainstVCF.class);

    @Argument(fullName = "active-region-summary", doc = "table output of active regions from haplotype caller or mutect")
    File input = null;

    @Argument(fullName = "output", doc = "output tsv to which to write the summary")
    String outputPath = null;
    private SimpleXSVWriter outputTableWriter;


    @Argument(fullName= "variants-to-count", doc="The set of alleles to force-call regardless of evidence", optional=true)
    public FeatureInput<VariantContext> overlappingVariantInput;

    @Override
    public void onTraversalStart() {
        super.onTraversalStart();
        try {
            this.outputTableWriter = new SimpleXSVWriter(IOUtils.getPath(outputPath), '\t');
            outputTableWriter.setHeaderLine(Lists.newArrayList("active_region", "kmers_assembled", "haplotypes_found", "variants_overlapping"));
        } catch (IOException e) {
            throw new GATKException("failed to open output writer");
        }
    }

    @Override
    protected boolean isAcceptableFeatureType(Class<? extends Feature> featureType) {
        return featureType.equals(TableFeature.class);
    }

    @Override
    public void apply(TableFeature feature, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        List<VariantContext> overlappingVariants = featureContext.getValues(overlappingVariantInput);

        outputTableWriter.getNewLineBuilder()
                .setColumn("active_region", feature.get("active_region"))
                .setColumn("kmers_assembled", feature.get("kmers_assembled"))
                .setColumn("haplotypes_found", feature.get("haplotypes_found"))
                .setColumn("variants_overlapping", Integer.toString(overlappingVariants.size()));
    }

    @Override
    public File getDrivingFeatureFile() {
        return input;
    }
}