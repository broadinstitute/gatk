package org.broadinstitute.hellbender.tools.walkers.variantutils;

import htsjdk.tribble.Feature;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.*;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.ReadThreadingAssembler;
import org.broadinstitute.hellbender.utils.codecs.table.TableFeature;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.tsv.SimpleXSVWriter;
import picard.cmdline.programgroups.VariantEvaluationProgramGroup;

import java.io.File;
import java.io.IOException;
import java.util.*;


/**
 * This tool exists to overcome a limitation of the HaplotypeCaller '--assembly-region-out' argument in that it doesn't
 * track the variants themselves. For some debugging/assembly engine evaluation work it is helpful to annotate the overlapping
 * regions based on some subsetted VCF output. THis tool allows for that to happen.
 */
@CommandLineProgramProperties(
        summary = "Annotates an AssemblyRegionOutput file from the haplotype caller with a count of overlapping variants.",
        oneLineSummary = "Annotate AssemblyRegionOutput",
        programGroup = VariantEvaluationProgramGroup.class
)

@ExperimentalFeature
public class SummarizeActiveRegionOutAgainstVCF extends FeatureWalker<TableFeature> {
    private static final Logger logger = LogManager.getLogger(SummarizeActiveRegionOutAgainstVCF.class);

    @Argument(fullName = "active-region-summary", doc = "table output of active regions from haplotype caller or mutect")
    File input = null;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "output tsv to which to write the summary")
    String outputPath = null;
    private SimpleXSVWriter outputTableWriter;

    @Argument(fullName = StandardArgumentDefinitions.VARIANT_LONG_NAME,
            shortName = StandardArgumentDefinitions.VARIANT_SHORT_NAME,
            doc="Variant file to use for annotating file", optional=true)
    public FeatureInput<VariantContext> overlappingVariantInput;

    @Override
    public void onTraversalStart() {
        super.onTraversalStart();
        try {
            this.outputTableWriter = new SimpleXSVWriter(IOUtils.getPath(outputPath), '\t');
            List<String> newHeader = new ArrayList<>(ReadThreadingAssembler.DEBUG_ACTIVE_REGION_OUT_HEADER_LINES);
            newHeader.add("variants_overlapping");
            outputTableWriter.setHeaderLine(newHeader);
        } catch (IOException e) {
            throw new GATKException("failed to open output writer");
        }
    }

    @Override
    public Object onTraversalSuccess() {
        Object ret = super.onTraversalSuccess();
        try {
            outputTableWriter.close();
        } catch (IOException outputException) {
            throw new GATKException("Failed when trying to close the output table", outputException);
        }
        return ret;
    }

    @Override
    protected boolean isAcceptableFeatureType(Class<? extends Feature> featureType) {
        return featureType.equals(TableFeature.class);
    }

    @Override
    // Copy all the values of the line and append our extra row.
    public void apply(TableFeature feature, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        List<VariantContext> overlappingVariants = featureContext.getValues(overlappingVariantInput);

        SimpleXSVWriter.LineBuilder builder = outputTableWriter.getNewLineBuilder();
        for (String key : ReadThreadingAssembler.DEBUG_ACTIVE_REGION_OUT_HEADER_LINES) {
            builder.setColumn(key, feature.get(key));
        }
        builder.setColumn("variants_overlapping", Integer.toString(overlappingVariants.size())).write();
    }

    @Override
    public File getDrivingFeatureFile() {
        return input;
    }
}