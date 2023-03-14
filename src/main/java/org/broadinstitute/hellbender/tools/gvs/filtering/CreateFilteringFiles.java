package org.broadinstitute.hellbender.tools.gvs.filtering;

import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ShortVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.gvs.common.ChromosomeEnum;
import org.broadinstitute.hellbender.tools.gvs.common.IngestConstants;
import org.broadinstitute.hellbender.tools.gvs.common.SchemaUtils;
import org.broadinstitute.hellbender.tools.gvs.ingest.CreateVariantIngestFiles;
import org.broadinstitute.hellbender.utils.tsv.SimpleXSVWriter;

import java.util.*;
import java.io.File;
import java.io.IOException;

/**
 * Ingest variant walker
 */
@CommandLineProgramProperties(
        summary = "Create files for ingesting VQSR Filtering data into BigQuery",
        oneLineSummary = "Filtering Data Ingest tool for BQJG",
        programGroup = ShortVariantDiscoveryProgramGroup.class,
        omitFromCommandLine = true
)
public final class CreateFilteringFiles extends VariantWalker {
    static final Logger logger = LogManager.getLogger(CreateVariantIngestFiles.class);

    private SimpleXSVWriter writer;

    private List<String> HEADER = 
        Arrays.asList("filter_set_name", "mode", "location", "ref", "alt", "vqslod", "culprit", "training_label", "yng");
    // TODO - not sure about culprit and training label.
    private List<String> HEADER_VQSR_LITE =
        Arrays.asList("filter_set_name", "mode", "location", "ref", "alt", "calibration_sensitivity", "score", "training_label", "yng");

    
    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, 
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, 
            doc = "Output file for filtering to be loaded into BigQuery", 
            optional = false)
    private File output;
    
    @Argument(
        fullName = "filter-set-name",
        doc = "Name to use in output files as the filter set name",
        optional = false)
    private String filterSetName;

    @Argument(
            fullName = "ref-version",
            doc = "Remove this option!!!! only for ease of testing. Valid options are 37 or 38",
            optional = true)
    private String refVersion = "37";

    @Argument(
        fullName = "mode",
        doc = "SNP or INDEL",
        optional = false)
    private String mode;

    @Argument(
            fullName = "classic",
            doc = "Whether or not this is using classic VQSR or the newer VQSR-Lite",
            optional = true)
    private Boolean usingOldVQSR = null;

    @Override
    public boolean requiresIntervals() {
        return false;
    }

    @Override
    public void onTraversalStart() {
        try {
            writer = new SimpleXSVWriter(output.toPath(), IngestConstants.SEPARATOR);
        } catch (IOException ioe) {
            throw new GATKException("Unable to initialize writer", ioe);
        }

        if (usingOldVQSR == null) { // default to using the old, or "classic" VQSR if the user specifies nothing
            usingOldVQSR = Boolean.TRUE;
        }

        if (usingOldVQSR) {
            writer.setHeaderLine(HEADER);
        } else {
            writer.setHeaderLine(HEADER_VQSR_LITE);
        }


        // Set reference version -- TODO remove this in the future, also, can we get ref version from the header?
        ChromosomeEnum.setRefVersion(refVersion);
    }

    @Override
    public void apply(final VariantContext variant, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext) {
        Long location = SchemaUtils.encodeLocation(variant.getContig(), variant.getStart());
        String ref = variant.getReference().getBaseString();

        if (variant.getAlternateAlleles().size() > 1) {
            throw new GATKException("Processing VCF with more than one alternate allele is disallowed: " + variant);
        }
        String alt = variant.getAlternateAllele(0).getBaseString();

        List<String> row;
        if (usingOldVQSR) {
            String vqslod = variant.getAttributeAsString("VQSLOD", "");
            String culprit = variant.getAttributeAsString("culprit", "");
            // TODO: check with Laura -- should NEGATIVES also be NAYs?
            String trainingLabel = variant.hasAttribute("POSITIVE_TRAIN_SITE") ? "POSITIVE" : (variant.hasAttribute("NEGATIVE_TRAIN_SITE") ? "NEGATIVE" : "");
            String yng = variant.hasAttribute("POSITIVE_TRAIN_SITE") ? "Y" : "G";
            row = Arrays.asList(
                    filterSetName,
                    mode,
                    location.toString(),
                    ref,
                    alt,
                    vqslod,
                    culprit,
                    trainingLabel,
                    yng
            );
        } else {
            // New VQSR-Lite has CALIBRATION_SENSITIVITY instead of vqslod
            String calibration_sensitivity = variant.getAttributeAsString("CALIBRATION_SENSITIVITY","");
            String score = variant.getAttributeAsString("SCORE","");
            String trainingLabel = variant.hasAttribute("training") ? "POSITIVE" : "";
            String yng = variant.hasAttribute("training") ? "Y" : "G";
            row = Arrays.asList(
                    filterSetName,
                    mode,
                    location.toString(),
                    ref,
                    alt,
                    calibration_sensitivity,
                    score,
                    trainingLabel,
                    yng
            );
        }

        writer.getNewLineBuilder().setRow(row).write();

    }

    @Override
    public void closeTool() {        
        if (writer != null) {
            try {
                writer.close();
            } catch (final IOException e) {
                throw new GATKException("Couldn't close writer", e);
            }
        }
    }

}
