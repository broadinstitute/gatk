package org.broadinstitute.hellbender.tools.variantdb.nextgen;

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
import org.broadinstitute.hellbender.tools.variantdb.*;
import org.broadinstitute.hellbender.tools.variantdb.IngestConstants;
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
        Arrays.asList("filter_set_name","mode","location","ref","alt","vqslod","culprit","training_label","yng");

    
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
        writer.setHeaderLine(HEADER);

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

        String vqslod = variant.getAttributeAsString("VQSLOD","");
        String culprit = variant.getAttributeAsString("culprit","");
        String trainingLabel = variant.hasAttribute("POSITIVE_TRAIN_SITE")?"POSITIVE":(variant.hasAttribute("NEGATIVE_TRAIN_SITE")?"NEGATIVE":"");

        // TODO: check with Laura -- should NEGATIVES also be NAYs?
        String yng = variant.hasAttribute("POSITIVE_TRAIN_SITE")?"Y":"G";

        List<String> row = Arrays.asList(
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
