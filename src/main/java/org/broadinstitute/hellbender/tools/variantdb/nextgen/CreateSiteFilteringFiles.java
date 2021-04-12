package org.broadinstitute.hellbender.tools.variantdb.nextgen;

import htsjdk.variant.variantcontext.VariantContext;
import io.grpc.netty.shaded.io.netty.util.internal.StringUtil;

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

@CommandLineProgramProperties(
        summary = "Create files for ingesting Site-level filtering data into BigQuery",
        oneLineSummary = "Site-level filtering Data Ingest tool for GVS",
        programGroup = ShortVariantDiscoveryProgramGroup.class,
        omitFromCommandLine = true
)
public final class CreateSiteFilteringFiles extends VariantWalker {
    static final Logger logger = LogManager.getLogger(CreateSiteFilteringFiles.class);

    private SimpleXSVWriter writer;

    private List<String> lastEntryProcessed = new ArrayList<>();

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
        writer.setHeaderLine(SchemaUtils.FILTER_SET_SITE_FIELDS);

        // Set reference version -- TODO remove this in the future, also, can we get ref version from the header?
        ChromosomeEnum.setRefVersion(refVersion);
    }

    @Override
    public void apply(final VariantContext variant, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext) {
        if (variant.isFiltered()) {
            Long location = SchemaUtils.encodeLocation(variant.getContig(), variant.getStart());
            String filters = StringUtil.join(",", variant.getFilters()).toString();
            
            List<String> row = Arrays.asList(
                filterSetName,
                location.toString(),
                filters
            );

            // the input file may have multiple lines for different alleles,
            // however the filters will be the same.  This simple prevents us from
            // writing out essentially duplicate rows to the TSV
            if (!row.equals(lastEntryProcessed)) {
                writer.getNewLineBuilder().setRow(row).write();
            }
            lastEntryProcessed = row;
        }
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
