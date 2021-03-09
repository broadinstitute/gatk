package org.broadinstitute.hellbender.tools.variantdb.arrays;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.ShortVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;
import org.broadinstitute.hellbender.tools.variantdb.*;
import org.broadinstitute.hellbender.utils.bigquery.TableReference;

import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.*;

/**
 * Ingest variant walker
 */
@CommandLineProgramProperties(
        summary = "Ingest tool for the Joint Genotyping in Big Query project",
        oneLineSummary = "Ingest tool for BQJG",
        programGroup = ShortVariantDiscoveryProgramGroup.class,
        omitFromCommandLine = true
)
public final class CreateImputedIngestFiles extends VariantWalker {
    static final Logger logger = LogManager.getLogger(CreateImputedIngestFiles.class);

    private Map<Integer, ImputedTsvCreator> tableToCreatorMap;
    private SampleList sampleNameMap;

    @Argument(fullName = "sample-list-table",
            shortName = "SLT",
            doc = "Fully qualified table name for the sample list table",
            optional = false)
    public String sampleListFQTablename;

    @Argument(
            fullName = "use-compressed-data",
            doc = "If true, use bit-packed fields for data",
            optional = true)
    private boolean useCompressedData = false;

    @Argument(
            fullName = "ref-version",
            doc = "Remove this option!!!! only for ease of testing. Valid options are 37 or 38",
            optional = true)
    private String refVersion = "37";

    @Argument(
            fullName = "print-debug-information",
            doc = "If true, print extra debugging output",
            optional = true)
    private boolean printDebugInformation = false;


    @Override
    public void onTraversalStart() {

        // Get sample name
        final VCFHeader inputVCFHeader = getHeaderForVariants();
        TableReference sampleTable = new TableReference(sampleListFQTablename, SchemaUtils.SAMPLE_FIELDS);
        sampleNameMap = new SampleList(sampleListFQTablename, null, null, printDebugInformation);
        tableToCreatorMap = new HashMap<>();

        Map<Integer, Set<String>> tableNumberToSampleList = new HashMap<>();

        // create the list of samples for each table number
        sampleNameMap.getMap().forEach((id, name) -> {
            int sampleTableNumber = IngestUtils.getTableNumber(id, IngestConstants.partitionPerTable);
            if (!tableNumberToSampleList.containsKey(sampleTableNumber)) {
                tableNumberToSampleList.put(sampleTableNumber, new HashSet<>());
            }
            tableNumberToSampleList.get(sampleTableNumber).add(name);
        });

        String runId = LocalDateTime.now().format(DateTimeFormatter.ofPattern("MMddyy_HHmm"));
        // create the writers for each table number
        tableNumberToSampleList.keySet().forEach(tableNumber -> tableToCreatorMap.put(tableNumber, new ImputedTsvCreator(getTableNumberPrefix(tableNumber), runId, tableNumberToSampleList.get(tableNumber), sampleNameMap, useCompressedData)));

        ChromosomeEnum.setRefVersion(refVersion);
    }


    public String getTableNumberPrefix(int number) {
        return String.format("%03d_", number);
    }

    @Override
    public void apply(final VariantContext variant, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext) {
        tableToCreatorMap.forEach((tableNumber, creator) -> {
            creator.apply(variant);
        });
    }

    @Override
    public Object onTraversalSuccess() {
        return 0;
    }

    @Override
    public void closeTool() {
        if (tableToCreatorMap == null || tableToCreatorMap.isEmpty()) {
            logger.info("tsvCreator is empty when closing tool");
        } else  {
            tableToCreatorMap.values().forEach(creator -> creator.closeTool());
        }
    }
}
