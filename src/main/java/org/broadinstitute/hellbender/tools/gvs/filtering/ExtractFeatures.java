package org.broadinstitute.hellbender.tools.gvs.filtering;

import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.programgroups.ShortVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.gvs.common.*;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.bigquery.TableReference;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;

import java.util.*;

@CommandLineProgramProperties(
        summary = "(\"ExtractFeatures\") - Extract features data from BQ to train filtering model.",
        oneLineSummary = "Tool to extract variants out of big query to train filtering model",
        programGroup = ShortVariantDiscoveryProgramGroup.class
)
@DocumentedFeature
public class ExtractFeatures extends ExtractTool {
    private static final Logger logger = LogManager.getLogger(ExtractFeatures.class);
    private ExtractFeaturesEngine engine;
    private SampleList sampleList;

    @Argument(
            fullName = "alt-allele-table",
            doc = "Fully qualified name of the table where the alternate allele info is",
            optional = false
    )
    protected String fqAltAlleleTable = null;

    @Argument(
        fullName = "use-batch-queries",
        doc = "If true, use batch (rather than interactive) priority queries in BigQuery",
        optional = true)
    protected boolean useBatchQueries = true;

    @Argument(
        fullName = "hq-genotype-gq-threshold",
        doc = "GQ threshold defining a high quality genotype",
        optional = true)
    protected int hqGenotypeGQThreshold = 20;

    @Argument(
        fullName = "hq-genotype-depth-threshold",
        doc = "Depth threshold defining a high quality genotype",
        optional = true)
    protected int hqGenotypeDepthThreshold = 10;

    @Argument(
        fullName = "hq-genotype-ab-threshold",
        doc = "Ab threshold defining a high quality genotype",
        optional = true)
    protected double hqGenotypeABThreshold = 0.2;

    @Argument(
        fullName = "excess-alleles-threshold",
        doc = "Non-reference alleles threshold above which a site will be filtered out",
        optional = true)
    protected int excessAllelesThreshold = CommonCode.EXCESS_ALLELES_THRESHOLD;

    @Argument(
        fullName = "query-labels",
        doc = "Key-value pairs to be added to the extraction BQ query. Ex: --query-labels key1=value1 --query-labels key2=value2",
        optional = true)
    protected List<String> queryLabels = new ArrayList<>();

    @Argument(
            fullName = "cost-observability-tablename",
            doc = "Name of the bigquery table in which to store cost observability metadata",
            optional = true)
    protected String costObservabilityTableName = null;

    @Argument(
            fullName = "call-set-identifier",
            doc = "Name of callset identifier, which is used to track costs in cost_observability table",
            optional = true)
    protected String callSetIdentifier = null;

    @Argument(
            fullName = "wdl-step",
            doc = "Name of the WDL step/task (used for cost observability)",
            optional = true)
    protected String wdlStep = null;

    @Argument(
            fullName = "wdl-call",
            doc = "Name of the call in the WDL step/task (used for cost observability)",
            optional = true)
    protected String wdlCall = null;

    @Argument(
            fullName = "shard-identifier",
            doc = "Identifier for which shard was used in the WDL (used for cost observability)",
            optional = true)
    protected String shardIdentifier = null;

    @Override
    public boolean requiresIntervals() {
        return false;
    }

    /**
     * Enforce that if cost information is being recorded to the cost-observability-tablename then *all* recorded
     * parameters are set
     */
    @Override
    protected String[] customCommandLineValidation() {
        final List<String> errors = new ArrayList<>();
        if (projectID != null || datasetID != null || costObservabilityTableName != null || callSetIdentifier != null || wdlStep != null || wdlCall != null || shardIdentifier != null) {
            if (projectID == null || datasetID == null || costObservabilityTableName == null || callSetIdentifier == null || wdlStep == null || wdlCall == null || shardIdentifier == null) {
                errors.add("Parameters 'project-id', 'dataset-id', 'cost-observability-tablename', 'call-set-identifier', 'wdl-step', 'wdl-call', and 'shardIdentifier' must either ALL be set or ALL NOT BE set");
            }
        }
        if (!errors.isEmpty()) {
            return errors.toArray(new String[0]);
        }
        return null;
    }

    @Override
    protected void onStartup() {
        super.onStartup();

        TableReference sampleTableRef = new TableReference(sampleTableName, SchemaUtils.SAMPLE_FIELDS);
        sampleList = new SampleList(sampleTableName, sampleFileName, projectID, printDebugInformation, "extract-features");

        Set<VCFHeaderLine> extraHeaderLines = new HashSet<>();
        extraHeaderLines.add(
            FilterSensitivityTools.getExcessAllelesHeader(excessAllelesThreshold, GATKVCFConstants.EXCESS_ALLELES));

        extraHeaderLines.add(GATKVCFHeaderLines.getFilterLine(GATKVCFConstants.LOW_QUAL_FILTER_NAME));

        VCFHeader header = CommonCode.generateVcfHeader(
            new HashSet<>(), reference.getSequenceDictionary(), extraHeaderLines);

        final List<SimpleInterval> traversalIntervals = getTraversalIntervals();

        if (minLocation == null && maxLocation == null && hasUserSuppliedIntervals()) {
            final SimpleInterval firstInterval = traversalIntervals.get(0);
            final SimpleInterval lastInterval = traversalIntervals.get(traversalIntervals.size() - 1);

            minLocation = SchemaUtils.encodeLocation(firstInterval.getContig(), firstInterval.getStart());
            maxLocation = SchemaUtils.encodeLocation(lastInterval.getContig(), lastInterval.getEnd());
        } else if ((minLocation != null || maxLocation != null) && hasUserSuppliedIntervals()) {
            throw new UserException("min-location and max-location should not be used together with intervals (-L).");
        }

        engine = new ExtractFeaturesEngine(
            projectID,
            datasetID,
            vcfWriter,
            header,
            annotationEngine,
            reference,
            fqAltAlleleTable,
            sampleTableRef,
            traversalIntervals,
            minLocation,
            maxLocation,
            localSortMaxRecordsInRam,
            printDebugInformation,
            useBatchQueries,
            progressMeter,
            sampleList.size(),
            hqGenotypeGQThreshold,
            hqGenotypeDepthThreshold,
            hqGenotypeABThreshold,
            excessAllelesThreshold,
            queryLabels);

        vcfWriter.writeHeader(header);
}

    @Override
    // maybe think about creating a BigQuery Row walker?
    public void traverse() {
        progressMeter.setRecordsBetweenTimeChecks(100L);

        engine.traverse();
    }

    @Override
    public Object onTraversalSuccess() {
        if (costObservabilityTableName != null) {
            CostObservability costObservability = new CostObservability(projectID, datasetID, costObservabilityTableName);
            costObservability.writeCostObservability(callSetIdentifier, wdlStep, wdlCall, shardIdentifier,
                    new Date(), new Date(), "BigQuery Query Scanned",
                    sampleList.getBigQueryQueryByteScanned() + engine.getBigQueryQueryBytesScanned());
            costObservability.writeCostObservability(callSetIdentifier, wdlStep, wdlCall, shardIdentifier,
                    new Date(), new Date(), "Storage API Scanned",
                    engine.getStorageAPIBytesScanned());
        }
        return null;
    }

    @Override
    protected void onShutdown() {
        super.onShutdown();

        // Close up our writer if we have to:
        if ( vcfWriter != null ) {
            vcfWriter.close();
        }
    }

}
