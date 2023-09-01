package org.broadinstitute.hellbender.tools.gvs.filtering;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ShortVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.gvs.common.*;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.gvs.bigquery.TableReference;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;

import java.util.*;

@SuppressWarnings("unused")
@CommandLineProgramProperties(
        summary = "(\"ExtractFeatures\") - Extract features data from BigQuery to train a filtering model.",
        oneLineSummary = "Tool to extract variants out of BigQuery to train a filtering model.",
        programGroup = ShortVariantDiscoveryProgramGroup.class
)
@DocumentedFeature
public class ExtractFeatures extends ExtractTool {
    private ExtractFeaturesEngine engine;
    protected ReferenceDataSource reference;
    private SampleList sampleList;

    private VariantContextWriter vcfWriter = null;

    @Argument(
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            doc = "Output VCF file to which annotated variants should be written."
    )
    private String outputVcfPathString = null;

    @Argument(
            fullName = "alt-allele-table",
            doc = "Fully qualified name of the table where the alternate allele info is"
    )
    private String fqAltAlleleTable = null;

    @Argument(
            fullName = "use-batch-queries",
            doc = "If true, use batch (rather than interactive) priority queries in BigQuery",
            optional = true)
    private boolean useBatchQueries = true;

    @Argument(
            fullName = "hq-genotype-gq-threshold",
            doc = "GQ threshold defining a high quality genotype",
            optional = true)
    private int hqGenotypeGQThreshold = 20;

    @Argument(
            fullName = "hq-genotype-depth-threshold",
            doc = "Depth threshold defining a high quality genotype",
            optional = true)
    private int hqGenotypeDepthThreshold = 10;

    @Argument(
            fullName = "hq-genotype-ab-threshold",
            doc = "Ab threshold defining a high quality genotype",
            optional = true)
    private double hqGenotypeABThreshold = 0.2;

    @Argument(
            fullName = "excess-alleles-threshold",
            doc = "Non-reference alleles threshold above which a site will be filtered out",
            optional = true)
    private int excessAllelesThreshold = CommonCode.EXCESS_ALLELES_THRESHOLD;

    @Argument(
            fullName = "query-labels",
            doc = "Key-value pairs to be added to the extraction BQ query. Ex: --query-labels key1=value1 --query-labels key2=value2",
            optional = true)
    private List<String> queryLabels = new ArrayList<>();

    @Argument(
            fullName = "cost-observability-tablename",
            doc = "Name of the bigquery table in which to store cost observability metadata",
            optional = true)
    private String costObservabilityTableName = null;

    @Argument(
            fullName = "call-set-identifier",
            doc = "Name of callset identifier, which is used to track costs in cost_observability table",
            optional = true)
    private String callSetIdentifier = null;

    @Argument(
            fullName = "wdl-step",
            doc = "Name of the WDL step/task (used for cost observability)",
            optional = true)
    private String wdlStep = null;

    @Argument(
            fullName = "wdl-call",
            doc = "Name of the call in the WDL step/task (used for cost observability)",
            optional = true)
    private String wdlCall = null;

    @Argument(
            fullName = "shard-identifier",
            doc = "Identifier for which shard was used in the WDL (used for cost observability)",
            optional = true)
    private String shardIdentifier = null;

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
        extraHeaderLines.add(GATKVCFHeaderLines.getFilterLine(GATKVCFConstants.LOW_QUAL_FILTER_NAME));

        reference = directlyAccessEngineReferenceDataSource();

        VCFHeader header = generateVcfHeader(
                new HashSet<>(), reference.getSequenceDictionary(), extraHeaderLines);

        vcfWriter = createVCFWriter(IOUtils.getPath(outputVcfPathString));

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
        if (vcfWriter != null) {
            vcfWriter.close();
        }
    }

    private static VCFHeader generateVcfHeader(Set<String> sampleNames,
                                              final SAMSequenceDictionary sequenceDictionary,
                                              final Set<VCFHeaderLine> extraHeaders) {
        final Set<VCFHeaderLine> headerLines = new HashSet<>();

        headerLines.addAll( getEvoquerVcfHeaderLines() );
        headerLines.addAll( extraHeaders );
        final VCFHeader header = new VCFHeader(headerLines, sampleNames);
        header.setSequenceDictionary(sequenceDictionary);

        return header;
    }

    private static Set<VCFHeaderLine> getEvoquerVcfHeaderLines() {
        final Set<VCFHeaderLine> headerLines = new HashSet<>();
        // Add standard VCF fields first:
        VCFStandardHeaderLines.addStandardInfoLines( headerLines, true,
                VCFConstants.STRAND_BIAS_KEY,
                VCFConstants.DEPTH_KEY,
                VCFConstants.RMS_MAPPING_QUALITY_KEY,
                VCFConstants.ALLELE_COUNT_KEY,
                VCFConstants.ALLELE_FREQUENCY_KEY,
                VCFConstants.ALLELE_NUMBER_KEY,
                VCFConstants.END_KEY
        );

        VCFStandardHeaderLines.addStandardFormatLines(headerLines, true,
                VCFConstants.GENOTYPE_KEY,
                VCFConstants.GENOTYPE_QUALITY_KEY,
                VCFConstants.DEPTH_KEY,
                VCFConstants.GENOTYPE_PL_KEY,
                VCFConstants.GENOTYPE_ALLELE_DEPTHS
        );

        headerLines.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.REFERENCE_GENOTYPE_QUALITY));

        headerLines.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_GT_KEY));
        headerLines.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_ID_KEY));

        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.AS_RAW_RMS_MAPPING_QUALITY_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.AS_RMS_MAPPING_QUALITY_KEY));

        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.AS_RAW_MAP_QUAL_RANK_SUM_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.AS_MAP_QUAL_RANK_SUM_KEY));

        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.AS_RAW_QUAL_APPROX_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.RAW_QUAL_APPROX_KEY));

        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.AS_RAW_READ_POS_RANK_SUM_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.AS_READ_POS_RANK_SUM_KEY));

        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.AS_SB_TABLE_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.SB_TABLE_KEY));

        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.AS_VARIANT_DEPTH_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.VARIANT_DEPTH_KEY));

        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.AS_STRAND_ODDS_RATIO_KEY));

        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.AS_FISHER_STRAND_KEY));

        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.AS_QUAL_BY_DEPTH_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.QUAL_BY_DEPTH_KEY));

        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.EXCESS_HET_KEY));

        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.SB_TABLE_KEY));

        headerLines.add(GATKVCFHeaderLines.getFilterLine(GATKVCFConstants.EXCESS_HET_KEY));
        headerLines.add(GATKVCFHeaderLines.getFilterLine(GATKVCFConstants.EXCESS_ALLELES));
        headerLines.add(GATKVCFHeaderLines.getFilterLine(GATKVCFConstants.NO_HQ_GENOTYPES));

        return headerLines;
    }

}
