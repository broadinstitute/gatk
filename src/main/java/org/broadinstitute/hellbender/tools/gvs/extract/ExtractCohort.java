package org.broadinstitute.hellbender.tools.gvs.extract;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.programgroups.ShortVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.gvs.common.*;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;

import java.util.*;



@CommandLineProgramProperties(
        summary = "(\"ExtractCohort\") - Filter and extract variants out of big query.",
        oneLineSummary = "Tool to extract variants out of big query for a subset of samples",
        programGroup = ShortVariantDiscoveryProgramGroup.class
)
@DocumentedFeature
public class ExtractCohort extends ExtractTool {
    private static final Logger logger = LogManager.getLogger(ExtractCohort.class);
    private ExtractCohortEngine engine;
    private SampleList sampleList;

    public enum VQSLODFilteringType { GENOTYPE, SITES, NONE }

    @Argument(
            fullName = "filter-set-info-table",
            doc = "Fully qualified name of the filtering set info table to use for cohort extraction",
            optional = true
    )
    private String filterSetInfoTableName = null;

    @Argument(
        fullName = "filter-set-site-table",
        doc = "Fully qualified name of the site filtering table to use for cohort extraction",
        optional = true
    )
    private String filterSetSiteTableName = null;

    @Argument(
            fullName = "tranches-table",
            doc = "Fully qualified name of the tranches table to use for cohort extraction",
            optional = true
    )
    private String tranchesTableName = null;

    @Argument(
            fullName = "cohort-extract-table",
            doc = "Fully qualified name of the table where the cohort data exists (already subsetted)",
            mutex = {"cohort-avro-file-name"},
            optional = true
    )
    private String cohortTable = null;

    @Argument(
            fullName = "cohort-avro-file-name",
            doc = "Path of the cohort avro file",
            mutex = {"cohort-extract-table"},
            optional = true
    )
    private GATKPath cohortAvroFileName = null;

    @Argument(
            fullName = "vet-ranges-fq-dataset",
            doc = "Fully qualified name for the dataset (<project>.<dataset>) that contains the VET and REF_RANGES data for extract",
            mutex = {"cohort-extract-table"},
            optional = true
    )
    private String vetRangesFQDataSet = null;


    @Argument(
            fullName = "vet-ranges-extract-fq-table",
            doc = "EXPERIMENTAL - ",
            optional = true
    )
    private String fqRangesExtractVetTable = null;

    @Argument(
            fullName = "ref-ranges-extract-fq-table",
            doc = "EXPERIMENTAL - ",
            optional = true
    )
    private String fqRangesExtractRefTable = null;

    @Argument(
            fullName = "vet-avro-file-name",
            doc = "Path to data from Vet table in Avro format",
            mutex = {"cohort-extract-table"},
            optional = true
    )
    private GATKPath vetAvroFileName = null;

    @Argument(
            fullName = "ref-ranges-avro-file-name",
            doc = "Path to data from Vet table in Avro format",
            mutex = {"cohort-extract-table"},
            optional = true
    )
    private GATKPath refRangesAvroFileName = null;

    @Argument(
            fullName = "presorted-avro-files",
            doc = "Indicates if Avro data is pre-sorted",
            mutex = {"cohort-extract-table"},
            optional = true
    )
    private boolean presortedAvroFiles = false;

    @Argument(
            fullName = "filter-set-name",
            doc = "Name in filter_set_name column of filtering table to use. Which training set should be applied in extract.",
            optional = true
    )
    private String filterSetName = null;

    @Argument(
            fullName = "emit-pls",
            doc = "Should PLs be emitted in output VCF",
            optional = true
    )
    private boolean emitPLs = false;

    @Argument(
            fullName = "emit-ads",
            doc = "Should allele depths be emitted in output VCF",
            optional = true
    )
    private boolean emitADs = false;

    // what if this was a flag input only?

    @Argument(
            fullName = "vqslod-filter-by-site",
            doc = "If VQSLOD filtering is applied, it should be at a site level. Default is false",
            optional = true
    )
    private boolean performSiteSpecificVQSLODFiltering = false;
    private VQSLODFilteringType vqslodfilteringType = VQSLODFilteringType.NONE;

    @Argument(
            fullName ="snps-truth-sensitivity-filter-level",
            mutex = {"snps-lod-score-cutoff"},
            doc = "The truth sensitivity level at which to start filtering SNPs",
            optional = true
    )
    private Double truthSensitivitySNPThreshold = null;

    @Argument(
            fullName = "indels-truth-sensitivity-filter-level",
            mutex = {"indels-lod-score-cutoff"},
            doc = "The truth sensitivity level at which to start filtering INDELs",
            optional = true
    )
    private Double truthSensitivityINDELThreshold = null;

    @Advanced
    @Argument(
            fullName = "snps-lod-score-cutoff",
            mutex = {"snps-truth-sensitivity-filter-level"},
            doc = "The VQSLOD score below which to start filtering SNPs",
            optional = true)
    private Double vqsLodSNPThreshold = null;

    @Advanced
    @Argument(
            fullName = "indels-lod-score-cutoff",
            mutex = {"indels-truth-sensitivity-filter-level"},
            doc = "The VQSLOD score below which to start filtering INDELs",
            optional = true)
    private Double vqsLodINDELThreshold = null;

    /**
    * If this flag is enabled, sites that have been marked as filtered (i.e. have anything other than '.' or 'PASS'
    * in the FILTER field) will be excluded from the output.
    */
    @Argument(
            fullName="exclude-filtered",
            doc="Don't include filtered sites in the final jointVCF",
            optional=true)
    private boolean excludeFilteredSites = false;

    @Argument(fullName = "inferred-reference-state",
            shortName = "irs",
            doc = "Reference state to be inferred from GVS, must match what was used during loading",
            optional = true)
    public GQStateEnum inferredReferenceState = GQStateEnum.SIXTY;

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

    protected static VCFHeader generateVcfHeader(Set<String> sampleNames,
                                                 final SAMSequenceDictionary sequenceDictionary,
                                                 final Set<VCFHeaderLine> extraHeaders) {
        final Set<VCFHeaderLine> headerLines = new HashSet<>();

        // Filter fields
        headerLines.add(GATKVCFHeaderLines.getFilterLine(GATKVCFConstants.LOW_QUAL_FILTER_NAME));
        headerLines.add(GATKVCFHeaderLines.getFilterLine(GATKVCFConstants.EXCESS_HET_KEY));
        headerLines.add(GATKVCFHeaderLines.getFilterLine(GATKVCFConstants.NO_HQ_GENOTYPES));

        // Info fields
        VCFStandardHeaderLines.addStandardInfoLines( headerLines, true,
                VCFConstants.ALLELE_COUNT_KEY,
                VCFConstants.ALLELE_FREQUENCY_KEY,
                VCFConstants.ALLELE_NUMBER_KEY
        );
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.AS_RAW_QUAL_APPROX_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.RAW_QUAL_APPROX_KEY));

        VCFStandardHeaderLines.addStandardFormatLines(headerLines, true,
                VCFConstants.GENOTYPE_KEY,
                VCFConstants.GENOTYPE_QUALITY_KEY
        );
        headerLines.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.REFERENCE_GENOTYPE_QUALITY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.AS_VQS_LOD_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.AS_YNG_STATUS_KEY));


        headerLines.addAll( extraHeaders );

        final VCFHeader header = new VCFHeader(headerLines, sampleNames);
        header.setSequenceDictionary(sequenceDictionary);

        return header;
    }

    /**
     * Enforce that if cost information is being recorded to the cost-observability-tablename then *all* recorded
     * parameters are set
     */
    @Override
    protected String[] customCommandLineValidation() {
        final List<String> errors = new ArrayList<>();
        if (costObservabilityTableName != null) {
            if (projectID == null || datasetID == null || callSetIdentifier == null || wdlStep == null || wdlCall == null || shardIdentifier == null) {
                errors.add("Parameters 'project-id', 'dataset-id', 'call-set-identifier', 'wdl-step', 'wdl-call', and 'shardIdentifier' must be set if 'cost-observability-tablename' is set.");
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

        Set<VCFHeaderLine> extraHeaderLines = new HashSet<>();

        if (filterSetInfoTableName != null) { // filter using vqslod-- default to GENOTYPE unless SITES specifically selected
            vqslodfilteringType = performSiteSpecificVQSLODFiltering ? VQSLODFilteringType.SITES : VQSLODFilteringType.GENOTYPE;
        }

        // filter at a site level (but not necesarily use vqslod)
        if ((filterSetSiteTableName != null && filterSetName == null) || (filterSetSiteTableName == null && filterSetName != null)) {
           throw new UserException("--filter-set-name and --filter-set-site-table are both necessary for any filtering related operations");
        }
        if (!vqslodfilteringType.equals(VQSLODFilteringType.NONE)) {
          if (filterSetInfoTableName == null || filterSetSiteTableName == null || filterSetName == null) {
            throw new UserException(" --filter-set-site-table, --filter-set-name and --filter-set-site-table are all necessary for any vqslod filtering operations");
          }
        }

        if (!vqslodfilteringType.equals(VQSLODFilteringType.NONE)) {
          FilterSensitivityTools.validateFilteringCutoffs(truthSensitivitySNPThreshold, truthSensitivityINDELThreshold, vqsLodSNPThreshold, vqsLodINDELThreshold, tranchesTableName);
          Map<String, Map<Double, Double>> trancheMaps = FilterSensitivityTools.getTrancheMaps(filterSetName, tranchesTableName, projectID);

          if (vqsLodSNPThreshold != null) { // we already have vqslod thresholds directly
            extraHeaderLines.add(FilterSensitivityTools.getVqsLodHeader(vqsLodSNPThreshold, GATKVCFConstants.SNP));
            extraHeaderLines.add(FilterSensitivityTools.getVqsLodHeader(vqsLodINDELThreshold, GATKVCFConstants.INDEL));
          } else { // using sensitivity threshold inputs; need to convert these to vqslod thresholds
            vqsLodSNPThreshold = FilterSensitivityTools.getVqslodThreshold(trancheMaps.get(GATKVCFConstants.SNP), truthSensitivitySNPThreshold, GATKVCFConstants.SNP);
            vqsLodINDELThreshold = FilterSensitivityTools.getVqslodThreshold(trancheMaps.get(GATKVCFConstants.INDEL), truthSensitivityINDELThreshold, GATKVCFConstants.INDEL);
            // set headers
            extraHeaderLines.add(FilterSensitivityTools.getTruthSensitivityHeader(truthSensitivitySNPThreshold, vqsLodSNPThreshold, GATKVCFConstants.SNP));
            extraHeaderLines.add(FilterSensitivityTools.getTruthSensitivityHeader(truthSensitivityINDELThreshold, vqsLodINDELThreshold, GATKVCFConstants.INDEL));
          }
        }

        if (vqslodfilteringType.equals(VQSLODFilteringType.GENOTYPE)) {
            extraHeaderLines.add(new VCFFormatHeaderLine("FT", 1, VCFHeaderLineType.String, "Genotype Filter Field"));
        }

        if (emitPLs) {
            VCFStandardHeaderLines.addStandardFormatLines(extraHeaderLines, true,
                    VCFConstants.GENOTYPE_PL_KEY
            );
        }

        if (emitADs) {
            VCFStandardHeaderLines.addStandardFormatLines(extraHeaderLines, true,
                    VCFConstants.GENOTYPE_ALLELE_DEPTHS
            );
        }

        sampleList = new SampleList(sampleTableName, sampleFileName, projectID, printDebugInformation, "extract-cohort");
        Map<Long, String> sampleIdToName = sampleList.getSampleIdToNameMap();

        VCFHeader header = generateVcfHeader(new HashSet<>(sampleIdToName.values()), reference.getSequenceDictionary(), extraHeaderLines);

        final List<SimpleInterval> traversalIntervals = getTraversalIntervals();

        if (minLocation == null && maxLocation == null && hasUserSuppliedIntervals()) {
            final SimpleInterval firstInterval = traversalIntervals.get(0);
            final SimpleInterval lastInterval = traversalIntervals.get(traversalIntervals.size() - 1);

            minLocation = SchemaUtils.encodeLocation(firstInterval.getContig(), firstInterval.getStart());
            maxLocation = SchemaUtils.encodeLocation(lastInterval.getContig(), lastInterval.getEnd());
        } else if ((minLocation != null || maxLocation != null) && hasUserSuppliedIntervals()) {
            throw new UserException("min-location and max-location should not be used together with intervals (-L).");
        }

        // if there is an avro file, the BQ specific parameters are unnecessary,
        // but they all are required if there is no avro file
        // KCIBUL: revisit!!!
//        if ((cohortAvroFileName == null && vetAvroFileName == null && refRangesAvroFileName == null) && (projectID == null || (cohortTable == null && vetRangesFQDataSet == null))) {
//            throw new UserException("Project id (--project-id) and cohort table (--cohort-extract-table) are required " +
//                "if no avro file (--cohort-avro-file-name or --vet-avro-file-name and --ref-ranges-avro-file-name) is provided.");
//        }

        // if there is a sample file, the BQ specific parameters are unnecessary,
        // but without a sample file, both a sample-table and a project-id are needed
        if (sampleFileName == null && (projectID == null || sampleTableName == null)) {
            throw new UserException("Project id (--project-id) and sample table (--sample-table) are required " +
                "if no sample file (--sample-file) is provided.");
        }

        engine = new ExtractCohortEngine(
                projectID,
                vcfWriter,
                header,
                annotationEngine,
                reference,
                sampleIdToName,
                cohortTable,
                cohortAvroFileName,
                vetRangesFQDataSet,
                fqRangesExtractVetTable,
                fqRangesExtractRefTable,
                vetAvroFileName,
                refRangesAvroFileName,
                traversalIntervals,
                minLocation,
                maxLocation,
                filterSetInfoTableName,
                filterSetSiteTableName,
                localSortMaxRecordsInRam,
                printDebugInformation,
                vqsLodSNPThreshold,
                vqsLodINDELThreshold,
                progressMeter,
                filterSetName,
                emitPLs,
                emitADs,
                vqslodfilteringType,
                excludeFilteredSites,
                inferredReferenceState,
                presortedAvroFiles);

        vcfWriter.writeHeader(header);
    }

    @Override
    // maybe think about creating a BigQuery Row walker?
    public void traverse() {
        progressMeter.setRecordsBetweenTimeChecks(100L);

        if ( filterSetInfoTableName == null || filterSetInfoTableName.equals("") ) {
            logger.warn("--filter-set-info-table is not specified, no filtering of cohort! ");
        }

        engine.traverse();
     }

    @Override
    public Object onTraversalSuccess() {
        if (costObservabilityTableName != null) {
            CostObservability costObservability = new CostObservability(projectID, datasetID, costObservabilityTableName);
            // Note - this is ONLY for the cost of querying the sample list.
            costObservability.writeCostObservability(callSetIdentifier, wdlStep, wdlCall, shardIdentifier,
                    new Date(), new Date(), "BigQuery Query Scanned",
                    sampleList.getBigQueryQueryByteScanned());
            costObservability.writeCostObservability(callSetIdentifier, wdlStep, wdlCall, shardIdentifier,
                    new Date(), new Date(), "Storage API Scanned",
                    engine.getTotalEstimatedBytesScanned());
        }
        return null;
    }

    @Override
    protected void onShutdown() {
        super.onShutdown();

        if ( engine != null ) {
            logger.info(String.format("***Processed %d total sites", engine.getTotalNumberOfSites()));
            logger.info(String.format("***Processed %d total variants", engine.getTotalNumberOfVariants()));
            logger.info(String.format("***Read API scanned %d bytes", engine.getTotalEstimatedBytesScanned()));
        }

        // Close up our writer if we have to:
        if ( vcfWriter != null ) {
            vcfWriter.close();
        }
    }

}