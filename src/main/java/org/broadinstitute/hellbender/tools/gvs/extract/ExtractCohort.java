package org.broadinstitute.hellbender.tools.gvs.extract;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.*;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.gvs.common.*;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;

import java.util.*;


@DocumentedFeature
public abstract class ExtractCohort extends ExtractTool {
    protected static final Logger logger = LogManager.getLogger(ExtractCohort.class);
    protected ExtractCohortEngine engine;
    protected ReferenceDataSource reference;
    private SampleList sampleList;

    public enum VQScoreFilteringType {GENOTYPE, SITES, NONE}

    protected VCFHeader header;

    @Argument(
            fullName = "filter-set-info-table",
            doc = "Fully qualified name of the filter set info table to use for cohort extraction",
            optional = true
    )
    private String filterSetInfoTableName = null;

    @Argument(
            fullName = "filter-set-site-table",
            doc = "Fully qualified name of the filter set site table to use for cohort extraction",
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
            fullName = "vet-ranges-fq-dataset",
            doc = "Fully qualified name for the dataset (<project>.<dataset>) that contains the VET and REF_RANGES data for extract",
            optional = true
    )
    private String vetRangesFQDataSet = null;


    @Argument(
            fullName = "vet-ranges-extract-fq-table",
            doc = "Fully qualified name for the VET table prepared for extract.",
            optional = true
    )
    private String fqRangesExtractVetTable = null;

    @Argument(fullName = "vet-ranges-extract-table-version",
            doc = "Version of the vet ranges extract table - for maintaining backwards-compatibility",
            optional = true)
    private VetRangesExtractVersionEnum vetRangesExtractTableVersion = VetRangesExtractVersionEnum.V2;


    @Argument(
            fullName = "ref-ranges-extract-fq-table",
            doc = "Fully qualified name for the REF_RANGES table prepared for extract.",
            optional = true
    )
    private String fqRangesExtractRefTable = null;

    @Argument(
            fullName = "vet-avro-file-name",
            doc = "Path to data from Vet table in Avro format",
            optional = true
    )
    private GATKPath vetAvroFileName = null;

    @Argument(
            fullName = "ref-ranges-avro-file-name",
            doc = "Path to data from Vet table in Avro format",
            optional = true
    )
    private GATKPath refRangesAvroFileName = null;

    @Argument(
            fullName = "presorted-avro-files",
            doc = "Indicates if Avro data is pre-sorted",
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

    @Argument(
            fullName = "use-vqsr-scoring",
            doc = "If true, use VQSR scoring (vqs Lod score). Otherwise use VETS scoring (calibration_sensitivity)",
            optional = true
    )
    private boolean isVQSR = false;

    @Argument(
            fullName = "vqs-score-filter-by-site",
            doc = "If Variant Quality Score filtering (either VETS or VQSR) is applied, it should be at a site level. Default is false",
            optional = true
    )
    // historical note that this parameter was previously named 'vqsr-score-filter-by-site', changed as it's not VQSR-specific
    private boolean performSiteSpecificVQScoreFiltering = false;
    private VQScoreFilteringType vqScoreFilteringType = VQScoreFilteringType.NONE;

    @Argument(
            fullName = "convert-filtered-genotypes-to-no-calls",
            doc = "Set filtered genotypes to no-calls. This option can only be used if Variant Quality Score filtering " +
                    "is applied at the genotype level",
            optional = true
    )
    private boolean convertFilteredGenotypesToNoCalls = false;

    @Argument(
            fullName = "maximum-alternate-alleles",
            doc = "The maximum number of alternate alleles a site can have before it is hard-filtered from the output",
            optional = true
    )
    private Long maximumAlternateAlleles = 0L;

    @Argument(
            fullName = "snps-truth-sensitivity-filter-level",
            mutex = {"snps-lod-score-cutoff"},
            doc = "The truth sensitivity level at which to start filtering SNPs (0-100)",
            optional = true
    )
    private Double truthSensitivitySNPThreshold = null;

    @Argument(
            fullName = "indels-truth-sensitivity-filter-level",
            mutex = {"indels-lod-score-cutoff"},
            doc = "The truth sensitivity level at which to start filtering INDELs (0-100)",
            optional = true
    )
    private Double truthSensitivityINDELThreshold = null;

    @Advanced
    @Argument(
            fullName = "snps-lod-score-cutoff",
            mutex = {"snps-truth-sensitivity-filter-level"},
            doc = "The VQSLOD score below which to start filtering SNPs. For VQSR Classic Mode ONLY.",
            optional = true)
    private Double vqsLodSNPThreshold = null;

    @Advanced
    @Argument(
            fullName = "indels-lod-score-cutoff",
            mutex = {"indels-truth-sensitivity-filter-level"},
            doc = "The VQSLOD score below which to start filtering INDELs. For VQSR Classic Mode ONLY.",
            optional = true)
    private Double vqsLodINDELThreshold = null;

    /**
     * If this flag is enabled, sites that have been marked as filtered (i.e. have anything other than '.' or 'PASS'
     * in the FILTER field) will be excluded from the output.
     */
    @Argument(
            fullName = "exclude-filtered",
            doc = "Don't include filtered sites in the final jointVCF",
            optional = true)
    protected boolean excludeFilteredSites = false;

    @Argument(fullName = "inferred-reference-state",
            shortName = "irs",
            doc = "Reference state to be inferred from GVS, must match what was used during loading",
            optional = true)
    private GQStateEnum inferredReferenceState = GQStateEnum.SIXTY;

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

    private static VCFHeader generateVcfHeader(Set<String> sampleNames,
                                                 final SAMSequenceDictionary sequenceDictionary,
                                                 final Set<VCFHeaderLine> extraHeaders) {
        final Set<VCFHeaderLine> headerLines = new HashSet<>();

        // Filter fields
        headerLines.add(GATKVCFHeaderLines.getFilterLine(GATKVCFConstants.LOW_QUAL_FILTER_NAME));
        headerLines.add(GATKVCFHeaderLines.getFilterLine(GATKVCFConstants.EXCESS_ALLELES));
        headerLines.add(GATKVCFHeaderLines.getFilterLine(GATKVCFConstants.EXCESS_HET_KEY));
        headerLines.add(GATKVCFHeaderLines.getFilterLine(GATKVCFConstants.NO_HQ_GENOTYPES));

        // Info fields
        VCFStandardHeaderLines.addStandardInfoLines(headerLines, true,
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

        headerLines.addAll(extraHeaders);

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
        if (truthSensitivitySNPThreshold != null) {
            if ((truthSensitivitySNPThreshold <= 0.0) || (truthSensitivityINDELThreshold >= 100.0)) {
                errors.add("Parameter 'snps-truth-sensitivity-filter-level' must be between > 0.0 and < 100.0");
            }
        }
        if (truthSensitivityINDELThreshold != null) {
            if ((truthSensitivityINDELThreshold <= 0.0) || (truthSensitivityINDELThreshold >= 100.0)) {
                errors.add("Parameter 'indels-truth-sensitivity-filter-level' must be between > 0.0 and < 100.0");
            }
        }
        if (!isVQSR) {
            if (tranchesTableName != null) {
                errors.add("Parameter 'tranches-table' is not allowed for VETS");
            }
            if ((vqsLodSNPThreshold != null) || (vqsLodINDELThreshold != null)) {
                errors.add("Parameters 'snps-lod-score-cutoff' and 'indels-lod-score-cutoff' cannot be used in VETS mode");
            }
        }
        if (!errors.isEmpty()) {
            return errors.toArray(new String[0]);
        }
        return null;
    }

    protected abstract void apply(VariantContext variantContext);

    @Override
    protected void onStartup() {
        super.onStartup();

        if (filterSetInfoTableName != null) { // filter using VQScore (vqslod or calibration_sensitivity) -- default to GENOTYPE unless SITES specifically selected
            vqScoreFilteringType = performSiteSpecificVQScoreFiltering ? VQScoreFilteringType.SITES : VQScoreFilteringType.GENOTYPE;
        }

        if (convertFilteredGenotypesToNoCalls && vqScoreFilteringType != VQScoreFilteringType.GENOTYPE) {
            throw new UserException("The option '--convert-filtered-genotypes-to-no-calls' can ONLY be used if you are filtering at the " +
                    "Genotype level (you have set '--filter-set-info-table' and NOT set '--vqs-score-filter-by-site')");
        }

        // filter at a site level (but not necessarily use vqslod)
        if ((filterSetSiteTableName != null && filterSetName == null) || (filterSetSiteTableName == null && filterSetName != null)) {
            throw new UserException("--filter-set-name and --filter-set-site-table are both necessary for any filtering related operations");
        }
        if (!vqScoreFilteringType.equals(VQScoreFilteringType.NONE)) {
            //noinspection ConstantValue
            if (filterSetInfoTableName == null || filterSetSiteTableName == null || filterSetName == null) {
                throw new UserException(" --filter-set-info-table, --filter-set-name and --filter-set-site-table are all necessary for any Variant Quality" +
                        " filtering operations");
            }
        }

        Set<VCFHeaderLine> extraHeaderLines = new HashSet<>();
        if (!vqScoreFilteringType.equals(VQScoreFilteringType.NONE)) {
            if (isVQSR) {
                extraHeaderLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.AS_VQS_LOD_KEY));
                FilterSensitivityTools.validateFilteringCutoffs(truthSensitivitySNPThreshold, truthSensitivityINDELThreshold, vqsLodSNPThreshold, vqsLodINDELThreshold, tranchesTableName);
                Map<String, Map<Double, Double>> trancheMaps = FilterSensitivityTools.getTrancheMaps(filterSetName, tranchesTableName, projectID);

                if (vqsLodSNPThreshold != null) { // we already have vqslod thresholds directly
                    extraHeaderLines.add(FilterSensitivityTools.getVqsLodHeader(vqsLodSNPThreshold, GATKVCFConstants.SNP));
                    extraHeaderLines.add(FilterSensitivityTools.getVqsLodHeader(vqsLodINDELThreshold, GATKVCFConstants.INDEL));
                } else { // using sensitivity threshold inputs; need to convert these to vqslod thresholds
                    vqsLodSNPThreshold = FilterSensitivityTools.getVqslodThreshold(trancheMaps.get(GATKVCFConstants.SNP), truthSensitivitySNPThreshold, GATKVCFConstants.SNP);
                    vqsLodINDELThreshold = FilterSensitivityTools.getVqslodThreshold(trancheMaps.get(GATKVCFConstants.INDEL), truthSensitivityINDELThreshold, GATKVCFConstants.INDEL);
                    // set headers

                    if (vqScoreFilteringType.equals(VQScoreFilteringType.SITES)) {
                        extraHeaderLines.add(FilterSensitivityTools.getTruthSensitivityFilterHeader(truthSensitivitySNPThreshold, vqsLodSNPThreshold, GATKVCFConstants.SNP));
                        extraHeaderLines.add(FilterSensitivityTools.getTruthSensitivityFilterHeader(truthSensitivityINDELThreshold, vqsLodINDELThreshold, GATKVCFConstants.INDEL));
                    }
                    else if (vqScoreFilteringType.equals(VQScoreFilteringType.GENOTYPE)) {
                        extraHeaderLines.add(FilterSensitivityTools.getTruthSensitivityHeader(truthSensitivitySNPThreshold, vqsLodSNPThreshold, GATKVCFConstants.SNP));
                        extraHeaderLines.add(FilterSensitivityTools.getTruthSensitivityHeader(truthSensitivityINDELThreshold, vqsLodINDELThreshold, GATKVCFConstants.INDEL));
                    }
                }
            } else {
                extraHeaderLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.SCORE_KEY));
                extraHeaderLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.CALIBRATION_SENSITIVITY_KEY));
                if (truthSensitivitySNPThreshold == null) {
                    truthSensitivitySNPThreshold = FilterSensitivityTools.DEFAULT_TRUTH_SENSITIVITY_THRESHOLD_SNPS;
                }
                truthSensitivitySNPThreshold /= 100.0;
                logger.info("Passing all SNP variants with " + GATKVCFConstants.CALIBRATION_SENSITIVITY_KEY + " < " + truthSensitivitySNPThreshold);

                if (truthSensitivityINDELThreshold == null) {
                    truthSensitivityINDELThreshold = FilterSensitivityTools.DEFAULT_TRUTH_SENSITIVITY_THRESHOLD_INDELS;
                }
                truthSensitivityINDELThreshold /= 100.0;
                logger.info("Passing all INDEL variants with " + GATKVCFConstants.CALIBRATION_SENSITIVITY_KEY + " < " + truthSensitivityINDELThreshold);

                if (vqScoreFilteringType.equals(VQScoreFilteringType.SITES)) {
                    extraHeaderLines.add(new VCFFilterHeaderLine(GATKVCFConstants.CALIBRATION_SENSITIVITY_FAILURE_SNP,
                            "Site failed SNP model calibration sensitivity cutoff (" + truthSensitivitySNPThreshold.toString() + ")"));
                    extraHeaderLines.add(new VCFFilterHeaderLine(GATKVCFConstants.CALIBRATION_SENSITIVITY_FAILURE_INDEL,
                            "Site failed INDEL model calibration sensitivity cutoff (" + truthSensitivityINDELThreshold.toString() + ")"));
                }
                else if (vqScoreFilteringType.equals(VQScoreFilteringType.GENOTYPE)) {
                    extraHeaderLines.add(new VCFHeaderLine(GATKVCFConstants.CALIBRATION_SENSITIVITY_FAILURE_SNP,
                            "Sample Genotype FT filter value indicating that the genotyped allele failed SNP model calibration sensitivity cutoff (" + truthSensitivitySNPThreshold.toString() + ")"));
                    extraHeaderLines.add(new VCFHeaderLine(GATKVCFConstants.CALIBRATION_SENSITIVITY_FAILURE_INDEL,
                            "Sample Genotype FT filter value indicating that the genotyped allele failed INDEL model calibration sensitivity cutoff (" + truthSensitivityINDELThreshold.toString() + ")"));
                }
            }
        }

        if (vqScoreFilteringType.equals(VQScoreFilteringType.GENOTYPE)) {
            extraHeaderLines.add(new VCFFormatHeaderLine("FT", 1, VCFHeaderLineType.String, "Sample Genotype Filter Field"));
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

        reference = directlyAccessEngineReferenceDataSource();

        if (vetRangesExtractTableVersion == VetRangesExtractVersionEnum.V2) {
            extraHeaderLines.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_GT_KEY));
            extraHeaderLines.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_ID_KEY));
            extraHeaderLines.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.PHASE_SET_KEY));
        }

        header = generateVcfHeader(new HashSet<>(sampleIdToName.values()), reference.getSequenceDictionary(), extraHeaderLines);

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

        OptionalLong optionalMaximumAlternateAlleles =
                maximumAlternateAlleles != null && maximumAlternateAlleles > 0L ?
                        OptionalLong.of(maximumAlternateAlleles) : OptionalLong.empty();

        engine = !isVQSR ?
                new ExtractCohortVETSEngine(
                        projectID,
                        header,
                        annotationEngine,
                        reference,
                        sampleIdToName,
                        vetRangesFQDataSet,
                        fqRangesExtractVetTable,
                        vetRangesExtractTableVersion,
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
                        truthSensitivitySNPThreshold,
                        truthSensitivityINDELThreshold,
                        filterSetName,
                        emitPLs,
                        emitADs,
                        vqScoreFilteringType,
                        convertFilteredGenotypesToNoCalls,
                        optionalMaximumAlternateAlleles,
                        inferredReferenceState,
                        presortedAvroFiles,
                        this::apply)
                :
                new ExtractCohortEngine(
                        projectID,
                        header,
                        annotationEngine,
                        reference,
                        sampleIdToName,
                        vetRangesFQDataSet,
                        fqRangesExtractVetTable,
                        vetRangesExtractTableVersion,
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
                        filterSetName,
                        emitPLs,
                        emitADs,
                        vqScoreFilteringType,
                        convertFilteredGenotypesToNoCalls,
                        optionalMaximumAlternateAlleles,
                        inferredReferenceState,
                        presortedAvroFiles,
                        this::apply);
    }

    @Override
    // maybe think about creating a BigQuery Row walker?
    public void traverse() {
        progressMeter.setRecordsBetweenTimeChecks(100L);

        if (filterSetInfoTableName == null || filterSetInfoTableName.equals("")) {
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

        if (engine != null) {
            logger.info(String.format("***Processed %d total sites", engine.getTotalNumberOfSites()));
            logger.info(String.format("***Processed %d total variants", engine.getTotalNumberOfVariants()));
            logger.info(String.format("***Read API scanned %d bytes", engine.getTotalEstimatedBytesScanned()));
        }
    }
}
