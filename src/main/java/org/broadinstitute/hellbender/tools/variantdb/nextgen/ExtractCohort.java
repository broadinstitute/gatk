package org.broadinstitute.hellbender.tools.variantdb.nextgen;

import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.apache.avro.generic.GenericRecord;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.programgroups.ShortVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.variantdb.CommonCode;
import org.broadinstitute.hellbender.tools.variantdb.SampleList;
import org.broadinstitute.hellbender.tools.variantdb.SchemaUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.bigquery.StorageAPIAvroReader;
import org.broadinstitute.hellbender.utils.bigquery.TableReference;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

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

   @Argument(
            fullName = "variant-filter-table",
            doc = "Fully qualified name of the filtering table to use for cohort extraction",
            optional = true
    )
    private String filteringFQTableName = null;

    @Argument(
            fullName = "tranches-table",
            doc = "Fully qualified name of the tranches table to use for cohort extraction",
            optional = true
    )
    private String tranchesTableName = null;

    @Argument(
            fullName = "cohort-extract-table",
            doc = "Fully qualified name of the table where the cohort data exists (already subsetted)",
            optional = false
    )
    private String cohortTable = null;

    @Argument(
            fullName = "cohort-avro-file-name",
            doc = "Path of the cohort avro file",
            optional = true
    )
    private String cohortAvroFileName = null;

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
            fullName="snps-truth-sensitivity-filter-level",
            doc="The truth sensitivity level at which to start filtering SNPs",
            optional=true
    )
    private Double truthSensitivitySNPThreshold = null;

    @Argument(
            fullName="indels-truth-sensitivity-filter-level",
            doc="The truth sensitivity level at which to start filtering INDELs",
            optional=true
    )
    private Double truthSensitivityINDELThreshold = null;

    @Advanced
    @Argument(
            fullName="snps-lod-score-cutoff",
            doc="The VQSLOD score below which to start filtering SNPs",
            optional=true)
    private Double vqsLodSNPThreshold = null;

    @Advanced
    @Argument(
            fullName="indels-lod-score-cutoff",
            doc="The VQSLOD score below which to start filtering INDELs",
            optional=true)
    private Double vqsLodINDELThreshold = null;

    private static double DEFAULT_TRUTH_SENSITIVITY_THRESHOLD_SNPS = 99.7;
    private static double DEFAULT_TRUTH_SENSITIVITY_THRESHOLD_INDELS = 99.0;

    private Map<Double, Double> snpTranches = new HashMap<Double, Double>();
    private Map<Double, Double> indelTranches = new HashMap<Double, Double>();

    @Override
    protected void onStartup() {
        super.onStartup();

        SampleList sampleList = new SampleList(sampleTableName, sampleFileName, projectID, printDebugInformation);
        Set<String> sampleNames = new HashSet<>(sampleList.getSampleNames());

        Set<VCFHeaderLine> vqsrHeaderLines = new HashSet<>();
        if (filteringFQTableName != null) {
            validateFilteringCutoffs();
            getTranches();
            vqsrHeaderLines = getVqsLodThresholdsAndHeaders();
        }

        VCFHeader header = CommonCode.generateVcfHeader(sampleNames, reference.getSequenceDictionary(), vqsrHeaderLines);

        final List<SimpleInterval> traversalIntervals = getTraversalIntervals();

        if (minLocation == null && maxLocation == null && hasUserSuppliedIntervals()) {
            final SimpleInterval firstInterval = traversalIntervals.get(0);
            final SimpleInterval lastInterval = traversalIntervals.get(traversalIntervals.size() - 1);

            minLocation = SchemaUtils.encodeLocation(firstInterval.getContig(), firstInterval.getStart());
            maxLocation = SchemaUtils.encodeLocation(lastInterval.getContig(), lastInterval.getEnd());
        } else if ((minLocation != null || maxLocation != null) && hasUserSuppliedIntervals()) {
            throw new UserException("min-location and max-location should not be used together with intervals (-L).");
        }

        engine = new ExtractCohortEngine(
                projectID,
                vcfWriter,
                header,
                annotationEngine,
                reference,
                sampleNames,
                mode,
                cohortTable,
                cohortAvroFileName,
                traversalIntervals,
                minLocation,
                maxLocation,
                filteringFQTableName,
                localSortMaxRecordsInRam,
                printDebugInformation,
                vqsLodSNPThreshold,
                vqsLodINDELThreshold,
                progressMeter,
                queryMode,
                filterSetName,
                emitPLs);
        vcfWriter.writeHeader(header);
    }

    private void validateFilteringCutoffs() {
        if (truthSensitivitySNPThreshold != null ^ truthSensitivityINDELThreshold != null) {
            throw new UserException("If one of (--snps-truth-sensitivity-filter-level, --indels-truth-sensitivity-filter-level) is provided, both must be provided.");
        } else if (truthSensitivitySNPThreshold != null) {  // at this point, we know that if SNP is defined, INDEL is also defined
            // if the user specifies both truth sensitivity thresholds and lod cutoffs then throw a user error
            if (vqsLodSNPThreshold != null || vqsLodINDELThreshold != null) {
                throw new UserException("Arguments --[snps/indels]-truth-sensitivity-filter-level and --[snps/indels]-lod-score-cutoff are mutually exclusive. Please only specify one set of options.");
            }
        } else if (vqsLodSNPThreshold != null ^ vqsLodINDELThreshold != null) {
            throw new UserException("If one of (--snps-lod-score-cutoff, --indels-lod-score-cutoff) is provided, both must be provided.");
        } else if (vqsLodSNPThreshold == null) {  // at this point, we know that all vqsr threshold inputs are null, so use defaults
            // defaults if no values are given
            truthSensitivitySNPThreshold = DEFAULT_TRUTH_SENSITIVITY_THRESHOLD_SNPS;
            truthSensitivityINDELThreshold = DEFAULT_TRUTH_SENSITIVITY_THRESHOLD_INDELS;
        }

        if (vqsLodSNPThreshold == null) {
            // user must supply tranches table to look up vqslod score to use for cutoff
            if (tranchesTableName == null) {
                throw new UserException("Unless using lod score cutoffs (advanced), you must provide a tranches table using the argument --tranches-table.");
            }
        }
    }

    private void getTranches() {
        // get tranches from BigQuery
        final String restrictionWithFilterSetName = SchemaUtils.FILTER_SET_NAME + " = '" + filterSetName + "'";
        final TableReference tranchesTableRef = new TableReference(tranchesTableName, SchemaUtils.TRANCHE_FIELDS);
        final StorageAPIAvroReader tranchesTableAvroReader = new StorageAPIAvroReader(tranchesTableRef, restrictionWithFilterSetName, projectID);

        // format list
        for ( final GenericRecord queryRow : tranchesTableAvroReader ) {
            switch (queryRow.get(SchemaUtils.TRANCHE_MODEL).toString()) {
                case "SNP":
                    double targetSnpTruthSensitivity = Double.parseDouble(queryRow.get(SchemaUtils.TARGET_TRUTH_SENSITIVITY).toString());
                    double minSnpVqslod = Double.parseDouble(queryRow.get(SchemaUtils.MIN_VQSLOD).toString());
                    snpTranches.put(targetSnpTruthSensitivity, minSnpVqslod);
                    break;
                case "INDEL":
                    double targetIndelTruthSensitivity = Double.parseDouble(queryRow.get(SchemaUtils.TARGET_TRUTH_SENSITIVITY).toString());
                    double minIndelVqslod = Double.parseDouble(queryRow.get(SchemaUtils.MIN_VQSLOD).toString());
                    indelTranches.put(targetIndelTruthSensitivity, minIndelVqslod);
                    break;
            }
        }

        tranchesTableAvroReader.close();
    }

    private Set<VCFHeaderLine> getVqsLodThresholdsAndHeaders() {
        Set<VCFHeaderLine> vqsrHeaderLines = new HashSet<>();

        if (vqsLodSNPThreshold != null) { // user provided lod cutoffs
            vqsrHeaderLines.add(new VCFFilterHeaderLine(GATKVCFConstants.VQSR_FAILURE_SNP, "Site failed SNP model VQSLOD cutoff of " + vqsLodSNPThreshold.toString()));
            vqsrHeaderLines.add(new VCFFilterHeaderLine(GATKVCFConstants.VQSR_FAILURE_INDEL, "Site failed INDEL model VQSLOD cutoff of " + vqsLodINDELThreshold.toString()));
        } else { // user provided sensitivity or no cutoffs
            vqsLodSNPThreshold = getVqslodThreshold(snpTranches, truthSensitivitySNPThreshold, "SNP");
            vqsLodINDELThreshold = getVqslodThreshold(indelTranches, truthSensitivityINDELThreshold, "INDEL");

            vqsrHeaderLines.add(new VCFFilterHeaderLine(GATKVCFConstants.VQSR_FAILURE_SNP,
                    "Site failed SNP model sensitivity cutoff (" + truthSensitivitySNPThreshold.toString() + "), corresponding with VQSLOD cutoff of " + vqsLodSNPThreshold.toString()));
            vqsrHeaderLines.add(new VCFFilterHeaderLine(GATKVCFConstants.VQSR_FAILURE_INDEL,
                    "Site failed INDEL model sensitivity cutoff (" + truthSensitivityINDELThreshold.toString() + "), corresponding with VQSLOD cutoff of " + vqsLodINDELThreshold.toString()));
        }

        return vqsrHeaderLines;
    }

    private Double getVqslodThreshold(Map<Double, Double> tranches, Double truthSensitivityThreshold, String variantMode) {
        logger.info("Retrieving the min vqslod threshold for " + variantMode + "s and truth sensitivity of " + truthSensitivityThreshold);

        // We want to find the tranche with the smallest target_truth_sensitivity that is
        // equal to or greater than our truthSensitivityThreshold.
        // e.g. if truthSensitivitySNPThreshold is 99.8 and we have tranches with target_truth_sensitivities
        // of 99.5, 99.7, 99.9, and 100.0, we want the 99.9 sensitivity tranche.

        Double effectiveSensitivity = null;

        List<Double> sortedSensitivities = new ArrayList<>(tranches.keySet());
        Collections.sort(sortedSensitivities);  // sorts in ASCENDING order
        if (sortedSensitivities.contains(truthSensitivityThreshold)) {
            effectiveSensitivity = truthSensitivityThreshold;
        } else {
            // find the first sensitivity that's greater than our target truthSensitivityThreshold
            for ( Double trancheSensitivity : sortedSensitivities ) {
                if ( trancheSensitivity > truthSensitivityThreshold ) {
                    effectiveSensitivity = trancheSensitivity;
                    break;
                }
            }
        }

        if (effectiveSensitivity == null) {
            throw new UserException("No " + variantMode + " tranches found with target_truth_sensitivity >= " + truthSensitivityThreshold);
        }

        Double minVqslod = tranches.get(effectiveSensitivity);

        logger.info("Found " + variantMode + " tranche defined by sensitivity " + effectiveSensitivity + " and VQSLOD >= " + minVqslod + "; keeping all variants in this tranche.");
        logger.info("Passing all " + variantMode + " variants with VQSLOD >= " + minVqslod);

        return minVqslod;
    }

    @Override
    // maybe think about creating a BigQuery Row walker?
    public void traverse() {
        progressMeter.setRecordsBetweenTimeChecks(100L);

        if ( filteringFQTableName == null || filteringFQTableName.equals("") ) {
            logger.warn("--variant-filter-table is not specified, no filtering of cohort! ");
        }

        engine.traverse();
     }

    @Override
    protected void onShutdown() {
        super.onShutdown();

        if ( engine != null ) {
            logger.info(String.format("***Processed %d total sites", engine.getTotalNumberOfSites()));
            logger.info(String.format("***Processed %d total variants", engine.getTotalNumberOfVariants()));
        }

        // Close up our writer if we have to:
        if ( vcfWriter != null ) {
            vcfWriter.close();
        }
    }

}
