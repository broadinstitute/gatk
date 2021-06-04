package org.broadinstitute.hellbender.tools.gvs.extract;

import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.programgroups.ShortVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.gvs.common.CommonCode;
import org.broadinstitute.hellbender.tools.gvs.common.SampleList;
import org.broadinstitute.hellbender.tools.gvs.common.SchemaUtils;
import org.broadinstitute.hellbender.tools.gvs.common.ExtractTool;
import org.broadinstitute.hellbender.tools.gvs.common.FilterSensitivityTools;
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
            fullName = "disable-gnarly",
            doc = "Disable use of GnarlyGenotyper",
            optional = true
    )
    private boolean disableGnarlyGenotyper = true;

    @Argument(
            fullName = "vqslod-filter-genotypes",
            doc = "Should VQSLOD filtering be applied at the genotype level",
            optional = true
    )
    private boolean performGenotypeVQSLODFiltering = true;

    @Argument(
            fullName ="snps-truth-sensitivity-filter-level",
            doc = "The truth sensitivity level at which to start filtering SNPs",
            optional = true
    )
    private Double truthSensitivitySNPThreshold = null;

    @Argument(
            fullName = "indels-truth-sensitivity-filter-level",
            doc = "The truth sensitivity level at which to start filtering INDELs",
            optional = true
    )
    private Double truthSensitivityINDELThreshold = null;

    @Advanced
    @Argument(
            fullName = "snps-lod-score-cutoff",
            doc = "The VQSLOD score below which to start filtering SNPs",
            optional = true)
    private Double vqsLodSNPThreshold = null;

    @Advanced
    @Argument(
            fullName = "indels-lod-score-cutoff",
            doc = "The VQSLOD score below which to start filtering INDELs",
            optional = true)
    private Double vqsLodINDELThreshold = null;

    /**
    * If this flag is enabled, sites that have been marked as filtered (i.e. have anything other than `.` or `PASS`
    * in the FILTER field) will be excluded from the output.
    */
    @Argument(
            fullName="exclude-filtered",
            doc="Don't include filtered sites in the final jointVCF",
            optional=true)
    private boolean excludeFilteredSites = false;


    @Override
    protected void onStartup() {
        super.onStartup();

        if ( (filterSetInfoTableName != null || filterSetSiteTableName != null) && (filterSetName == null || filterSetName.equals(""))) {
            throw new UserException("--filter-set-name must be specified if any filtering related operations are requested");
        }

        Set<VCFHeaderLine> extraHeaderLines = new HashSet<>();
        if (filterSetInfoTableName != null) {
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

        extraHeaderLines.add(GATKVCFHeaderLines.getFilterLine(GATKVCFConstants.LOW_QUAL_FILTER_NAME));

        if (performGenotypeVQSLODFiltering) {
            extraHeaderLines.add(new VCFFormatHeaderLine("FT", 1, VCFHeaderLineType.String, "Genotype Filter Field"));
        }

        SampleList sampleList = new SampleList(sampleTableName, sampleFileName, projectID, printDebugInformation, "extract-cohort");
        Map<Long, String> sampleIdToName = sampleList.getMap();

        VCFHeader header = CommonCode.generateVcfHeader(new HashSet<>(sampleIdToName.values()), reference.getSequenceDictionary(), extraHeaderLines);

        final List<SimpleInterval> traversalIntervals = getTraversalIntervals();

        if (minLocation == null && maxLocation == null && hasUserSuppliedIntervals()) {
            final SimpleInterval firstInterval = traversalIntervals.get(0);
            final SimpleInterval lastInterval = traversalIntervals.get(traversalIntervals.size() - 1);

            minLocation = SchemaUtils.encodeLocation(firstInterval.getContig(), firstInterval.getStart());
            maxLocation = SchemaUtils.encodeLocation(lastInterval.getContig(), lastInterval.getEnd());
        } else if ((minLocation != null || maxLocation != null) && hasUserSuppliedIntervals()) {
            throw new UserException("min-location and max-location should not be used together with intervals (-L).");
        }

        // if there is a avro file, the BQ specific parameters are unnecessary,
        // but they all are required if there is no avro file
        if (cohortAvroFileName == null && (projectID == null || cohortTable == null)) {
            throw new UserException("Project id (--project-id) and cohort table (--cohort-extract-table) are required " +
                "if no avro file (--cohort-avro-file-name) is provided.");
        }

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
                mode,
                cohortTable,
                cohortAvroFileName,
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
                disableGnarlyGenotyper,
                performGenotypeVQSLODFiltering,
                excludeFilteredSites);

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
