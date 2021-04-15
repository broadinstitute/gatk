package org.broadinstitute.hellbender.tools.variantdb.nextgen;

import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
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
            mutex={"cohort-avro-file-name"},
            optional = true
    )
    private String cohortTable = null;

    @Argument(
            fullName = "cohort-avro-file-name",
            doc = "Path of the cohort avro file",
            mutex={"cohort-extract-table"},
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


    @Override
    protected void onStartup() {
        super.onStartup();

        Set<VCFHeaderLine> vqsrHeaderLines = new HashSet<>();
        if (filteringFQTableName != null) {
            FilterSensitivityTools.validateFilteringCutoffs(truthSensitivitySNPThreshold, truthSensitivityINDELThreshold, vqsLodSNPThreshold, vqsLodINDELThreshold, tranchesTableName);
            Map<String, Map<Double, Double>> trancheMaps = FilterSensitivityTools.getTrancheMaps(filterSetName, tranchesTableName, projectID);

            if (vqsLodSNPThreshold != null) { // we already have vqslod thresholds directly
                vqsrHeaderLines.add(FilterSensitivityTools.getVqsLodHeader(vqsLodSNPThreshold, GATKVCFConstants.SNP));
                vqsrHeaderLines.add(FilterSensitivityTools.getVqsLodHeader(vqsLodINDELThreshold, GATKVCFConstants.INDEL));
            } else { // using sensitivity threshold inputs; need to convert these to vqslod thresholds
                vqsLodSNPThreshold = FilterSensitivityTools.getVqslodThreshold(trancheMaps.get(GATKVCFConstants.SNP), truthSensitivitySNPThreshold, GATKVCFConstants.SNP);
                vqsLodINDELThreshold = FilterSensitivityTools.getVqslodThreshold(trancheMaps.get(GATKVCFConstants.INDEL), truthSensitivityINDELThreshold, GATKVCFConstants.INDEL);
                // set headers
                vqsrHeaderLines.add(FilterSensitivityTools.getTruthSensitivityHeader(truthSensitivitySNPThreshold, vqsLodSNPThreshold, GATKVCFConstants.SNP));
                vqsrHeaderLines.add(FilterSensitivityTools.getTruthSensitivityHeader(truthSensitivityINDELThreshold, vqsLodINDELThreshold, GATKVCFConstants.INDEL));
            }
        }

        SampleList sampleList = new SampleList(sampleTableName, sampleFileName, projectID, printDebugInformation);
        Set<String> sampleNames = new HashSet<>(sampleList.getSampleNames());

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

        // if there is a avro file, the BQ specific parameters are unnecessary,
        // but they all are required if there is no avro file
        if (cohortAvroFileName == null && (projectID == null || cohortTable == null)) {
            throw new UserException("a project id and cohort table are required if no avro file is provided");
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
