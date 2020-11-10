package org.broadinstitute.hellbender.tools.variantdb.arrays;

import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import java.io.File;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ShortVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.tools.variantdb.CommonCode;
import org.broadinstitute.hellbender.tools.variantdb.SampleList;
import org.broadinstitute.hellbender.tools.variantdb.arrays.tables.ProbeInfo;
import org.broadinstitute.hellbender.tools.variantdb.arrays.tables.ProbeQcMetrics;
import org.broadinstitute.hellbender.tools.variantdb.nextgen.ExtractCohort;
import org.broadinstitute.hellbender.tools.walkers.annotator.Annotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.StandardAnnotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.AS_StandardAnnotation;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.util.*;


@CommandLineProgramProperties(
        summary = "(\"ExtractCohort\") - Filter and extract arrayvariants out of big query.",
        oneLineSummary = "Tool to extract variants out of big query for a subset of samples",
        programGroup = ShortVariantDiscoveryProgramGroup.class
)
@DocumentedFeature
public class ArrayExtractCohort extends GATKTool {
    private static final Logger logger = LogManager.getLogger(ExtractCohort.class);
    public static final int DEFAULT_LOCAL_SORT_MAX_RECORDS_IN_RAM = 1000000;
    private VariantContextWriter vcfWriter = null;
    private ArrayExtractCohortEngine engine;

    public enum QueryMode {
        LOCAL_SORT,
        QUERY
    }

    @Argument(
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName  = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            doc = "Output VCF file to which annotated variants should be written.",
            optional = false
    )
    private String outputVcfPathString = null;

    @Argument(
            fullName = "read-project-id",
            doc = "ID of the Google Cloud project to use (bill) when reading the microarray data tables",
            optional = true
    )
    private String readProjectID = null;

    @Argument(
            fullName = "cohort-sample-table",
            doc = "Fully qualified name of a bigquery table containing the sample_id and sample_name for the samples in the cohort you are extracting ",
            optional = true
    )
    private String sampleTableName = null;

   @Argument(
       fullName = "cohort-sample-file",
       doc = "CSV of sample_id,sample_name map in the cohort",
       optional = true
   )
    private File cohortSampleFile = null;

    @Argument(
            fullName = "probe-info-table",
            doc = "Fully qualified name of a bigquery table containing probe information",
            optional = true
    )
    private String probeTableName = null;

    @Argument(
        fullName = "probe-info-csv",
        doc = "Filepath to CSV export of probe-info table",
        optional = true
    )
    private String probeCsvExportFile = null;

    @Argument(
            fullName = "min-probe-id",
            doc = "When extracting data, only include probes with id >= this value",
            optional = true
    )
    private Integer minProbeId = null;

    @Argument(
        fullName = "max-probe-id",
        doc = "When extracting data, only include probes with id <= this value",
        optional = true
    )
    private Integer maxProbeId = null;

    @Argument(
        fullName = "qc-metrics-table",
        doc = "Fully qualified name of a bigquery table containing probe qc information",
        optional = true
    )
    private String qcMetricsTableName = null;

    @Argument(
            fullName = "cohort-extract-table",
            doc = "Fully qualified name of the table where the cohort data exists (already subsetted)",
            optional = false
    )
    private String cohortTable = null;

    @Argument(
            fullName = "gt-only",
            doc = "If true, only get the genotype info. Otherwise include NORMX, NORMY, BAF, and LRR",
            optional = true)
    private boolean gtDataOnly = false;

    @Argument(
        fullName = "remove-filtered-variants",
        doc = "If true, don't emit variants that are filtered out",
        optional = true)
    private boolean removeFilteredVariants = false;

    @Argument(
        fullName = "hwe-phred-scaled-threshold",
        doc = "Filter variants with HWE phred-scaled p-value greater than this value",
        optional = true)
    private float hwePvalThreshold = 60.0f;

    @Argument(
        fullName = "call-rate-threshold",
        doc = "Filter variants with a call rate less than this value",
        optional = true)
    private float callRateThreshold = 0.95f;
    
    @Argument(
        fullName = "filter-invariants",
        doc = "Filter variants with no called variants in the qc metrics",
        optional = true)
    private boolean filterInvariants = false;

    @Argument(
            fullName = "use-compressed-data",
            doc = "If true, use bit-packed fields for data",
            optional = true)
    private boolean useCompressedData = false;

    @Argument(
            fullName = "print-debug-information",
            doc = "If true, print extra debugging output",
            optional = true)
    private boolean printDebugInformation = false;

    @Argument(
            fullName = "local-sort-max-records-in-ram",
            doc = "When doing local sort, store at most this many records in memory at once",
            optional = true
    )
    private int localSortMaxRecordsInRam = DEFAULT_LOCAL_SORT_MAX_RECORDS_IN_RAM;

    // TODO remove before production
    @Argument(
            fullName = "use-legacy-gt-encoding",
            doc = "If the GT encodoing was AA, AB, BB",
            optional = true
    )
    private Boolean useLegacyGTEncoding = false;

    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    public boolean useVariantAnnotations() { return true; }

    @Override
    public List<Class<? extends Annotation>> getDefaultVariantAnnotationGroups() {
        return Arrays.asList(
                StandardAnnotation.class, AS_StandardAnnotation.class
        );
    }

    @Override
    protected void onStartup() {
        super.onStartup();

        if (minProbeId != null ^ maxProbeId != null) {
            throw new IllegalArgumentException("Either both or neither of min-probe-id and max-probe-id must be specified");
        }

        //TODO verify what we really need here
        final VariantAnnotatorEngine annotationEngine = new VariantAnnotatorEngine(makeVariantAnnotations(), null, Collections.emptyList(), false, false);

        vcfWriter = createVCFWriter(IOUtils.getPath(outputVcfPathString));

        SampleList sampleIdMap = new SampleList(sampleTableName, cohortSampleFile, printDebugInformation);
//        Map<Integer, String> sampleIdMap;
        VCFHeader header = CommonCode.generateRawArrayVcfHeader(new HashSet<>(sampleIdMap.getSampleNames()), reference.getSequenceDictionary());

        Map<Long, ProbeInfo> probeIdMap;
        if (probeCsvExportFile == null) {
            probeIdMap = ProbeInfo.getProbeIdMapWithStorageAPI(probeTableName, printDebugInformation, readProjectID);
        } else {
            probeIdMap = ProbeInfo.getProbeIdMapFromExport(probeCsvExportFile);
        }

        // if we have a qcMetrics table, augment the probeInfo map with that information
        Map<Long, ProbeQcMetrics> probeQcMetricsMap = null;
        if (qcMetricsTableName != null) {
            probeQcMetricsMap = ProbeQcMetrics.getProbeQcMetricsWithStorageAPI(qcMetricsTableName, readProjectID);
        }

        //ChromosomeEnum.setRefVersion(refVersion);

        engine = new ArrayExtractCohortEngine(
            readProjectID,
                vcfWriter,
                header,
                annotationEngine,
                reference,
                sampleIdMap,
                probeIdMap,
                probeQcMetricsMap,
                cohortTable,
                gtDataOnly,
                minProbeId,
                maxProbeId,
                localSortMaxRecordsInRam,
                useCompressedData,
                printDebugInformation,
                progressMeter,
                useLegacyGTEncoding,
                removeFilteredVariants,
                hwePvalThreshold,
                callRateThreshold,
                filterInvariants);
        vcfWriter.writeHeader(header);
    }

    @Override
    // maybe think about creating a BigQuery Row walker?
    public void traverse() {
        progressMeter.setRecordsBetweenTimeChecks(100L);
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
