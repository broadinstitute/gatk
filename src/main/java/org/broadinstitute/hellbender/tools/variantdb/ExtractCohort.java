package org.broadinstitute.hellbender.tools.variantdb;

import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ShortVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.tools.walkers.annotator.Annotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.StandardAnnotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.AS_StandardAnnotation;
import org.broadinstitute.hellbender.utils.bigquery.TableReference;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.util.*;


@CommandLineProgramProperties(
        summary = "(\"ExtractCohort\") - Filter and extract variants out of big query.",
        oneLineSummary = "Tool to extract variants out of big query for a subset of samples",
        programGroup = ShortVariantDiscoveryProgramGroup.class
)
@DocumentedFeature
public class ExtractCohort extends GATKTool {
    private static final Logger logger = LogManager.getLogger(ExtractCohort.class);
    public static final int DEFAULT_LOCAL_SORT_MAX_RECORDS_IN_RAM = 1000000;
    private VariantContextWriter vcfWriter = null;
    private ExtractCohortEngine engine;
    public enum Mode {
        ARRAYS,
        EXOMES,
        GENOMES;
    };

    @Argument(
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName  = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            doc = "Output VCF file to which annotated variants should be written.",
            optional = false
    )
    private String outputVcfPathString = null;

    @Argument(
            fullName = "project-id",
            doc = "ID of the Google Cloud project to use when executing queries",
            optional = false
    )
    private String projectID = null;

    @Argument(
            fullName = "sample-table",
            doc = "Fully qualified name of a bigquery table containing a single column `sample` that describes the full list of samples to evoque",
            optional = true
    )
    private String sampleTableName = null;


    @Argument(
            fullName = "variant-filter-table",
            doc = "Fully qualified name of the filtering table to use for cohort extraction",
            optional = true
    )
    private String filteringFQTableName = null;

    @Argument(
            fullName = "cohort-extract-table",
            doc = "Fully qualified name of the table where the cohort data exists (already subsetted)",
            optional = false
    )
    private String cohortTable = null;

    @Argument(
            fullName = "print-debug-information",
            doc = "If true, print extra debugging output",
            optional = true)
    private boolean printDebugInformation = false;

    @Argument(
            fullName = "vqslog-SNP-threshold",
            doc = "The minimum value required for a SNP to pass.",
            optional = true)
    private double vqsLodSNPThreshold = 99.95;

    @Argument(
            fullName = "vqslog-INDEL-threshold",
            doc = "The minimum value required for an INDEL to pass.",
            optional = true)
    private double vqsLodINDELThreshold = 99.4;

    @Argument(
            fullName = "local-sort-max-records-in-ram",
            doc = "When doing local sort, store at most this many records in memory at once",
            optional = true
    )
    private int localSortMaxRecordsInRam = DEFAULT_LOCAL_SORT_MAX_RECORDS_IN_RAM;

    @Argument(
            fullName = "mode",
            doc = "Source of genomic data. Valid options are one of ARRAYS, EXOMES, GENOMES",
            optional = false
    )
    private Mode mode = Mode.EXOMES;

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

        //TODO verify what we really need here
        final VariantAnnotatorEngine annotationEngine = new VariantAnnotatorEngine(makeVariantAnnotations(), null, Collections.emptyList(), false, false);

        vcfWriter = createVCFWriter(IOUtils.getPath(outputVcfPathString));

        TableReference sampleTableRef = new TableReference(sampleTableName, SchemaUtils.SAMPLE_FIELDS);
        Set<String> sampleNames = ExtractCohortBQ.populateSampleNames(sampleTableRef, printDebugInformation);

        VCFHeader header = CommonCode.generateVcfHeader(sampleNames, reference.getSequenceDictionary());

        engine = new ExtractCohortEngine(
                projectID,
                vcfWriter,
                header,
                annotationEngine,
                reference,
                sampleNames,
                mode,
                cohortTable,
                filteringFQTableName,
                localSortMaxRecordsInRam,
                printDebugInformation,
                vqsLodSNPThreshold,
                vqsLodINDELThreshold,
                progressMeter);
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
