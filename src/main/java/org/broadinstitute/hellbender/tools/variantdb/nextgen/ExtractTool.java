package org.broadinstitute.hellbender.tools.variantdb.nextgen;

import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.tools.variantdb.ChromosomeEnum;
import org.broadinstitute.hellbender.tools.variantdb.CommonCode;
import org.broadinstitute.hellbender.tools.variantdb.SchemaUtils;
import org.broadinstitute.hellbender.tools.variantdb.arrays.ExtractCohortBQ;
import org.broadinstitute.hellbender.tools.walkers.annotator.Annotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.StandardAnnotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.AS_StandardAnnotation;
import org.broadinstitute.hellbender.utils.bigquery.TableReference;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Set;

public abstract class ExtractTool extends GATKTool {
    private static final Logger logger = LogManager.getLogger(ExtractTool.class);
    public static final int DEFAULT_LOCAL_SORT_MAX_RECORDS_IN_RAM = 1000000;
    protected VariantContextWriter vcfWriter = null;
    protected VariantAnnotatorEngine annotationEngine;

    public enum QueryMode {
        LOCAL_SORT,
        QUERY
    }

    @Argument(
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            doc = "Output VCF file to which annotated variants should be written.",
            optional = false
    )
    protected String outputVcfPathString = null;

    @Argument(
            fullName = "project-id",
            doc = "ID of the Google Cloud project to use when executing queries",
            optional = false
    )
    protected String projectID = null;

    @Argument(
            fullName = "sample-table",
            doc = "Fully qualified name of a bigquery table containing a single column `sample` that describes the full list of samples to extract",
            optional = false
    )
    protected String sampleTableName = null;


    @Argument(
            fullName = "print-debug-information",
            doc = "If true, print extra debugging output",
            optional = true)
    protected boolean printDebugInformation = false;

    @Argument(
            fullName = "vqslog-SNP-threshold",
            doc = "The minimum value required for a SNP to pass.",
            optional = true)
    protected double vqsLodSNPThreshold = 99.95;

    @Argument(
            fullName = "vqslog-INDEL-threshold",
            doc = "The minimum value required for an INDEL to pass.",
            optional = true)
    protected double vqsLodINDELThreshold = 99.4;

    @Argument(
            fullName = "local-sort-max-records-in-ram",
            doc = "When doing local sort, store at most this many records in memory at once",
            optional = true
    )
    protected int localSortMaxRecordsInRam = DEFAULT_LOCAL_SORT_MAX_RECORDS_IN_RAM;

    @Argument(
            fullName = "mode",
            doc = "Source of genomic data. Valid options are one of ARRAYS, EXOMES, GENOMES",
            optional = false
    )
    protected CommonCode.ModeEnum mode = CommonCode.ModeEnum.EXOMES;

    @Argument(
            fullName = "query-mode",
            doc = "Source of genomic data. Valid options are one of GROUP_BY, LOCAL_SORT, QUERY",
            optional = false
    )
    protected ExtractTool.QueryMode queryMode = ExtractTool.QueryMode.QUERY;

    @Argument(
            fullName = "ref-version",
            doc = "Remove this option!!!! only for ease of testing. Valid options are 37 or 38",
            optional = true
    )
    protected String refVersion = "37";

    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    public boolean useVariantAnnotations() {
        return true;
    }

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
        annotationEngine = new VariantAnnotatorEngine(makeVariantAnnotations(), null, Collections.emptyList(), false, false);

        vcfWriter = createVCFWriter(IOUtils.getPath(outputVcfPathString));

        ChromosomeEnum.setRefVersion(refVersion);

    }
}
