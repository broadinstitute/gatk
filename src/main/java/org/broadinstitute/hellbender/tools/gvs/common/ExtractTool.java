package org.broadinstitute.hellbender.tools.gvs.common;

import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.tools.walkers.annotator.Annotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.StandardAnnotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.AS_StandardAnnotation;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.File;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public abstract class ExtractTool extends GATKTool {
    private static final Logger logger = LogManager.getLogger(ExtractTool.class);
    public static final int DEFAULT_LOCAL_SORT_MAX_RECORDS_IN_RAM = 1000000;
    protected VariantContextWriter vcfWriter = null;
    protected VariantAnnotatorEngine annotationEngine;

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
            optional = true
    )
    protected String projectID = null;

    @Argument(
            fullName = "sample-table",
            doc = "Fully qualified name of a bigquery table containing a single column `sample` that describes the full list of samples to extract",
            optional = true,
            mutex={"sample-file"}
    )
    protected String sampleTableName = null;

    @Argument(
            fullName = "sample-file",
            doc = "Alternative to `sample-table`. Pass in a (sample_id,sample_name) CSV that describes the full list of samples to extract. No header",
            optional = true,
            mutex={"sample-table"}

    )
    protected File sampleFileName = null;

    @Argument(
            fullName = "print-debug-information",
            doc = "If true, print extra debugging output",
            optional = true)
    protected boolean printDebugInformation = false;

    @Argument(
            fullName = "local-sort-max-records-in-ram",
            doc = "When doing local sort, store at most this many records in memory at once",
            optional = true
    )
    protected int localSortMaxRecordsInRam = DEFAULT_LOCAL_SORT_MAX_RECORDS_IN_RAM;

    @Argument(
            fullName = "mode",
            doc = "Reference representation mode.  Valid values are PET or RANGES",
            optional = false
    )
    protected CommonCode.ModeEnum mode = CommonCode.ModeEnum.PET;

    @Argument(
            fullName = "ref-version",
            doc = "Remove this option!!!! only for ease of testing. Valid options are 37 or 38",
            optional = true
    )
    protected String refVersion = "37";

    @Argument(
        fullName = "min-location",
        doc = "When extracting data, only include locations >= this value",
        optional = true
    )
    protected Long minLocation = null;

    @Argument(
        fullName = "max-location",
        doc = "When extracting data, only include locations <= this value",
        optional = true
    )
    protected Long maxLocation = null;

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
