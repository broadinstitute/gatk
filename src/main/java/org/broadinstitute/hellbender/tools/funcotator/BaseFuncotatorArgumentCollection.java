package org.broadinstitute.hellbender.tools.funcotator;

import com.google.common.annotations.VisibleForTesting;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.Hidden;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation;

import java.io.File;
import java.io.Serializable;
import java.nio.file.Path;
import java.util.*;
import java.util.stream.Collectors;

public abstract class BaseFuncotatorArgumentCollection implements Serializable {
    private static final long serialVersionUID = 1L;

    /** String representing the b37 version of the homo sapiens reference. */
    protected static final String FuncotatorReferenceVersionB37 = "b37";
    /**
     * String representing the hg19 version of the homo sapiens reference.
     * This variable is necessary to resolve the differences between b37 and hg19 when
     * dealing with Homo Sapiens samples.
     */
    @VisibleForTesting
    public static final String FuncotatorReferenceVersionHg19 = "hg19";
    /** String representing the hg38 version of the homo sapiens reference. */
    @VisibleForTesting
    public static final String FuncotatorReferenceVersionHg38 = "hg38";

    @Argument(
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName  = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            doc = "Output file to which annotated variants should be written.")
    public File outputFile;

    @Argument(
            fullName =  FuncotatorArgumentDefinitions.REFERENCE_VERSION_LONG_NAME,
            doc = "The version of the Human Genome reference to use (e.g. hg19, hg38, etc.).  This will correspond to a sub-folder of each data source corresponding to that data source for the given reference."
    )
    public String referenceVersion;

    @Argument(
            fullName =  FuncotatorArgumentDefinitions.DATA_SOURCES_PATH_LONG_NAME,
            doc = "The path to a data source folder for Funcotator.  May be specified more than once to handle multiple data source folders."
    )
    public List<String> dataSourceDirectories;

    @Argument(
            fullName =  FuncotatorArgumentDefinitions.OUTPUT_FORMAT_LONG_NAME,
            doc = "The output file format.  Either VCF, MAF, or SEG.  Please note that MAF output for germline use case VCFs is unsupported.  SEG will generate two output files: a simple tsv and a gene list."
    )
    public FuncotatorArgumentDefinitions.OutputFormatType outputFormatType;

    //-----------------------------------------------------
    // Optional args:

    @Argument(
            fullName  = FuncotatorArgumentDefinitions.EXCLUSION_FIELDS_LONG_NAME,
            optional = true,
            doc = "Fields that should not be rendered in the final output.  Only exact name matches will be excluded."
    )
    public Set<String> excludedFields = new HashSet<>();

    @Advanced
    @Hidden
    @Argument(
            fullName = FuncotatorArgumentDefinitions.FORCE_B37_TO_HG19_REFERENCE_CONTIG_CONVERSION,
            optional = true,
            doc = "(Advanced / DO NOT USE*) If you select this flag, Funcotator will force a conversion of variant contig names from b37 to hg19.  *This option is useful in integration tests (written by devs) only."
    )
    public boolean forceB37ToHg19ContigNameConversion = false;

    @Argument(
            fullName  = FuncotatorArgumentDefinitions.TRANSCRIPT_SELECTION_MODE_LONG_NAME,
            optional = true,
            doc = "Method of detailed transcript selection.  This will select the transcript for detailed annotation (CANONICAL, ALL, or BEST_EFFECT)."
    )
    public TranscriptSelectionMode transcriptSelectionMode = FuncotatorArgumentDefinitions.TRANSCRIPT_SELECTION_MODE_DEFAULT_VALUE;

    @Argument(
            fullName  = FuncotatorArgumentDefinitions.TRANSCRIPT_LIST_LONG_NAME,
            optional = true,
            doc = "File to use as a list of transcripts (one transcript ID per line, version numbers are ignored) OR A set of transcript IDs to use for annotation to override selected transcript."
    )
    public Set<String> userTranscriptIdSet = new HashSet<>();

    @Argument(
            fullName  = FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_LONG_NAME,
            optional = true,
            doc = "Annotations to include in all annotated variants if the annotation is not specified in the data sources (in the format <ANNOTATION>:<VALUE>).  This will add the specified annotation to every annotated variant if it is not already present."
    )
    public List<String> annotationDefaults = new ArrayList<>();

    @Argument(
            fullName  = FuncotatorArgumentDefinitions.ANNOTATION_OVERRIDES_LONG_NAME,
            optional = true,
            doc = "Override values for annotations (in the format <ANNOTATION>:<VALUE>).  Replaces existing annotations of the given name with given values."
    )
    public List<String> annotationOverrides = new ArrayList<>();

    @Argument(
            fullName = FuncotatorArgumentDefinitions.LOOKAHEAD_CACHE_IN_BP_NAME,
            optional = true,
            minValue = 0,
            doc = "Number of base-pairs to cache when querying variants.  Can be overridden in individual data source configuration files."
    )
    public int lookaheadFeatureCachingInBp = FuncotatorArgumentDefinitions.LOOKAHEAD_CACHE_IN_BP_DEFAULT_VALUE;

    @Advanced
    @Argument(
            fullName = FuncotatorArgumentDefinitions.MIN_NUM_BASES_FOR_SEGMENT_FUNCOTATION,
            optional = true,
            doc = "The minimum number of bases for a variant to be annotated as a segment.  Recommended to be changed only for use with FuncotateSegments.  Defaults to " + FuncotatorUtils.DEFAULT_MIN_NUM_BASES_FOR_VALID_SEGMENT
    )
    public int minNumBasesForValidSegment = FuncotatorUtils.DEFAULT_MIN_NUM_BASES_FOR_VALID_SEGMENT;

    @Argument(
            fullName = FuncotatorArgumentDefinitions.CUSTOM_VARIANT_CLASS_ORDER_FILE,
            optional = true,
            doc = "TSV File containing custom Variant Classification severity map of the form: VARIANT_CLASSIFICATION\tSEV.  VARIANT_CLASSIFICATION must match one of the VariantClassification names (" + GencodeFuncotation.VariantClassification.ALL_VC_NAMES + ").  SEV is an unsigned integer, where lower is sorted first."
    )
    public GATKPath customVariantClassificationOrderFile = null;
}
