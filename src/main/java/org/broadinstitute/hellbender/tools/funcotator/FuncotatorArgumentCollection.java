package org.broadinstitute.hellbender.tools.funcotator;

import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.Hidden;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;

import java.io.File;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Arguments to be be used by the {@link Funcotator} {@link org.broadinstitute.hellbender.engine.GATKTool},
 * which are specific to {@link Funcotator}.
 * Created by jonn on 9/12/18.
 */
public class FuncotatorArgumentCollection implements Serializable {
    private static final long serialVersionUID = 1L;

    //-----------------------------------------------------
    // Required args:

    @Argument(
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName  = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            doc = "Output VCF file to which annotated variants should be written.")
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
            doc = "The output file format.  Either VCF or MAF.  Please note that MAF output for germline use case VCFs is unsupported."
    )
    public FuncotatorArgumentDefinitions.OutputFormatType outputFormatType;

    //-----------------------------------------------------
    // Optional args:

    @Argument(
            fullName = FuncotatorArgumentDefinitions.REMOVE_FILTERED_VARIANTS_LONG_NAME,
            optional = true,
            doc = "Ignore/drop variants that have been filtered in the input.  These variants will not appear in the output file."
    )
    public boolean removeFilteredVariants = false;

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
            fullName = FuncotatorArgumentDefinitions.FIVE_PRIME_FLANK_SIZE_NAME,
            optional = true,
            doc = "Variants within this many bases of the 5' end of a transcript (and not overlapping any part of the transcript itself) will be annotated as being in the 5' flanking region of that transcript"
    )
    public int fivePrimeFlankSize = FuncotatorArgumentDefinitions.FIVE_PRIME_FLANK_SIZE_DEFAULT_VALUE;

    @Argument(
            fullName = FuncotatorArgumentDefinitions.THREE_PRIME_FLANK_SIZE_NAME,
            optional = true,
            doc = "Variants within this many bases of the 3' end of a transcript (and not overlapping any part of the transcript itself) will be annotated as being in the 3' flanking region of that transcript"
    )
    public int threePrimeFlankSize = FuncotatorArgumentDefinitions.THREE_PRIME_FLANK_SIZE_DEFAULT_VALUE;

    @Argument(
            fullName = FuncotatorArgumentDefinitions.LOOKAHEAD_CACHE_IN_BP_NAME,
            optional = true,
            minValue = 0,
            doc = "Number of base-pairs to cache when querying variants."
    )
    public int lookaheadFeatureCachingInBp = FuncotatorArgumentDefinitions.LOOKAHEAD_CACHE_IN_BP_DEFAULT_VALUE;

    @Advanced
    @Hidden
    @Argument(
            fullName = FuncotatorArgumentDefinitions.FORCE_B37_TO_HG19_REFERENCE_CONTIG_CONVERSION,
            optional = true,
            doc = "(Advanced / DO NOT USE*) If you select this flag, Funcotator will force a conversion of variant contig names from b37 to hg19.  *This option is useful in integration tests (written by devs) only."
    )
    public boolean forceB37ToHg19ContigNameConversion = false;

    @Argument(
            fullName  = FuncotatorArgumentDefinitions.EXCLUSION_FIELDS_LONG_NAME,
            optional = true,
            doc = "Fields that should not be rendered in the final output.  Only exact name matches will be excluded."
    )
    public Set<String> excludedFields = new HashSet<>();
}
