package org.broadinstitute.hellbender.tools.exome.plotting;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.*;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.*;
import org.broadinstitute.hellbender.utils.R.RScriptExecutor;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.io.Resource;
import org.broadinstitute.hellbender.utils.reference.ReferenceUtils;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Plots segmented coverage results.
 *
 * <p>The order of plotting follows the contig ordering within the required reference sequence dictionary. </p>
 *
 * <h3>Examples</h3>
 *
 * <p>The --output parameter specifies a directory.</p>
 *
 * <pre>
 * java -Xmx4g -jar $gatk_jar PlotSegmentedCopyRatio \
 *   --tangentNormalized tn_coverage.tn.tsv \
 *   --preTangentNormalized pre_tn_coverage.ptn.tsv \
 *   --segments called_segments.seg \
 *   -SD ref_fasta_dict.dict \
 *   --output output_dir \
 *   --outputPrefix $entity_id
 * </pre>
 */
@CommandLineProgramProperties(
        summary = "Create plots of denoised and segmented copy ratio.",
        oneLineSummary = "Create plots of denoised and segmented copy ratio",
        programGroup = CopyNumberProgramGroup.class
)
@DocumentedFeature
public final class PlotSegmentedCopyRatio extends CommandLineProgram {
    private static final String CNV_PLOTTING_R_LIBRARY = "CNVPlottingLibrary.R";
    private static final String COPY_RATIO_PLOTTING_R_SCRIPT = "CopyRatioPlotting.R";

    private static final String CONTIG_DELIMITER = "CONTIG_DELIMITER";  //used to delimit contig names and lengths passed to the R script
    private static final int DEFAULT_MINIMUM_CONTIG_LENGTH = 1000000;   //can be used to filter out mitochondrial contigs, unlocalized contigs, etc.

    //CLI arguments
    protected static final String OUTPUT_PREFIX_LONG_NAME = "outputPrefix";
    protected static final String OUTPUT_PREFIX_SHORT_NAME = "pre";

    protected static final String MINIMUM_CONTIG_LENGTH_LONG_NAME = "minimumContigLength";
    protected static final String MINIMUM_CONTIG_LENGTH_SHORT_NAME = "minContigLength";

    @Argument(
            doc = "File containing tangent-normalized coverage at targets, produced by NormalizeSomaticReadCounts (tn).",
            fullName =  ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME,
            optional = false
    )
    protected File tangentFile;

    @Argument(
            doc = "File containing pre-tangent-normalized coverage at targets, produced by NormalizeSomaticReadCounts (preTN).",
            fullName =  ExomeStandardArgumentDefinitions.PRE_TANGENT_NORMALIZED_COUNTS_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.PRE_TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME,
            optional = false
    )
    protected File preTangentFile;

    @Argument(
            doc = "File containing segmented regions of the genome, produced by PerformSegmentation.",
            fullName =  ExomeStandardArgumentDefinitions.SEGMENT_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.SEGMENT_FILE_SHORT_NAME,
            optional = false
    )
    protected File segmentsFile;

    @Argument(
            doc = "File containing the reference sequence dictionary (used to determine relative contig lengths). " +
                    "Contigs will be plotted in the order given. " +
                    "Contig names should not include \"" + CONTIG_DELIMITER + "\"." +
                    " Only data for contigs represented by the dictionary will be plotted. " +
                    "Data for contigs absent in the dictionary generate a warning.",
            shortName = StandardArgumentDefinitions.SEQUENCE_DICTIONARY_SHORT_NAME,
            optional = false
    )
    protected File sequenceDictionaryFile;

    @Argument(
            doc = "Threshold length (in bp) for contigs to be plotted. " +
                    "Contigs with lengths less than this threshold will not be plotted. " +
                    "This can be used to filter out mitochondrial contigs, unlocalized contigs, etc.",
            fullName =  MINIMUM_CONTIG_LENGTH_LONG_NAME,
            shortName = MINIMUM_CONTIG_LENGTH_SHORT_NAME,
            optional = false
    )
    protected int minContigLength = DEFAULT_MINIMUM_CONTIG_LENGTH;

    @Argument(
            doc = "Prefix for output filenames.",
            fullName =  OUTPUT_PREFIX_LONG_NAME,
            shortName = OUTPUT_PREFIX_SHORT_NAME,
            optional = false
    )
    protected String outputPrefix;

    @Argument(
            doc = "Output directory.",
            fullName =  StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            optional = false
    )
    protected String outputDir;

    @Argument(
            doc = "Set to true if input coverage data has had a log2 transform applied.",
            fullName =  ExomeStandardArgumentDefinitions.LOG2_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.LOG2_SHORT_NAME,
            optional = true
    )
    protected Boolean isLog2Input = true;

    @Override
    protected Object doWork() {
        checkRegularReadableUserFiles();

        //get sample name from input files (consistency check is performed)
        final String sampleName = getSampleName();

        //load contig names and lengths from the sequence dictionary into a LinkedHashMap
        final SAMSequenceDictionary sequenceDictionary = ReferenceUtils.loadFastaDictionary(sequenceDictionaryFile);
        Utils.validateArg(sequenceDictionary.getSequences().stream().map(SAMSequenceRecord::getSequenceName).noneMatch(n -> n.contains(CONTIG_DELIMITER)),
                String.format("Contig names cannot contain \"%s\".", CONTIG_DELIMITER));
        final Map<String, Integer> contigLengthMap = sequenceDictionary.getSequences().stream()
                .filter(s -> s.getSequenceLength() >= minContigLength)
                .collect(Collectors.toMap(SAMSequenceRecord::getSequenceName, SAMSequenceRecord::getSequenceLength,
                        (c, l) -> {
                            throw new IllegalArgumentException(String.format("Duplicate contig in sequence dictionary: %s", c));
                        },
                        LinkedHashMap::new));
        Utils.validateArg(contigLengthMap.size() > 0,
                "There must be at least one contig above the threshold length in the sequence dictionary.");
        logger.info("Contigs above length threshold: " + contigLengthMap.toString());

        //check that contigs in input files are present in sequence dictionary and that data points are valid given lengths
        validateContigs(contigLengthMap);

        //generate the plots
        final List<String> contigNames = new ArrayList<>(contigLengthMap.keySet());
        final List<Integer> contigLengths = new ArrayList<>(contigLengthMap.values());
        writeSegmentedCopyRatioPlots(sampleName, contigNames, contigLengths);

        return "SUCCESS";
    }

    private void checkRegularReadableUserFiles() {
        IOUtils.canReadFile(tangentFile);
        IOUtils.canReadFile(preTangentFile);
        IOUtils.canReadFile(segmentsFile);
        IOUtils.canReadFile(sequenceDictionaryFile);
        if (!new File(outputDir).exists()) {
            throw new UserException(String.format("Output directory %s does not exist.", outputDir));
        }
    }

    private String getSampleName() {
        final String tangentSampleName = ReadCountCollectionUtils.getSampleNameForCLIsFromReadCountsFile(tangentFile);
        final String preTangentSampleName = ReadCountCollectionUtils.getSampleNameForCLIsFromReadCountsFile(preTangentFile);
        final String segmentsSampleName = SegmentUtils.getSampleNameForCLIsFromSegmentFile(segmentsFile);
        Utils.validateArg(new HashSet<>(Arrays.asList(tangentSampleName, preTangentSampleName, segmentsSampleName)).size() == 1,
                "Sample names in input files do not all match.");
        return tangentSampleName;
    }

    private void validateContigs(final Map<String, Integer> contigLengthMap) {
        final Set<String> contigNames = contigLengthMap.keySet();

        //validate contig names and lengths in tangent file
        final ReadCountCollection tangent;
        try {
            tangent = ReadCountCollectionUtils.parse(tangentFile);
        } catch (final IOException e) {
            throw new UserException.CouldNotReadInputFile(tangentFile, e);
        }
        final Set<String> tangentContigNames = tangent.targets().stream().map(Target::getContig).collect(Collectors.toSet());
        if (!contigNames.containsAll(tangentContigNames)) {
            logger.warn("Contigs present in the tangent-normalized coverage file are missing from the sequence dictionary and will not be plotted.");
        }
        final Map<String, Integer> tangentContigMaxPositionMap = tangent.targets().stream().filter(t -> contigNames.contains(t.getContig()))
                .collect(Collectors.toMap(Target::getContig, Target::getEnd, Integer::max));
        tangentContigMaxPositionMap.keySet().forEach(c -> Utils.validateArg(tangentContigMaxPositionMap.get(c) <= contigLengthMap.get(c),
               "Position present in the tangent-normalized coverage file exceeds contig length in the sequence dictionary."));

        //validate contig names and lengths in pre-tangent file
        final ReadCountCollection preTangent;
        try {
            preTangent = ReadCountCollectionUtils.parse(preTangentFile);
        } catch (final IOException e) {
            throw new UserException.CouldNotReadInputFile(preTangentFile, e);
        }
        final Set<String> preTangentContigNames = preTangent.targets().stream().map(Target::getContig).collect(Collectors.toSet());
        if (!contigNames.containsAll(preTangentContigNames)) {
            logger.warn("Contigs present in the pre-tangent-normalized coverage file are missing from the sequence dictionary and will not be plotted.");
        }
        final Map<String, Integer> preTangentContigMaxPositionMap = preTangent.targets().stream().filter(t -> contigNames.contains(t.getContig()))
                .collect(Collectors.toMap(Target::getContig, Target::getEnd, Integer::max));
        preTangentContigMaxPositionMap.keySet().forEach(c -> Utils.validateArg(preTangentContigMaxPositionMap.get(c) <= contigLengthMap.get(c),
                "Position present in the pre-tangent-normalized coverage file exceeds contig length in the sequence dictionary."));

        //validate contig names and lengths in segments file
        final List<ModeledSegment> segments = SegmentUtils.readModeledSegmentsFromSegmentFile(segmentsFile);
        final Set<String> segmentsContigNames = segments.stream().map(ModeledSegment::getContig).collect(Collectors.toSet());
        if (!contigNames.containsAll(segmentsContigNames)) {
            logger.warn("Contigs present in the segments file are missing from the sequence dictionary and will not be plotted.");
        }
        final Map<String, Integer> segmentsContigMaxPositionMap = segments.stream().filter(s -> contigNames.contains(s.getContig()))
                .collect(Collectors.toMap(ModeledSegment::getContig, ModeledSegment::getEnd, Integer::max));
        segmentsContigMaxPositionMap.keySet().forEach(c -> Utils.validateArg(segmentsContigMaxPositionMap.get(c) <= contigLengthMap.get(c),
                "Position present in the segments file exceeds contig length in the sequence dictionary."));
    }

    /**
     * @param sampleName Sample name derived from input files
     * @param contigNames List containing contig names
     * @param contigLengths List containing contig lengths (same order as contigNames)
     */
    private void writeSegmentedCopyRatioPlots(final String sampleName,
                                              final List<String> contigNames,
                                              final List<Integer> contigLengths) {
        final String contigNamesArg = contigNames.stream().collect(Collectors.joining(CONTIG_DELIMITER));                            //names separated by delimiter
        final String contigLengthsArg = contigLengths.stream().map(Object::toString).collect(Collectors.joining(CONTIG_DELIMITER));  //names separated by delimiter
        final String outputDirArg = addTrailingSlashIfNecessary(outputDir);
        final String isLog2InputArg = isLog2Input ? "TRUE" : "FALSE";

        final RScriptExecutor executor = new RScriptExecutor();

        //this runs the R statement "source("CNVPlottingLibrary.R")" before the main script runs
        executor.addScript(new Resource(CNV_PLOTTING_R_LIBRARY, PlotSegmentedCopyRatio.class));
        executor.addScript(new Resource(COPY_RATIO_PLOTTING_R_SCRIPT, PlotSegmentedCopyRatio.class));
        //--args is needed for Rscript to recognize other arguments properly
        executor.addArgs("--args",
                "--sample_name=" + sampleName,
                "--tangent_file=" + tangentFile,
                "--pre_tangent_file=" + preTangentFile,
                "--segments_file=" + segmentsFile,
                "--contig_names=" + contigNamesArg,
                "--contig_lengths=" + contigLengthsArg,
                "--output_dir=" + outputDirArg,
                "--output_prefix=" + outputPrefix,
                "--is_log2_input=" + isLog2InputArg);
        executor.exec();
    }

    private static String addTrailingSlashIfNecessary(final String outputDir) {
        return outputDir.endsWith(File.separator) ? outputDir : outputDir + File.separator;
    }
}
