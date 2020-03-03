package org.broadinstitute.hellbender.tools.copynumber.plotting;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.tools.copynumber.DenoiseReadCounts;
import org.broadinstitute.hellbender.tools.copynumber.ModelSegments;
import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberArgumentValidationUtils;
import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberStandardArgument;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.AllelicCountCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.CopyRatioCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.ModeledSegmentCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.AllelicCount;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CopyRatio;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.ModeledSegment;
import org.broadinstitute.hellbender.utils.R.RScriptExecutor;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.Resource;
import org.broadinstitute.hellbender.utils.reference.ReferenceUtils;

import java.io.File;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * Creates plots of denoised and segmented copy-ratio and minor-allele-fraction estimates.
 *
 * <h3>Inputs</h3>
 *
 * <ul>
 *     <li>
 *         Modeled-segments file from {@link ModelSegments}.
 *     </li>
 *     <li>
 *         (Optional) Denoised-copy-ratios file from  {@link DenoiseReadCounts}.
 *         If allelic counts are not provided, then this is required.
 *     </li>
 *     <li>
 *         (Optional) Allelic-counts file containing the counts at sites genotyped as heterozygous (.hets.tsv output of {@link ModelSegments}).
 *         If denoised copy ratios are not provided, then this is required.
 *     </li>
 *     <li>
 *         Sequence-dictionary file.
 *         This determines the order and representation of contigs in the plot.
 *     </li>
 *     <li>
 *         Output prefix.
 *         This is used as the basename for output files.
 *     </li>
 *     <li>
 *         Output directory.
 *         This will be created if it does not exist.
 *     </li>
 * </ul>
 *
 * <h3>Output</h3>
 *
 * <ul>
 *     <li>
 *         Modeled-segments-plot file.
 *         This shows the input denoised copy ratios and/or alternate-allele fractions as points, as well as box plots
 *         for the available posteriors in each segment.  The colors of the points alternate with the segmentation.
 *         Copy ratios are only plotted up to the maximum value specified by the argument {@code maximum-copy-ratio}.
 *         Point sizes can be specified by the arguments {@code point-size-copy-ratio} and {@code point-size-allele-fraction}.
 *     </li>
 * </ul>
 *
 * <h3>Usage examples</h3>
 *
 * <pre>
 *     gatk PlotModeledSegments \
 *          --denoised-copy-ratios tumor.denoisedCR.tsv \
 *          --allelic-counts tumor.hets.tsv \
 *          --segments tumor.modelFinal.seg \
 *          --sequence-dictionary contigs_to_plot.dict \
 *          --output-prefix tumor \
 *          -O output_dir
 * </pre>
 *
 * <pre>
 *     gatk PlotModeledSegments \
 *          --denoised-copy-ratios tumor.denoisedCR.tsv \
 *          --segments tumor.modelFinal.seg \
 *          --sequence-dictionary contigs_to_plot.dict \
 *          --output-prefix tumor \
 *          -O output_dir
 * </pre>
 *
 * <pre>
 *     gatk PlotModeledSegments \
 *          --allelic-counts normal.hets.tsv \
 *          --segments normal.modelFinal.seg \
 *          --sequence-dictionary contigs_to_plot.dict \
 *          --output-prefix normal \
 *          -O output_dir
 * </pre>
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Creates plots of denoised and segmented copy-ratio and minor-allele-fraction estimates",
        oneLineSummary = "Creates plots of denoised and segmented copy-ratio and minor-allele-fraction estimates",
        programGroup = CopyNumberProgramGroup.class
)
@DocumentedFeature
public final class PlotModeledSegments extends CommandLineProgram {
    private static final String PLOT_MODELED_SEGMENTS_R_SCRIPT = "PlotModeledSegments.R";
    private static final String MODELED_SEGMENTS_PLOT_FILE_SUFFIX = ".modeled.png";

    @Argument(
            doc = "Input file containing denoised copy ratios (output of DenoiseReadCounts).",
            fullName = CopyNumberStandardArgument.DENOISED_COPY_RATIOS_FILE_LONG_NAME,
            optional = true
    )
    private File inputDenoisedCopyRatiosFile;

    @Argument(
            doc = "Input file containing allelic counts at heterozygous sites (.hets.tsv output of ModelSegments).",
            fullName =  CopyNumberStandardArgument.ALLELIC_COUNTS_FILE_LONG_NAME,
            optional = true
    )
    private File inputAllelicCountsFile;

    @Argument(
            doc = "Input file containing modeled segments (output of ModelSegments).",
            fullName = CopyNumberStandardArgument.SEGMENTS_FILE_LONG_NAME
    )
    private File inputModeledSegmentsFile;

    @Argument(
            doc = PlottingUtils.SEQUENCE_DICTIONARY_DOC_STRING,
            fullName = StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME,
            shortName = StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME
    )
    private File inputSequenceDictionaryFile;

    @Argument(
            doc = PlottingUtils.MINIMUM_CONTIG_LENGTH_DOC_STRING,
            fullName =  PlottingUtils.MINIMUM_CONTIG_LENGTH_LONG_NAME,
            minValue = 0,
            optional = true
    )
    private int minContigLength = PlottingUtils.DEFAULT_MINIMUM_CONTIG_LENGTH;

    @Argument(
            doc = PlottingUtils.MAXIMUM_COPY_RATIO_DOC_STRING,
            fullName =  PlottingUtils.MAXIMUM_COPY_RATIO_LONG_NAME,
            minValue = 0,
            optional = true
    )
    private double maxCopyRatio = PlottingUtils.DEFAULT_MAXIMUM_COPY_RATIO;

    @Argument(
            doc = PlottingUtils.POINT_SIZE_COPY_RATIO_DOC_STRING,
            fullName =  PlottingUtils.POINT_SIZE_COPY_RATIO_LONG_NAME,
            minValue = 0,
            optional = true
    )
    private double pointSizeCopyRatio = PlottingUtils.DEFAULT_POINT_SIZE_COPY_RATIO;

    @Argument(
            doc = PlottingUtils.POINT_SIZE_ALLELE_FRACTION_DOC_STRING,
            fullName =  PlottingUtils.POINT_SIZE_ALLELE_FRACTION_LONG_NAME,
            minValue = 0,
            optional = true
    )
    private double pointSizeAlleleFraction = PlottingUtils.DEFAULT_POINT_SIZE_ALLELE_FRACTION;

    @Argument(
            doc = "Prefix for output filenames.",
            fullName =  CopyNumberStandardArgument.OUTPUT_PREFIX_LONG_NAME
    )
    private String outputPrefix;

    @Argument(
            doc = "Output directory.  This will be created if it does not exist.",
            fullName =  StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private File outputDir;

    //read input files
    private CopyRatioCollection denoisedCopyRatios;
    private AllelicCountCollection allelicCounts;
    private ModeledSegmentCollection modeledSegments;

    @Override
    protected Object doWork() {
        validateArguments();

        logger.info("Reading and validating input files...");
        denoisedCopyRatios = inputDenoisedCopyRatiosFile == null ? null : new CopyRatioCollection(inputDenoisedCopyRatiosFile);
        allelicCounts = inputAllelicCountsFile == null ? null : new AllelicCountCollection(inputAllelicCountsFile);
        modeledSegments = new ModeledSegmentCollection(inputModeledSegmentsFile);

        //get sample name from input files (consistency check is performed)
        final String sampleName = getSampleName();

        //perform a basic check on the total number of copy-ratio intervals/allele-fraction sites per contig
        //(we do not check the number of intervals/sites overlapping each segment, in order to be robust against future
        //changes to the segmentation that might assign intervals to more than one segment)
        validateNumPointsPerContig();

        //validate sequence dictionaries and load contig names and lengths into a LinkedHashMap
        final SAMSequenceDictionary sequenceDictionary = modeledSegments.getMetadata().getSequenceDictionary();
        final SAMSequenceDictionary sequenceDictionaryToPlot = ReferenceUtils.loadFastaDictionary(inputSequenceDictionaryFile);
        if (denoisedCopyRatios != null) {
            Utils.validateArg(CopyNumberArgumentValidationUtils.isSameDictionary(denoisedCopyRatios.getMetadata().getSequenceDictionary(), sequenceDictionary),
                    "Sequence dictionary contained in denoised-copy-ratios file does not match that contained in other input files.");
        }
        if (allelicCounts != null) {
            Utils.validateArg(CopyNumberArgumentValidationUtils.isSameDictionary(allelicCounts.getMetadata().getSequenceDictionary(), sequenceDictionary),
                    "Sequence dictionary contained in allelic-counts file does not match that contained in other input files.");
        }
        PlottingUtils.validateSequenceDictionarySubset(sequenceDictionary, sequenceDictionaryToPlot);
        final Map<String, Integer> contigLengthMap = PlottingUtils.getContigLengthMap(sequenceDictionaryToPlot, minContigLength, logger);

        //check that contigs in input files are present in sequence dictionary and that data points are valid given lengths
        PlottingUtils.validateContigs(contigLengthMap, denoisedCopyRatios, inputDenoisedCopyRatiosFile, logger);
        PlottingUtils.validateContigs(contigLengthMap, allelicCounts, inputAllelicCountsFile, logger);
        PlottingUtils.validateContigs(contigLengthMap, modeledSegments, inputModeledSegmentsFile, logger);

        final List<String> contigNames = new ArrayList<>(contigLengthMap.keySet());
        final List<Integer> contigLengths = new ArrayList<>(contigLengthMap.values());
        final File outputFile = new File(outputDir, outputPrefix + MODELED_SEGMENTS_PLOT_FILE_SUFFIX);

        logger.info(String.format("Writing plot to %s...", outputFile.getAbsolutePath()));
        writeModeledSegmentsPlot(sampleName, contigNames, contigLengths, outputFile);

        logger.info(String.format("%s complete.", getClass().getSimpleName()));

        return null;
    }

    private void validateArguments() {
        Utils.validateArg(!(inputDenoisedCopyRatiosFile == null && inputAllelicCountsFile == null),
                "Must provide at least a denoised-copy-ratios file or an allelic-counts file.");

        CopyNumberArgumentValidationUtils.validateInputs(
                inputDenoisedCopyRatiosFile,
                inputAllelicCountsFile,
                inputModeledSegmentsFile,
                inputSequenceDictionaryFile);
        Utils.nonEmpty(outputPrefix);
        CopyNumberArgumentValidationUtils.validateAndPrepareOutputDirectories(outputDir);
    }

    private String getSampleName() {
        if (inputDenoisedCopyRatiosFile != null) {
            Utils.validateArg(denoisedCopyRatios.getMetadata().equals(modeledSegments.getMetadata()),
                    "Metadata in input files do not all match.");
        }
        if (inputAllelicCountsFile != null) {
            Utils.validateArg(allelicCounts.getMetadata().equals(modeledSegments.getMetadata()),
                    "Metadata in input files do not all match.");
        }
        return modeledSegments.getMetadata().getSampleName();
    }

    private void validateNumPointsPerContig() {
        if (inputDenoisedCopyRatiosFile != null) {
            final Map<String, Integer> modeledSegmentsContigToNumPointsMap = modeledSegments.getRecords().stream()
                    .collect(Collectors.groupingBy(ModeledSegment::getContig, LinkedHashMap::new, Collectors.summingInt(ModeledSegment::getNumPointsCopyRatio)));
            final Map<String, Integer> denoisedCopyRatiosContigToNumPointsMap = denoisedCopyRatios.getRecords().stream()
                    .collect(Collectors.groupingBy(CopyRatio::getContig, LinkedHashMap::new, Collectors.summingInt(x -> 1)));
            Utils.validateArg(modeledSegmentsContigToNumPointsMap.keySet().stream()
                            .allMatch(c -> (!denoisedCopyRatiosContigToNumPointsMap.containsKey(c) && modeledSegmentsContigToNumPointsMap.get(c) == 0) ||
                                    (denoisedCopyRatiosContigToNumPointsMap.containsKey(c) && modeledSegmentsContigToNumPointsMap.get(c).equals(denoisedCopyRatiosContigToNumPointsMap.get(c)))),
                    "Number of denoised-copy-ratio points in input modeled-segments file is inconsistent with that in input denoised-copy-ratios file.");
        }
        if (inputAllelicCountsFile != null) {
            final Map<String, Integer> modeledSegmentsContigToNumPointsMap = modeledSegments.getRecords().stream()
                    .collect(Collectors.groupingBy(ModeledSegment::getContig, LinkedHashMap::new, Collectors.summingInt(ModeledSegment::getNumPointsAlleleFraction)));
            final Map<String, Integer> allelicCountsContigToNumPointsMap = allelicCounts.getRecords().stream()
                    .collect(Collectors.groupingBy(AllelicCount::getContig, LinkedHashMap::new, Collectors.summingInt(x -> 1)));
            Utils.validateArg(modeledSegmentsContigToNumPointsMap.keySet().stream()
                            .allMatch(c -> (!allelicCountsContigToNumPointsMap.containsKey(c) && modeledSegmentsContigToNumPointsMap.get(c) == 0) ||
                                    (allelicCountsContigToNumPointsMap.containsKey(c) && modeledSegmentsContigToNumPointsMap.get(c).equals(allelicCountsContigToNumPointsMap.get(c)))),
                    "Number of allelic-count points in input modeled-segments file is inconsistent with that in input heterozygous allelic-counts file.");
        }
    }

    /**
     * @param sampleName Sample name derived from input files
     * @param contigNames List containing contig names
     * @param contigLengths List containing contig lengths (same order as contigNames)
     */
    private void writeModeledSegmentsPlot(final String sampleName,
                                          final List<String> contigNames,
                                          final List<Integer> contigLengths,
                                          final File outputFile) {
        final String contigNamesArg = String.join(PlottingUtils.CONTIG_DELIMITER, contigNames); //names separated by delimiter
        final String contigLengthsArg = contigLengths.stream().map(Object::toString).collect(Collectors.joining(PlottingUtils.CONTIG_DELIMITER));  //lengths separated by delimiter
        final RScriptExecutor executor = new RScriptExecutor();

        //this runs the R statement "source("CNVPlottingLibrary.R")" before the main script runs
        executor.addScript(new Resource(PlottingUtils.CNV_PLOTTING_R_LIBRARY, PlotModeledSegments.class));
        executor.addScript(new Resource(PLOT_MODELED_SEGMENTS_R_SCRIPT, PlotModeledSegments.class));
        //--args is needed for Rscript to recognize other arguments properly
        executor.addArgs("--args",
                "--sample_name=" + sampleName,
                "--denoised_copy_ratios_file=" + (inputDenoisedCopyRatiosFile == null ? null : CopyNumberArgumentValidationUtils.getCanonicalPath(inputDenoisedCopyRatiosFile)),
                "--allelic_counts_file=" + (inputAllelicCountsFile == null ? null : CopyNumberArgumentValidationUtils.getCanonicalPath(inputAllelicCountsFile)),
                "--modeled_segments_file=" + CopyNumberArgumentValidationUtils.getCanonicalPath(inputModeledSegmentsFile),
                "--contig_names=" + contigNamesArg,
                "--contig_lengths=" + contigLengthsArg,
                "--maximum_copy_ratio=" + maxCopyRatio,
                "--point_size_copy_ratio=" + pointSizeCopyRatio,
                "--point_size_allele_fraction=" + pointSizeAlleleFraction,
                "--output_file=" + CopyNumberArgumentValidationUtils.getCanonicalPath(outputFile));
        executor.exec();
    }
}
