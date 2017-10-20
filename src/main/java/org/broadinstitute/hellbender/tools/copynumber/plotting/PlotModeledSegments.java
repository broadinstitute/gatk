package org.broadinstitute.hellbender.tools.copynumber.plotting;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.allelic.alleliccount.AllelicCount;
import org.broadinstitute.hellbender.tools.copynumber.allelic.alleliccount.AllelicCountCollection;
import org.broadinstitute.hellbender.tools.copynumber.coverage.copyratio.CopyRatio;
import org.broadinstitute.hellbender.tools.copynumber.coverage.copyratio.CopyRatioCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.CopyNumberStandardArgument;
import org.broadinstitute.hellbender.tools.copynumber.multidimensional.model.ModeledSegment;
import org.broadinstitute.hellbender.tools.copynumber.multidimensional.model.ModeledSegmentCollection;
import org.broadinstitute.hellbender.utils.R.RScriptExecutor;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.io.Resource;

import java.io.File;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * Plot segmented copy-ratio and minor-allele-fraction results.
 *
 * <p>The order and representation of contigs in plots follows the contig ordering within the required reference sequence dictionary. </p>
 *
 * <h3>Example</h3>
 *
 * <pre>
 * gatk-launch --javaOptions "-Xmx4g" PlotACNVResults \
 *   --hets tumor.hets.tsv \
 *   --denoisedNormalized tn_coverage.tn.tsv \
 *   --segments acnv_segments.seg \
 *   -SD ref_fasta.dict \
 *   --output folder_name \
 *   --outputPrefix basename
 * </pre>
 *
 * <p>The --output parameter specifies a directory.</p>
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Create plots of denoised and segmented copy-ratio and minor-allele-fraction estimates.",
        oneLineSummary = "Create plots of denoised and segmented copy-ratio and minor-allele-fraction estimates.",
        programGroup = CopyNumberProgramGroup.class
)
@DocumentedFeature
public final class PlotModeledSegments extends CommandLineProgram {
    private static final String PLOT_MODELED_SEGMENTS_R_SCRIPT = "PlotModeledSegments.R";

    @Argument(
            doc = "Input file containing denoised copy-ratio profile (output of DenoiseReadCounts).",
            fullName = CopyNumberStandardArgument.DENOISED_COPY_RATIOS_FILE_LONG_NAME,
            shortName = CopyNumberStandardArgument.DENOISED_COPY_RATIOS_FILE_SHORT_NAME,
            optional = true
    )
    private File inputDenoisedCopyRatiosFile;

    @Argument(
            doc = "Input file containing allelic counts at heterozygous sites (output of ModelSegments).",
            fullName =  CopyNumberStandardArgument.ALLELIC_COUNTS_FILE_LONG_NAME,
            shortName = CopyNumberStandardArgument.ALLELIC_COUNTS_FILE_SHORT_NAME,
            optional = true
    )
    private File inputAllelicCountsFile;

    @Argument(
            doc = "Input file containing modeled segments (output of ModelSegments).",
            fullName = CopyNumberStandardArgument.SEGMENTS_FILE_LONG_NAME,
            shortName = CopyNumberStandardArgument.SEGMENTS_FILE_SHORT_NAME
    )
    private File inputModeledSegmentsFile;

    @Argument(
            doc = "File containing the reference sequence dictionary (used to determine relative contig lengths). " +
                    "Contigs will be plotted in the order given. " +
                    "Contig names should not include the string \"" + PlottingUtils.CONTIG_DELIMITER + "\". " +
                    "The tool only considers contigs in the given dictionary for plotting, and " +
                    "data for contigs absent in the dictionary generate only a warning. In other words, you may " +
                    "modify a reference dictionary for use with this tool to include only contigs for which plotting is desired, " +
                    "and sort the contigs to the order in which the plots should display the contigs.",
            shortName = StandardArgumentDefinitions.SEQUENCE_DICTIONARY_SHORT_NAME
    )
    private File sequenceDictionaryFile;

    @Argument(
            doc = "Threshold length (in bp) for contigs to be plotted. " +
                    "Contigs with lengths less than this threshold will not be plotted. " +
                    "This can be used to filter out mitochondrial contigs, unlocalized contigs, etc.",
            fullName =  PlottingUtils.MINIMUM_CONTIG_LENGTH_LONG_NAME,
            shortName = PlottingUtils.MINIMUM_CONTIG_LENGTH_SHORT_NAME
    )
    private int minContigLength = PlottingUtils.DEFAULT_MINIMUM_CONTIG_LENGTH;

    @Argument(
            doc = "Prefix for output filenames.",
            fullName =  CopyNumberStandardArgument.OUTPUT_PREFIX_LONG_NAME,
            shortName = CopyNumberStandardArgument.OUTPUT_PREFIX_SHORT_NAME
    )
    protected String outputPrefix;

    @Argument(
            doc = "Output directory.",
            fullName =  StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    protected String outputDir;

    //read input files
    private CopyRatioCollection denoisedCopyRatios;
    private AllelicCountCollection allelicCounts;
    private ModeledSegmentCollection modeledSegments;

    @Override
    protected Object doWork() {
        checkRegularReadableUserFiles();

        logger.info("Reading and validating input files...");
        denoisedCopyRatios = inputDenoisedCopyRatiosFile == null ? null : new CopyRatioCollection(inputDenoisedCopyRatiosFile);
        allelicCounts = inputAllelicCountsFile == null ? null : new AllelicCountCollection(inputAllelicCountsFile);
        modeledSegments = new ModeledSegmentCollection(inputModeledSegmentsFile);

        //get sample name from input files (consistency check is performed)
        final String sampleName = getSampleName();

        //check that the number of copy-ratio/allele-fraction points in each segment are correct
        validateNumPoints();

        //load contig names and lengths from the sequence dictionary into a LinkedHashMap
        final Map<String, Integer> contigLengthMap = PlottingUtils.getContigLengthMap(sequenceDictionaryFile, minContigLength, logger);

        //check that contigs in input files are present in sequence dictionary and that data points are valid given lengths
        PlottingUtils.validateContigs(contigLengthMap, denoisedCopyRatios, inputDenoisedCopyRatiosFile, logger);
        PlottingUtils.validateContigs(contigLengthMap, allelicCounts, inputAllelicCountsFile, logger);
        PlottingUtils.validateContigs(contigLengthMap, modeledSegments, inputModeledSegmentsFile, logger);

        logger.info("Generating plot...");
        final List<String> contigNames = new ArrayList<>(contigLengthMap.keySet());
        final List<Integer> contigLengths = new ArrayList<>(contigLengthMap.values());
        writeModeledSegmentsPlot(sampleName, contigNames, contigLengths);
        return "SUCCESS";
    }

    private void checkRegularReadableUserFiles() {
        Utils.nonNull(outputPrefix);
        Utils.validateArg(!(inputDenoisedCopyRatiosFile == null && inputAllelicCountsFile == null),
                "Must provide at least a denoised copy-ratio profile file or an allelic-counts file.");
        if (inputDenoisedCopyRatiosFile != null) {
            IOUtils.canReadFile(inputDenoisedCopyRatiosFile);
        }
        if (inputAllelicCountsFile != null) {
            IOUtils.canReadFile(inputAllelicCountsFile);
        }
        IOUtils.canReadFile(inputModeledSegmentsFile);
        IOUtils.canReadFile(sequenceDictionaryFile);
        if (!new File(outputDir).exists()) {
            throw new UserException(String.format("Output directory %s does not exist.", outputDir));
        }
    }

    private String getSampleName() {
        if (inputDenoisedCopyRatiosFile != null) {
            Utils.validateArg(denoisedCopyRatios.getSampleName().equals(modeledSegments.getSampleName()),
                    "Sample names in input files do not all match.");
        }
        if (inputAllelicCountsFile != null) {
            Utils.validateArg(allelicCounts.getSampleName().equals(modeledSegments.getSampleName()),
                    "Sample names in input files do not all match.");
        }
        return modeledSegments.getSampleName();
    }

    private void validateNumPoints() {
        if (inputDenoisedCopyRatiosFile != null) {
            final Map<String, Integer> modeledSegmentsContigToNumPointsMap = modeledSegments.getRecords().stream()
                    .collect(Collectors.groupingBy(ModeledSegment::getContig, LinkedHashMap::new, Collectors.summingInt(ModeledSegment::getNumPointsCopyRatio)));
            final Map<String, Integer> denoisedCopyRatiosContigToNumPointsMap = denoisedCopyRatios.getRecords().stream()
                    .collect(Collectors.groupingBy(CopyRatio::getContig, LinkedHashMap::new, Collectors.summingInt(x -> 1)));
            Utils.validateArg(modeledSegmentsContigToNumPointsMap.keySet().stream()
                    .allMatch(c -> (modeledSegmentsContigToNumPointsMap.get(c) == 0 && !denoisedCopyRatiosContigToNumPointsMap.containsKey(c)) ||
                            (modeledSegmentsContigToNumPointsMap.get(c).equals(denoisedCopyRatiosContigToNumPointsMap.get(c)))),
                    "Number of denoised copy-ratio points in input modeled-segments file is inconsistent with that in input denoised copy-ratio file.");
        }
        if (inputAllelicCountsFile != null) {
            final Map<String, Integer> modeledSegmentsContigToNumPointsMap = modeledSegments.getRecords().stream()
                    .collect(Collectors.groupingBy(ModeledSegment::getContig, LinkedHashMap::new, Collectors.summingInt(ModeledSegment::getNumPointsAlleleFraction)));
            final Map<String, Integer> allelicCountsContigToNumPointsMap = allelicCounts.getRecords().stream()
                    .filter(ac -> modeledSegmentsContigToNumPointsMap.get(ac.getContig()) != 0)
                    .collect(Collectors.groupingBy(AllelicCount::getContig, LinkedHashMap::new, Collectors.summingInt(x -> 1)));
            Utils.validateArg(modeledSegmentsContigToNumPointsMap.keySet().stream()
                            .allMatch(c -> (modeledSegmentsContigToNumPointsMap.get(c) == 0 && !allelicCountsContigToNumPointsMap.containsKey(c)) ||
                                    (modeledSegmentsContigToNumPointsMap.get(c).equals(allelicCountsContigToNumPointsMap.get(c)))),
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
                                          final List<Integer> contigLengths) {
        final String contigNamesArg = contigNames.stream().collect(Collectors.joining(PlottingUtils.CONTIG_DELIMITER));                            //names separated by delimiter
        final String contigLengthsArg = contigLengths.stream().map(Object::toString).collect(Collectors.joining(PlottingUtils.CONTIG_DELIMITER));  //lengths separated by delimiter
        final String outputDirArg = PlottingUtils.addTrailingSlashIfNecessary(outputDir);
        final RScriptExecutor executor = new RScriptExecutor();

        //this runs the R statement "source("CNVPlottingLibrary.R")" before the main script runs
        executor.addScript(new Resource(PlottingUtils.CNV_PLOTTING_R_LIBRARY, PlotModeledSegments.class));
        executor.addScript(new Resource(PLOT_MODELED_SEGMENTS_R_SCRIPT, PlotModeledSegments.class));
        //--args is needed for Rscript to recognize other arguments properly
        executor.addArgs("--args",
                "--sample_name=" + sampleName,
                "--denoised_copy_ratios_file=" + inputDenoisedCopyRatiosFile,
                "--allelic_counts_file=" + inputAllelicCountsFile,
                "--modeled_segments_file=" + inputModeledSegmentsFile,
                "--contig_names=" + contigNamesArg,
                "--contig_lengths=" + contigLengthsArg,
                "--output_dir=" + outputDirArg,
                "--output_prefix=" + outputPrefix);
        executor.exec();
    }
}
