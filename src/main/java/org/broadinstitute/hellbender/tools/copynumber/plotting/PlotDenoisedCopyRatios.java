package org.broadinstitute.hellbender.tools.copynumber.plotting;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.DenoiseReadCounts;
import org.broadinstitute.hellbender.tools.copynumber.coverage.copyratio.CopyRatioCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.CopyNumberStandardArgument;
import org.broadinstitute.hellbender.utils.R.RScriptExecutor;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.io.Resource;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * Plots the results of {@link DenoiseReadCounts}.
 *
 * <p>The order and representation of contigs in plots follows the contig ordering within the required reference sequence dictionary. </p>
 *
 * <h3>Examples</h3>
 *
 * <pre>
 * gatk-launch --javaOptions "-Xmx4g" PlotDenoisedCopyRatios \
 *   --standardizedCopyRatios tumor.standardizedCR.tsv \
 *   --denoisedCopyRatios tumor.denoisedCR.tsv \
 *   -SD ref_fasta.dict \
 *   --output output_dir \
 *   --outputPrefix tumor
 * </pre>
 *
 * <p>The --output parameter specifies a pre-existing directory.</p>
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Create plots of denoised copy ratios.",
        oneLineSummary = "Create plots of denoised copy ratios.",
        programGroup = CopyNumberProgramGroup.class
)
@DocumentedFeature
public final class PlotDenoisedCopyRatios extends CommandLineProgram {
    private static final String PLOT_DENOISED_COPY_RATIOS_R_SCRIPT = "PlotDenoisedCopyRatios.R";

    @Argument(
            doc = "Input file containing standardized copy-ratio profile (output of DenoiseReadCounts).",
            fullName = CopyNumberStandardArgument.STANDARDIZED_COPY_RATIOS_FILE_LONG_NAME,
            shortName = CopyNumberStandardArgument.STANDARDIZED_COPY_RATIOS_FILE_SHORT_NAME
    )
    private File inputStandardizedCopyRatiosFile;

    @Argument(
            doc = "Input file containing denoised copy-ratio profile (output of DenoiseReadCounts).",
            fullName = CopyNumberStandardArgument.DENOISED_COPY_RATIOS_FILE_LONG_NAME,
            shortName = CopyNumberStandardArgument.DENOISED_COPY_RATIOS_FILE_SHORT_NAME
    )
    private File inputDenoisedCopyRatiosFile;

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
            shortName = PlottingUtils.MINIMUM_CONTIG_LENGTH_SHORT_NAME,
            optional = true
    )
    private int minContigLength = PlottingUtils.DEFAULT_MINIMUM_CONTIG_LENGTH;

    @Argument(
            doc = "Prefix for output filenames.",
            fullName =  CopyNumberStandardArgument.OUTPUT_PREFIX_LONG_NAME,
            shortName = CopyNumberStandardArgument.OUTPUT_PREFIX_SHORT_NAME
    )
    private String outputPrefix;

    @Argument(
            doc = "Output directory.",
            fullName =  StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private String outputDir;

    @Override
    protected Object doWork() {
        checkRegularReadableUserFiles();

        logger.info("Reading and validating input files...");
        final CopyRatioCollection standardizedCopyRatios = new CopyRatioCollection(inputStandardizedCopyRatiosFile);
        final CopyRatioCollection denoisedCopyRatios = new CopyRatioCollection(inputDenoisedCopyRatiosFile);
        Utils.validateArg(standardizedCopyRatios.getIntervals().equals(denoisedCopyRatios.getIntervals()),
                "Intervals in input files must be identical.");

        //get sample name from input files (consistency check is performed)
        final String sampleName = getSampleName(standardizedCopyRatios, denoisedCopyRatios);

        //load contig names and lengths from the sequence dictionary into a LinkedHashMap
        final Map<String, Integer> contigLengthMap = PlottingUtils.getContigLengthMap(sequenceDictionaryFile, minContigLength, logger);

        //check that contigs in input files are present in sequence dictionary and that data points are valid given lengths
        PlottingUtils.validateContigs(contigLengthMap, standardizedCopyRatios, inputStandardizedCopyRatiosFile, logger);
        PlottingUtils.validateContigs(contigLengthMap, denoisedCopyRatios, inputDenoisedCopyRatiosFile, logger);

        logger.info("Generating plots...");
        final List<String> contigNames = new ArrayList<>(contigLengthMap.keySet());
        final List<Integer> contigLengths = new ArrayList<>(contigLengthMap.values());
        writeDenoisingPlots(sampleName, contigNames, contigLengths);

        return "SUCCESS";
    }

    private void checkRegularReadableUserFiles() {
        Utils.nonNull(outputPrefix);
        IOUtils.canReadFile(inputStandardizedCopyRatiosFile);
        IOUtils.canReadFile(inputDenoisedCopyRatiosFile);
        IOUtils.canReadFile(sequenceDictionaryFile);
        if (!new File(outputDir).exists()) {
            throw new UserException(String.format("Output directory %s does not exist.", outputDir));
        }
    }

    private String getSampleName(final CopyRatioCollection standardizedCopyRatios,
                                 final CopyRatioCollection denoisedCopyRatios) {
        final String standardizedSampleName = standardizedCopyRatios.getSampleName();
        final String denoisedSampleName = denoisedCopyRatios.getSampleName();
        Utils.validateArg(standardizedSampleName.equals(denoisedSampleName),"Sample names in input files must be identical.");
        return standardizedSampleName;
    }

    /**
     * @param sampleName Sample name derived from input files
     * @param contigNames List containing contig names
     * @param contigLengths List containing contig lengths (same order as contigNames)
     */
    private void writeDenoisingPlots(final String sampleName,
                                     final List<String> contigNames,
                                     final List<Integer> contigLengths) {
        final String contigNamesArg = contigNames.stream().collect(Collectors.joining(PlottingUtils.CONTIG_DELIMITER));                            //names separated by delimiter
        final String contigLengthsArg = contigLengths.stream().map(Object::toString).collect(Collectors.joining(PlottingUtils.CONTIG_DELIMITER));  //names separated by delimiter
        final String outputDirArg = PlottingUtils.addTrailingSlashIfNecessary(outputDir);

        final RScriptExecutor executor = new RScriptExecutor();

        //this runs the R statement "source("CNVPlottingLibrary.R")" before the main script runs
        executor.addScript(new Resource(PlottingUtils.CNV_PLOTTING_R_LIBRARY, PlotDenoisedCopyRatios.class));
        executor.addScript(new Resource(PLOT_DENOISED_COPY_RATIOS_R_SCRIPT, PlotDenoisedCopyRatios.class));
        //--args is needed for Rscript to recognize other arguments properly
        executor.addArgs("--args",
                "--sample_name=" + sampleName,
                "--standardized_copy_ratios_file=" + inputStandardizedCopyRatiosFile,
                "--denoised_copy_ratios_file=" + inputDenoisedCopyRatiosFile,
                "--contig_names=" + contigNamesArg,
                "--contig_lengths=" + contigLengthsArg,
                "--output_dir=" + outputDirArg,
                "--output_prefix=" + outputPrefix);
        executor.exec();
    }

}
