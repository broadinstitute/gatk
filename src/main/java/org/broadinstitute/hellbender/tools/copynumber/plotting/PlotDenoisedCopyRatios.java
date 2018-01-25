package org.broadinstitute.hellbender.tools.copynumber.plotting;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.DenoiseReadCounts;
import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberStandardArgument;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.CopyRatioCollection;
import org.broadinstitute.hellbender.utils.R.RScriptExecutor;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.io.Resource;
import org.broadinstitute.hellbender.utils.reference.ReferenceUtils;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * Creates plots of denoised copy ratios.  The tool also generates various denoising metrics.
 *
 * <h3>Inputs</h3>
 *
 * <ul>
 *     <li>
 *         Standardized-copy-ratios from {@link DenoiseReadCounts}.
 *     </li>
 *     <li>
 *         Denoised-copy-ratios from {@link DenoiseReadCounts}.
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
 *         This must be a pre-existing directory.
 *     </li>
 * </ul>
 *
 * <h3>Outputs</h3>
 *
 * <ul>
 *     <li>
 *         Denoised-plot files.
 *         Two versions of a plot showing both the standardized and denoised copy ratios are output;
 *         the first covers the entire range of the copy ratios, while the second is limited to copy ratios within [0, 4].
 *     </li>
 *      <li>
 *         Median-absolute-deviation files.
 *         These files contain the median absolute deviation (MAD) for both the standardized (.standardizedMAD.txt)
 *         and denoised (.denoisedMAD.txt) copy ratios, the change between the two (.deltaMAD.txt),
 *         and the change between the two scaled by the standardized MAD (.deltaScaledMAD.txt).
 *     </li>
 * </ul>
 *
 * <h3>Usage example</h3>
 *
 * <pre>
 *     gatk PlotDenoisedCopyRatios \
 *          --standardized-copy-ratios tumor.standardizedCR.tsv \
 *          --denoised-copy-ratios tumor.denoisedCR.tsv \
 *          --sequence-dictionary contigs_to_plot.dict \
 *          --output-prefix tumor \
 *          -O output_dir
 * </pre>
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Creates plots of denoised copy ratios",
        oneLineSummary = "Creates plots of denoised copy ratios",
        programGroup = CopyNumberProgramGroup.class
)
@DocumentedFeature
@BetaFeature
public final class PlotDenoisedCopyRatios extends CommandLineProgram {
    private static final String PLOT_DENOISED_COPY_RATIOS_R_SCRIPT = "PlotDenoisedCopyRatios.R";

    @Argument(
            doc = "Input file containing standardized copy ratios (output of DenoiseReadCounts).",
            fullName = CopyNumberStandardArgument.STANDARDIZED_COPY_RATIOS_FILE_LONG_NAME
    )
    private File inputStandardizedCopyRatiosFile;

    @Argument(
            doc = "Input file containing denoised copy ratios (output of DenoiseReadCounts).",
            fullName = CopyNumberStandardArgument.DENOISED_COPY_RATIOS_FILE_LONG_NAME
    )
    private File inputDenoisedCopyRatiosFile;

    @Argument(
            doc = PlottingUtils.SEQUENCE_DICTIONARY_DOC_STRING,
            fullName = StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME,
            shortName = StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME
    )
    private File inputSequenceDictionaryFile;

    @Argument(
            doc = PlottingUtils.MINIMUM_CONTIG_LENGTH_DOC_STRING,
            fullName =  PlottingUtils.MINIMUM_CONTIG_LENGTH_LONG_NAME,
            optional = true
    )
    private int minContigLength = PlottingUtils.DEFAULT_MINIMUM_CONTIG_LENGTH;

    @Argument(
            doc = "Prefix for output filenames.",
            fullName =  CopyNumberStandardArgument.OUTPUT_PREFIX_LONG_NAME
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

        //get sample name from input files (consistency check of metadata is performed)
        final String sampleName = getSampleName(standardizedCopyRatios, denoisedCopyRatios);

        //validate sequence dictionaries and load contig names and lengths into a LinkedHashMap
        final SAMSequenceDictionary sequenceDictionary = standardizedCopyRatios.getMetadata().getSequenceDictionary();
        final SAMSequenceDictionary sequenceDictionaryToPlot = ReferenceUtils.loadFastaDictionary(inputSequenceDictionaryFile);
        PlottingUtils.validateSequenceDictionarySubset(sequenceDictionary, sequenceDictionaryToPlot);
        final Map<String, Integer> contigLengthMap = PlottingUtils.getContigLengthMap(sequenceDictionaryToPlot, minContigLength, logger);

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
        IOUtils.canReadFile(inputSequenceDictionaryFile);
        if (!new File(outputDir).exists()) {
            throw new UserException(String.format("Output directory %s does not exist.", outputDir));
        }
    }

    private String getSampleName(final CopyRatioCollection standardizedCopyRatios,
                                 final CopyRatioCollection denoisedCopyRatios) {
        Utils.validateArg(standardizedCopyRatios.getMetadata().equals(denoisedCopyRatios.getMetadata()),
                "Metadata in input files must be identical.");
        return standardizedCopyRatios.getMetadata().getSampleName();
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
