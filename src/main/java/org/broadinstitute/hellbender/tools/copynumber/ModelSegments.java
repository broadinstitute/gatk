package org.broadinstitute.hellbender.tools.copynumber;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableSet;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.OverlapDetector;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberArgumentValidationUtils;
import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberStandardArgument;
import org.broadinstitute.hellbender.tools.copynumber.arguments.SomaticGenotypingArgumentCollection;
import org.broadinstitute.hellbender.tools.copynumber.arguments.SomaticModelingArgumentCollection;
import org.broadinstitute.hellbender.tools.copynumber.arguments.SomaticSegmentationArgumentCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.AbstractLocatableCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.AbstractRecordCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.AllelicCountCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.CopyRatioCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.CopyRatioSegmentCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.LegacySegmentCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.ModeledSegmentCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.ParameterDecileCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.SimpleIntervalCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SimpleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CopyRatio;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CopyRatioSegment;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.LegacySegment;
import org.broadinstitute.hellbender.tools.copynumber.models.AlleleFractionModeller;
import org.broadinstitute.hellbender.tools.copynumber.models.AlleleFractionPrior;
import org.broadinstitute.hellbender.tools.copynumber.models.CopyRatioModeller;
import org.broadinstitute.hellbender.tools.copynumber.models.MultidimensionalModeller;
import org.broadinstitute.hellbender.tools.copynumber.segmentation.MultisampleMultidimensionalKernelSegmenter;
import org.broadinstitute.hellbender.tools.copynumber.utils.genotyping.NaiveHeterozygousPileupGenotypingUtils;
import org.broadinstitute.hellbender.tools.copynumber.utils.segmentation.KernelSegmenter;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Objects;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

/**
 * Models segmented copy ratios from denoised copy ratios and segmented minor-allele fractions from allelic counts.
 *
 * <p>
 *     Possible data inputs are: 1) denoised copy ratios for the case sample, 2) allelic counts for the case sample,
 *     and 3) allelic counts for a matched-normal sample.  All available inputs will be used to to perform
 *     segmentation and model inference.
 * </p>
 *
 * <h4>Genotyping step</h4>
 *
 * <p>
 *     If allelic counts are available, the first step in the inference process is to genotype heterozygous sites,
 *     as the allelic counts at these sites will subsequently be modeled to infer segmented minor-allele fraction.
 *     We perform a relatively simple and naive genotyping based on the allele counts (i.e., pileups), which is
 *     controlled by a small number of parameters ({@code minimum-total-allele-count},
 *     {@code genotyping-homozygous-log-ratio-threshold}, and {@code genotyping-homozygous-log-ratio-threshold}).
 *     If the matched normal is available, its allelic counts will be used to genotype the sites, and
 *     we will simply assume these genotypes are the same in the case sample.  (This can be critical, for example,
 *     for determining sites with loss of heterozygosity in high purity case samples; such sites will be genotyped as
 *     homozygous if the matched-normal sample is not available.)
 * </p>
 *
 * <h4>Segmentation step</h4>
 *
 * <p>
 *     Next, we segment, if available, the denoised copy ratios and the alternate-allele fractions at the
 *     genotyped heterozygous sites.  This is done using kernel segmentation (see {@link KernelSegmenter}).
 *     Various segmentation parameters control the sensitivity of the segmentation and should be selected
 *     appropriately for each analysis.  If a Picard interval-list file has been specified by the {@code segments}
 *     argument, the corresponding segmentation will be used instead and this step will be skipped.
 * </p>
 *
 * <p>
 *     If both copy ratios and allele fractions are available, we perform segmentation using a combined kernel
 *     that is sensitive to changes that occur not only in either of the two but also in both.  However, in this case,
 *     we simply discard all allele fractions at sites that lie outside of the available copy-ratio intervals
 *     (rather than imputing the missing copy-ratio data); these sites are filtered out during the genotyping step
 *     discussed above.  This can have implications for analyses involving the sex chromosomes;
 *     see comments in {@link CreateReadCountPanelOfNormals}.
 * </p>
 *
 * <h4>Modeling step</h4>
 *
 * <p>
 *     After segmentation is complete, we run Markov-chain Monte Carlo (MCMC) to determine posteriors for
 *     segmented models for the log2 copy ratio and the minor-allele fraction; see {@link CopyRatioModeller}
 *     and {@link AlleleFractionModeller}, respectively.  After the first run of MCMC is complete,
 *     smoothing of the segmented posteriors is performed by merging adjacent segments whose posterior
 *     credible intervals sufficiently overlap according to specified segmentation-smoothing parameters.
 *     Then, additional rounds of segmentation smoothing (with intermediate MCMC optionally performed in between rounds)
 *     are performed until convergence, at which point a final round of MCMC is performed.
 * </p>
 *
 * <h3>Inputs</h3>
 *
 * <ul>
 *     <li>
 *         (Optional) Denoised-copy-ratios file from {@link DenoiseReadCounts}.
 *         If allelic counts are not provided, then this is required.
 *     </li>
 *     <li>
 *         (Optional) Allelic-counts file from {@link CollectAllelicCounts}.
 *         If denoised copy ratios are not provided, then this is required.
 *     </li>
 *     <li>
 *         (Optional) Matched-normal allelic-counts file from {@link CollectAllelicCounts}.
 *         This can only be provided if allelic counts for the case sample are also provided.
 *     </li>
 *     <li>
 *         (Optional, Advanced) Picard interval-list file containing a multisample segmentation output by
 *         a previous run of {@link ModelSegments} in multisample mode.
 *         Segmentation step will not be performed.
 *         See description of multisample mode below.
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
 * <h3>Outputs</h3>
 *
 * <ul>
 *     <li>
 *         Modeled-segments .modelBegin.seg and .modelFinal.seg files.
 *         These are tab-separated values (TSV) files with a SAM-style header containing a read group sample name, a sequence dictionary,
 *         a row specifying the column headers contained in {@link ModeledSegmentCollection.ModeledSegmentTableColumn},
 *         and the corresponding entry rows.
 *         The initial result before segmentation smoothing is output to the .modelBegin.seg file
 *         and the final result after segmentation smoothing is output to the .modelFinal.seg file.
 *     </li>
 *     <li>
 *         Allele-fraction-model global-parameter files (.modelBegin.af.param and .modelFinal.af.param).
 *         These are tab-separated values (TSV) files with a SAM-style header containing a read group sample name,
 *         a row specifying the column headers contained in {@link ParameterDecileCollection.ParameterTableColumn},
 *         and the corresponding entry rows.
 *         The initial result before segmentation smoothing is output to the .modelBegin.af.param file
 *         and the final result after segmentation smoothing is output to the .modelFinal.af.param file.
 *     </li>
 *     <li>
 *         Copy-ratio-model global-parameter files (.modelBegin.cr.param and .modelFinal.cr.param).
 *         These are tab-separated values (TSV) files with a SAM-style header containing a read group sample name,
 *         a row specifying the column headers contained in {@link ParameterDecileCollection.ParameterTableColumn},
 *         and the corresponding entry rows.
 *         The initial result before segmentation smoothing is output to the .modelBegin.cr.param file
 *         and the final result after segmentation smoothing is output to the .modelFinal.cr.param file.
 *     </li>
 *     <li>
 *         Copy-ratio segments file (.cr.seg).
 *         This is a tab-separated values (TSV) file with a SAM-style header containing a read group sample name, a sequence dictionary,
 *         a row specifying the column headers contained in {@link CopyRatioSegmentCollection.CopyRatioSegmentTableColumn},
 *         and the corresponding entry rows.
 *         It contains the segments from the .modelFinal.seg file converted to a format suitable for input to {@link CallCopyRatioSegments}.
 *     </li>
 *     <li>
 *         CBS-formatted .cr.igv.seg and .af.igv.seg files (compatible with IGV).
 *         These are tab-separated values (TSV) files with CBS-format column headers
 *         (see <a href="http://software.broadinstitute.org/cancer/software/genepattern/file-formats-guide#CBS">
 *             http://software.broadinstitute.org/cancer/software/genepattern/file-formats-guide#CBS</a>)
 *         and the corresponding entry rows that can be plotted using IGV (see
 *         <a href="https://software.broadinstitute.org/software/igv/SEG">
 *             https://software.broadinstitute.org/software/igv/SEG</a>).
 *         The posterior medians of the log2 copy ratio and minor-allele fraction are given in the SEGMENT_MEAN
 *         columns in the .cr.igv.seg and .af.igv.seg files, respectively.
 *     </li>
 *     <li>
 *         (Optional) Allelic-counts file containing the counts at sites genotyped as heterozygous in the case sample (.hets.tsv).
 *         This is a tab-separated values (TSV) file with a SAM-style header containing a read group sample name, a sequence dictionary,
 *         a row specifying the column headers contained in {@link AllelicCountCollection.AllelicCountTableColumn},
 *         and the corresponding entry rows.
 *         This is only output if allelic counts are provided as input.
 *     </li>
 *     <li>
 *         (Optional) Allelic-counts file containing the counts at sites genotyped as heterozygous in the matched-normal sample (.hets.normal.tsv).
 *         This is a tab-separated values (TSV) file with a SAM-style header containing a read group sample name, a sequence dictionary,
 *         a row specifying the column headers contained in {@link AllelicCountCollection.AllelicCountTableColumn},
 *         and the corresponding entry rows.
 *         This is only output if matched-normal allelic counts are provided as input.
 *     </li>
 * </ul>
 *
 * <h3>Usage examples</h3>
 *
 * <pre>
 *     gatk ModelSegments \
 *          --denoised-copy-ratios tumor.denoisedCR.tsv \
 *          --allelic-counts tumor.allelicCounts.tsv \
 *          --normal-allelic-counts normal.allelicCounts.tsv \
 *          --output-prefix tumor \
 *          -O output_dir
 * </pre>
 *
 * <pre>
 *     gatk ModelSegments \
 *          --denoised-copy-ratios normal.denoisedCR.tsv \
 *          --allelic-counts normal.allelicCounts.tsv \
 *          --output-prefix normal \
 *          -O output_dir
 * </pre>
 *
 * <pre>
 *     gatk ModelSegments \
 *          --allelic-counts tumor.allelicCounts.tsv \
 *          --normal-allelic-counts normal.allelicCounts.tsv \
 *          --output-prefix tumor \
 *          -O output_dir
 * </pre>
 *
 * <pre>
 *     gatk ModelSegments \
 *          --denoised-copy-ratios normal.denoisedCR.tsv \
 *          --output-prefix normal \
 *          -O output_dir
 * </pre>
 *
 * <pre>
 *     gatk ModelSegments \
 *          --allelic-counts tumor.allelicCounts.tsv \
 *          --output-prefix tumor \
 *          -O output_dir
 * </pre>
 *
 * <h2>Multisample Mode (ADVANCED / EXPERIMENTAL)</h2>
 *
 * Multisample mode is activated when inputs for more than one case sample are specified via the
 * {@code denoised-copy-ratios} and/or {@code allelic-counts} arguments.  In this mode, {@link ModelSegments}
 * finds common segments across multiple case samples using denoised copy ratios and allelic counts.
 * This segmentation can be used as input to subsequent, individual runs of {@link ModelSegments}
 * on each of the case samples.
 *
 * <p>
 *     Possible data inputs are: 1) denoised copy ratios for the case samples, 2) allelic counts for the case samples,
 *     and 3) allelic counts for a matched-normal sample.  All available inputs will be used to to perform
 *     segmentation.
 * </p>
 *
 * <p>
 *     As in single-sample mode, the first step is to genotype heterozygous sites.  If allelic counts from
 *     a matched normal are available, heterozygous sites in all samples are defined using the normal, as usual;
 *     otherwise, each case sample is genotyped individually.  The intersection of heterozygous sites across
 *     case samples after filtering on total count and overlap with available copy-ratio intervals is then used for
 *     segmentation.  (When the matched normal is available, this intersection is inconsequential if
 *     {@code minimum-total-allele-count-case} is appropriately set to zero.  As in single-sample mode,
 *     determining sites with loss of heterozygosity in high purity case samples will be difficult if the
 *     matched normal is not available.)
 * </p>
 *
 * <p>
 *     Next, we jointly segment, if available, the denoised copy ratios and the alternate-allele fractions at the
 *     genotyped heterozygous sites across all case samples (which can include the matched normal, if desired).
 *     The same caveats discussed above also apply to segmentation in multisample mode; as in single-sample mode,
 *     segmentation parameters should be tuned appropriately for each analysis.  Moreover, parameters may also need
 *     to be tuned as a function of the number of samples used.
 * </p>
 *
 * <p>
 *     The final output of multisample mode is a Picard interval-list file specifying the joint segmentation.
 *     This can be provided to subsequent, individual runs of {@link ModelSegments} via the {@code segments} argument
 *     to perform modeling for each of the case samples; the segmentation step will be skipped in these runs.
 * </p>
 *
 * <p>
 *     Note that the genotyping step will be repeated in these runs, so filters identical to those used in the
 *     multisample-mode run should be used.  If allelic counts from a matched normal are available, the resulting set of
 *     heterozygous sites used for modeling in these runs should then be identical to that used for segmentation in
 *     the multisample-mode run if {@code minimum-total-allele-count-case} is appropriately set to zero for all runs.
 *     (However, if this is set to a non-zero value or if no matched normal is available, because the intersection
 *     of heterozygous sites across samples performed in the multisample-mode run cannot be replicated in the
 *     single-sample mode runs, this may ultimately yield sets of sites for modeling that differ across samples.)
 * </p>
 *
 * <p>
 *     See below for usage examples illustrating a run in multisample mode followed by multiple runs in single-sample mode.
 * </p>
 *
 * <h3>Inputs</h3>
 *
 * <ul>
 *     <li>
 *         (Optional) List of more than one denoised-copy-ratios files from {@link DenoiseReadCounts}.
 *         If a list of allelic counts is not provided, then this is required.
 *         If a list of allelic counts is also provided, then sample order must match across both lists.
 *     </li>
 *     <li>
 *         (Optional) List of more than one allelic-counts files from {@link CollectAllelicCounts}.
 *         If a list of denoised copy ratios is not provided, then this is required.
 *         If a list of denoised copy ratios is also provided, then sample order must match across both lists.
 *     </li>
 *     <li>
 *         (Optional) Matched-normal allelic-counts file from {@link CollectAllelicCounts}.
 *         This can only be provided if a list of allelic counts for the case samples are also provided.
 *     </li>
 * </ul>
 *
 * <h3>Outputs</h3>
 *
 * <ul>
 *     <li>
 *         Multisample-segments .interval_list Picard interval-list file.
 *         This segmentation can be used as input to subsequent, individual runs of {@link ModelSegments} on each of
 *         the case samples.
 *     </li>
 * </ul>
 *
 * <h3>Usage examples</h3>
 *
 * <h4>Multisample-mode run (outputs {@code multisample-segmentation.interval_list} containing the multisample segmentation) </h4>
 *
 * <pre>
 *     gatk ModelSegments \
 *          --denoised-copy-ratios normal.denoisedCR.tsv \
 *          --denoised-copy-ratios tumor-1.denoisedCR.tsv \
 *          ...
 *          --denoised-copy-ratios tumor-N.denoisedCR.tsv \
 *          --allelic-counts normal.allelicCounts.tsv \
 *          --allelic-counts tumor-1.allelicCounts.tsv \
 *          ...
 *          --allelic-counts tumor-N.allelicCounts.tsv \
 *          --normal-allelic-counts normal.allelicCounts.tsv \
 *          --output-prefix multisample-segmentation
 *          -O output_dir
 * </pre>
 *
 * <h4>Single-sample mode runs (segmentation is taken from {@code multisample-segmentation.interval_list},
 * so the segmentation step is skipped in each run)</h4>
 *
 * <pre>
 *     gatk ModelSegments \
 *          --segments multisample-segmentation.interval_list
 *          --denoised-copy-ratios normal.denoisedCR.tsv \
 *          --allelic-counts normal.allelicCounts.tsv \
 *          --normal-allelic-counts normal.allelicCounts.tsv \
 *          --output-prefix normal
 *          -O output_dir
 * </pre>
 *
 * <pre>
 *     gatk ModelSegments \
 *          --segments multisample-segmentation.interval_list
 *          --denoised-copy-ratios tumor-1.denoisedCR.tsv \
 *          --allelic-counts tumor-1.allelicCounts.tsv \
 *          --normal-allelic-counts normal.allelicCounts.tsv \
 *          --output-prefix tumor-1
 *          -O output_dir
 * </pre>
 *
 * <pre>
 *     ...
 * </pre>
 *
 * <pre>
 *     gatk ModelSegments \
 *          --segments multisample-segmentation.interval_list
 *          --denoised-copy-ratios tumor-N.denoisedCR.tsv \
 *          --allelic-counts tumor-N.allelicCounts.tsv \
 *          --normal-allelic-counts normal.allelicCounts.tsv \
 *          --output-prefix tumor-N
 *          -O output_dir
 * </pre>
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Models segmented copy ratios from denoised copy ratios and segmented minor-allele fractions from allelic counts; " +
                "if multiple samples are specified, finds a joint segmentation that can be used in subsequent runs to perform modeling of each sample",
        oneLineSummary = "Models segmented copy ratios from denoised copy ratios and segmented minor-allele fractions from allelic counts",
        programGroup = CopyNumberProgramGroup.class
)
@DocumentedFeature
public final class ModelSegments extends CommandLineProgram {
    public enum RunMode {
        MULTIPLE_SAMPLE, SINGLE_SAMPLE
    }

    private enum DataMode {
        COPY_RATIO_ONLY, ALLELE_FRACTION_ONLY, COPY_RATIO_AND_ALLELE_FRACTION
    }

    //filename tags for output
    public static final String HET_ALLELIC_COUNTS_FILE_SUFFIX = ".hets.tsv";
    public static final String NORMAL_HET_ALLELIC_COUNTS_FILE_SUFFIX = ".hets.normal.tsv";
    public static final String SEGMENTS_FILE_SUFFIX = ".seg";
    public static final String PICARD_INTERVAL_LIST_FILE_SUFFIX = ".interval_list";
    public static final String BEGIN_FIT_FILE_TAG = ".modelBegin";
    public static final String FINAL_FIT_FILE_TAG = ".modelFinal";
    public static final String COPY_RATIO_MODEL_PARAMETER_FILE_SUFFIX = ".cr.param";
    public static final String ALLELE_FRACTION_MODEL_PARAMETER_FILE_SUFFIX = ".af.param";
    public static final String COPY_RATIO_SEGMENTS_FOR_CALLER_FILE_SUFFIX = ".cr" + SEGMENTS_FILE_SUFFIX;
    public static final String COPY_RATIO_LEGACY_SEGMENTS_FILE_SUFFIX = ".cr.igv" + SEGMENTS_FILE_SUFFIX;
    public static final String ALLELE_FRACTION_LEGACY_SEGMENTS_FILE_SUFFIX = ".af.igv" + SEGMENTS_FILE_SUFFIX;

    @Argument(
            doc = "Input files containing denoised copy ratios (output of DenoiseReadCounts).  " +
                    "If multiple samples are specified, multisample kernel segmentation will be performed but modeling will be skipped; " +
                    "sample order must match that of input allelic-counts files.",
            fullName = CopyNumberStandardArgument.DENOISED_COPY_RATIOS_FILE_LONG_NAME,
            optional = true,
            minElements = 1
    )
    private List<File> inputDenoisedCopyRatiosFiles = null;

    @Argument(
            doc = "Input files containing allelic counts (output of CollectAllelicCounts).  " +
                    "If multiple samples are specified, multisample kernel segmentation will be performed but modeling will be skipped; " +
                    "sample order must match that of input denoised-copy-ratios files.",
            fullName = CopyNumberStandardArgument.ALLELIC_COUNTS_FILE_LONG_NAME,
            optional = true,
            minElements = 1
    )
    private List<File> inputAllelicCountsFiles = null;

    @Argument(
            doc = "Input file containing allelic counts for a matched normal (output of CollectAllelicCounts).  " +
                    "If specified, these allelic counts will be used to perform genotyping but will not be used for multisample kernel segmentation; " +
                    "if the latter is desired, additionally specify this file as one of the arguments to --allelic-counts.",
            fullName = CopyNumberStandardArgument.NORMAL_ALLELIC_COUNTS_FILE_LONG_NAME,
            optional = true
    )
    private File inputNormalAllelicCountsFile = null;

    @Advanced
    @Argument(
            doc = "Input Picard interval-list file specifying segments.  " +
                    "If provided, kernel segmentation will be skipped.",
            fullName = CopyNumberStandardArgument.SEGMENTS_FILE_LONG_NAME,
            optional = true
    )
    private File inputSegmentsFile = null;

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

    @ArgumentCollection
    private SomaticGenotypingArgumentCollection genotypingArguments = new SomaticGenotypingArgumentCollection();

    @ArgumentCollection
    private SomaticSegmentationArgumentCollection segmentationArguments = new SomaticSegmentationArgumentCollection();

    @ArgumentCollection
    private SomaticModelingArgumentCollection modelingArguments = new SomaticModelingArgumentCollection();

    private RunMode runMode;
    private DataMode dataMode;

    private void logHeapUsage(final String phase) {
        final int mb = 1024 * 1024;
        final Runtime runtime = Runtime.getRuntime();
        logger.info("Used memory (MB) after " + phase + ": " + (runtime.totalMemory() - runtime.freeMemory()) / mb);
    }

    @Override
    protected Object doWork() {
        logHeapUsage("initializing engine");

        setModesAndValidateArguments();

        //read input files, validate data, and perform genotyping
        final ModelSegmentsData modelSegmentsData = new ModelSegmentsData();

        if (runMode == RunMode.MULTIPLE_SAMPLE) {
            //multisample mode, only perform segmentation
            final SimpleIntervalCollection segments = new MultisampleMultidimensionalKernelSegmenter(
                    modelSegmentsData.denoisedCopyRatiosPerSample, modelSegmentsData.genotypingResult.getHetAllelicCountsPerSample())
                    .findSegmentation(
                            segmentationArguments.maxNumSegmentsPerChromosome,
                            segmentationArguments.kernelVarianceCopyRatio,
                            segmentationArguments.kernelVarianceAlleleFraction,
                            segmentationArguments.kernelScalingAlleleFraction,
                            segmentationArguments.kernelApproximationDimension,
                            ImmutableSet.copyOf(segmentationArguments.windowSizes).asList(),
                            segmentationArguments.numChangepointsPenaltyFactor,
                            segmentationArguments.numChangepointsPenaltyFactor);
            logHeapUsage("segmentation");

            final File segmentsIntervalListFile = new File(outputDir, outputPrefix + PICARD_INTERVAL_LIST_FILE_SUFFIX);
            logger.info(String.format("Writing segments as Picard interval list to %s...", segmentsIntervalListFile.getAbsolutePath()));
            final IntervalList segmentsIntervalList = new IntervalList(segments.getMetadata().getSequenceDictionary());
            segments.getIntervals().forEach(i -> segmentsIntervalList.add(new Interval(i)));
            segmentsIntervalList.write(segmentsIntervalListFile);
        } else {
            //single-sample mode

            //both denoisedCopyRatios and hetAllelicCounts are guaranteed to be non-null at this point;
            //missing data has already been imputed as an empty collection with the appropriate metadata
            final CopyRatioCollection denoisedCopyRatios = modelSegmentsData.denoisedCopyRatiosPerSample.get(0);
            final AllelicCountCollection hetAllelicCounts = modelSegmentsData.genotypingResult.getHetAllelicCountsPerSample().get(0);
            final SampleLocatableMetadata metadata = denoisedCopyRatios.getMetadata();

            //write allelic-counts files containing hets for the case and the matched-normal when available
            if (dataMode != DataMode.COPY_RATIO_ONLY) {
                final File hetAllelicCountsFile = new File(outputDir, outputPrefix + HET_ALLELIC_COUNTS_FILE_SUFFIX);
                if (inputNormalAllelicCountsFile == null) {
                    //case-only mode
                    logger.info(String.format("Writing heterozygous allelic counts to %s...", hetAllelicCountsFile.getAbsolutePath()));
                } else {
                    //matched-normal mode
                    final File hetNormalAllelicCountsFile = new File(outputDir, outputPrefix + NORMAL_HET_ALLELIC_COUNTS_FILE_SUFFIX);
                    logger.info(String.format("Writing heterozygous allelic counts for matched normal to %s...", hetNormalAllelicCountsFile.getAbsolutePath()));
                    final AllelicCountCollection hetNormalAllelicCounts = modelSegmentsData.genotypingResult.getHetNormalAllelicCounts();
                    hetNormalAllelicCounts.write(hetNormalAllelicCountsFile);

                    logger.info(String.format("Writing allelic counts for case sample at heterozygous sites in matched normal to %s...", hetAllelicCountsFile.getAbsolutePath()));
                }
                hetAllelicCounts.write(hetAllelicCountsFile);
            }

            final SimpleIntervalCollection segments;
            if (inputSegmentsFile == null) {
                segments = new MultisampleMultidimensionalKernelSegmenter(
                        Collections.singletonList(denoisedCopyRatios), Collections.singletonList(hetAllelicCounts))
                        .findSegmentation(
                                segmentationArguments.maxNumSegmentsPerChromosome,
                                segmentationArguments.kernelVarianceCopyRatio,
                                segmentationArguments.kernelVarianceAlleleFraction,
                                segmentationArguments.kernelScalingAlleleFraction,
                                segmentationArguments.kernelApproximationDimension,
                                ImmutableSet.copyOf(segmentationArguments.windowSizes).asList(),
                                segmentationArguments.numChangepointsPenaltyFactor,
                                segmentationArguments.numChangepointsPenaltyFactor);
                logHeapUsage("segmentation");
            } else {
                final IntervalList segmentsIntervalList = IntervalList.fromFile(inputSegmentsFile);
                ParamUtils.isPositive(segmentsIntervalList.size(),
                        "Segments file must contain at least one segment.");
                logger.info(String.format("Using input segmentation from %s containing %d segments...",
                        inputSegmentsFile, segmentsIntervalList.size()));
                segments = new SimpleIntervalCollection(
                        metadata,
                        segmentsIntervalList.getIntervals().stream()
                                .map(i -> new SimpleInterval(i.getContig(), i.getStart(), i.getEnd()))
                                .collect(Collectors.toList()));
            }

            logger.info("Modeling available denoised copy ratios and heterozygous allelic counts...");
            //initial MCMC model fitting performed by MultidimensionalModeller constructor
            final AlleleFractionPrior alleleFractionPrior = new AlleleFractionPrior(modelingArguments.minorAlleleFractionPriorAlpha);
            final MultidimensionalModeller modeller = new MultidimensionalModeller(
                    segments, denoisedCopyRatios, hetAllelicCounts, alleleFractionPrior,
                    modelingArguments.numSamplesCopyRatio, modelingArguments.numBurnInCopyRatio,
                    modelingArguments.numSamplesAlleleFraction, modelingArguments.numBurnInAlleleFraction);

            //write initial segments and parameters to file
            writeModeledSegmentsAndParameterFiles(modeller, BEGIN_FIT_FILE_TAG);

            //segmentation smoothing
            modeller.smoothSegments(
                    modelingArguments.maxNumSmoothingIterations, modelingArguments.numSmoothingIterationsPerFit,
                    modelingArguments.smoothingCredibleIntervalThresholdCopyRatio, modelingArguments.smoothingCredibleIntervalThresholdAlleleFraction);

            //write final segments and parameters to file
            writeModeledSegmentsAndParameterFiles(modeller, FINAL_FIT_FILE_TAG);

            logHeapUsage("modeling");

            //write final segments for copy-ratio caller (TODO remove this and MEAN_LOG2_COPY_RATIO column when new caller is available)
            final OverlapDetector<CopyRatio> copyRatioMidpointOverlapDetector = denoisedCopyRatios.getMidpointOverlapDetector();
            final CopyRatioSegmentCollection copyRatioSegmentsFinal = new CopyRatioSegmentCollection(
                    modeller.getModeledSegments().getMetadata(),
                    modeller.getModeledSegments().getIntervals().stream()
                            .map(s -> new CopyRatioSegment(s, new ArrayList<>(copyRatioMidpointOverlapDetector.getOverlaps(s))))
                            .collect(Collectors.toList()));
            writeSegments(copyRatioSegmentsFinal, COPY_RATIO_SEGMENTS_FOR_CALLER_FILE_SUFFIX);

            //write IGV-compatible files
            final LegacySegmentCollection copyRatioLegacySegments = new LegacySegmentCollection(
                    metadata,
                    modeller.getModeledSegments().getRecords().stream()
                            .map(s -> new LegacySegment(
                                    metadata.getSampleName(),
                                    s.getInterval(),
                                    s.getNumPointsCopyRatio(),
                                    s.getLog2CopyRatioSimplePosteriorSummary().getDecile50()))
                            .collect(Collectors.toList()));
            final LegacySegmentCollection alleleFractionLegacySegments = new LegacySegmentCollection(
                    metadata,
                    modeller.getModeledSegments().getRecords().stream()
                            .map(s -> new LegacySegment(
                                    metadata.getSampleName(),
                                    s.getInterval(),
                                    s.getNumPointsAlleleFraction(),
                                    s.getMinorAlleleFractionSimplePosteriorSummary().getDecile50()))
                            .collect(Collectors.toList()));
            writeSegments(copyRatioLegacySegments, COPY_RATIO_LEGACY_SEGMENTS_FILE_SUFFIX);
            writeSegments(alleleFractionLegacySegments, ALLELE_FRACTION_LEGACY_SEGMENTS_FILE_SUFFIX);
        }

        logger.info(String.format("%s complete.", getClass().getSimpleName()));

        return null;
    }

    /**
     * Performs reading and validation of data, as well as genotyping (if allele-fraction data is available).
     *
     * Regardless of {@code dataMode}, the respective lists for the copy-ratio and allele-fraction data for the
     * case samples will be of the same length and populated; empty collections will be used to impute missing data
     * using the same sample order as the non-missing data. This is done so that a
     * {@link MultisampleMultidimensionalKernelSegmenter} can be used in all scenarios.
     *
     * In contrast, {@code normalAllelicCounts} and {@code segments} will be set to {@code null} if the respective
     * inputs are missing and subsequent control-flow code  for toggling genotyping and segmentation must be aware of
     * this convention.
     */
    private final class ModelSegmentsData {
        final ImmutableList<CopyRatioCollection> denoisedCopyRatiosPerSample;
        final ImmutableList<AllelicCountCollection> allelicCountsPerSample;
        final AllelicCountCollection normalAllelicCounts;
        final SimpleIntervalCollection segments;
        final NaiveHeterozygousPileupGenotypingUtils.NaiveHeterozygousPileupGenotypingResult genotypingResult;

        private ModelSegmentsData() {
            Utils.nonNull(dataMode);

            //read available data and impute metadata for missing data
            switch (dataMode) {
                case COPY_RATIO_ONLY:
                    denoisedCopyRatiosPerSample = ImmutableList.copyOf(inputDenoisedCopyRatiosFiles.stream()
                            .map(CopyRatioCollection::new)
                            .iterator());
                    allelicCountsPerSample = ImmutableList.copyOf(denoisedCopyRatiosPerSample.stream()
                            .map(CopyRatioCollection::getMetadata)
                            .map(m -> new AllelicCountCollection(m, Collections.emptyList()))
                            .iterator());
                    break;
                case ALLELE_FRACTION_ONLY:
                    allelicCountsPerSample = ImmutableList.copyOf(inputAllelicCountsFiles.stream()
                            .map(AllelicCountCollection::new)
                            .iterator());
                    denoisedCopyRatiosPerSample = ImmutableList.copyOf(allelicCountsPerSample.stream()
                            .map(AllelicCountCollection::getMetadata)
                            .map(m -> new CopyRatioCollection(m, Collections.emptyList()))
                            .iterator());
                    break;
                case COPY_RATIO_AND_ALLELE_FRACTION:
                    denoisedCopyRatiosPerSample = ImmutableList.copyOf(inputDenoisedCopyRatiosFiles.stream()
                            .map(f -> readOptionalFileOrNull(f, CopyRatioCollection::new))
                            .iterator());
                    allelicCountsPerSample = ImmutableList.copyOf(inputAllelicCountsFiles.stream()
                            .map(f -> readOptionalFileOrNull(f, AllelicCountCollection::new))
                            .iterator());
                    //check that sample metadata matches across copy-ratio and allele-fraction data
                    IntStream.range(0, inputDenoisedCopyRatiosFiles.size()).boxed()
                            .forEach(i -> CopyNumberArgumentValidationUtils.getValidatedMetadata(
                                    denoisedCopyRatiosPerSample.get(i),
                                    allelicCountsPerSample.get(i)));
                    break;
                default:
                    throw new GATKException.ShouldNeverReachHereException("Unknown DataMode.");
            }

            normalAllelicCounts = readOptionalFileOrNull(inputNormalAllelicCountsFile, AllelicCountCollection::new);

            final IntervalList segmentsIntervalList = readOptionalFileOrNull(inputSegmentsFile, IntervalList::fromFile);
            segments = segmentsIntervalList == null
                    ? null
                    : new SimpleIntervalCollection(
                            new SimpleLocatableMetadata(segmentsIntervalList.getHeader().getSequenceDictionary()),
                            segmentsIntervalList.getIntervals().stream()
                                    .map(i -> new SimpleInterval(i.getContig(), i.getStart(), i.getEnd()))
                                    .collect(Collectors.toList()));

            logHeapUsage("reading files");

            final SAMSequenceDictionary sequenceDictionary = CopyNumberArgumentValidationUtils.getValidatedSequenceDictionary(
                    Stream.of(denoisedCopyRatiosPerSample, allelicCountsPerSample, Arrays.asList(normalAllelicCounts, segments))
                            .flatMap(Collection::stream)
                            .toArray(AbstractLocatableCollection[]::new));

            Utils.validateArg((int) denoisedCopyRatiosPerSample.stream()
                            .map(CopyRatioCollection::getIntervals)
                            .distinct()
                            .count() == 1,
                    "Copy-ratio intervals must be identical across all case samples.");

            Utils.validateArg((int) Stream.of(allelicCountsPerSample, Collections.singletonList(normalAllelicCounts))
                            .flatMap(Collection::stream)
                            .filter(Objects::nonNull)
                            .map(AllelicCountCollection::getIntervals)
                            .distinct()
                            .count() == 1,
                    "Allelic-count sites must be identical across all samples.");

            logHeapUsage("validating data");

            final SimpleIntervalCollection copyRatioIntervals = new SimpleIntervalCollection(
                    new SimpleLocatableMetadata(sequenceDictionary),
                    denoisedCopyRatiosPerSample.get(0).getIntervals());
            genotypingResult = NaiveHeterozygousPileupGenotypingUtils.genotypeHets(
                    allelicCountsPerSample, normalAllelicCounts, copyRatioIntervals,
                    genotypingArguments.minTotalAlleleCountCase,
                    genotypingArguments.minTotalAlleleCountNormal,
                    genotypingArguments.genotypingHomozygousLogRatioThreshold,
                    genotypingArguments.genotypingBaseErrorRate);

            logHeapUsage("genotyping");
        }
    }

    private void setModesAndValidateArguments() {
        CopyNumberArgumentValidationUtils.validateInputs(
                Stream.of(
                        inputDenoisedCopyRatiosFiles,
                        inputAllelicCountsFiles,
                        Collections.singletonList(inputNormalAllelicCountsFile),
                        Collections.singletonList(inputSegmentsFile))
                        .flatMap(Collection::stream)
                        .toArray(File[]::new));
        Utils.nonEmpty(outputPrefix);
        CopyNumberArgumentValidationUtils.validateAndPrepareOutputDirectories(outputDir);

        Utils.validateArg(!(inputDenoisedCopyRatiosFiles.isEmpty() && inputAllelicCountsFiles.isEmpty()),
                "Must provide at least one denoised-copy-ratios file or allelic-counts file.");
        Utils.validateArg(!(inputAllelicCountsFiles.isEmpty() && inputNormalAllelicCountsFile != null),
                "Must provide an allelic-counts file for the case sample to run in matched-normal mode.");

        runMode = (inputDenoisedCopyRatiosFiles.size() > 1 || inputAllelicCountsFiles.size() > 1)
                ? RunMode.MULTIPLE_SAMPLE
                : RunMode.SINGLE_SAMPLE;

        if (runMode == RunMode.MULTIPLE_SAMPLE) {
            if (!inputDenoisedCopyRatiosFiles.isEmpty() && !inputAllelicCountsFiles.isEmpty()) {
                Utils.validateArg(inputDenoisedCopyRatiosFiles.size() == inputAllelicCountsFiles.size(),
                        "Number of denoised-copy-ratios and allelic-counts files for the case samples must be equal " +
                                "if both input types are specified in multisample mode.");
            }
            Utils.validateArg(inputSegmentsFile == null,
                    "Segments file cannot be specified in multisample mode.");
        }

        if (inputAllelicCountsFiles.isEmpty()) {
            dataMode = DataMode.COPY_RATIO_ONLY;
        } else if (inputDenoisedCopyRatiosFiles.isEmpty()) {
            dataMode = DataMode.ALLELE_FRACTION_ONLY;
        } else {
            dataMode = DataMode.COPY_RATIO_AND_ALLELE_FRACTION;
        }

        modelingArguments.validateArguments();
    }

    private <T> T readOptionalFileOrNull(final File file,
                                         final Function<File, T> read) {
        if (file == null) {
            return null;
        }
        logger.info(String.format("Reading file (%s)...", file));
        return read.apply(file);
    }

    private void writeModeledSegmentsAndParameterFiles(final MultidimensionalModeller modeller,
                                                       final String fileTag) {
        final ModeledSegmentCollection modeledSegments = modeller.getModeledSegments();
        writeSegments(modeledSegments, fileTag + SEGMENTS_FILE_SUFFIX);
        final File copyRatioParameterFile = new File(outputDir, outputPrefix + fileTag + COPY_RATIO_MODEL_PARAMETER_FILE_SUFFIX);
        final File alleleFractionParameterFile = new File(outputDir, outputPrefix + fileTag + ALLELE_FRACTION_MODEL_PARAMETER_FILE_SUFFIX);
        modeller.writeModelParameterFiles(copyRatioParameterFile, alleleFractionParameterFile);
    }

    private void writeSegments(final AbstractRecordCollection<?, ?> segments,
                               final String fileSuffix) {
        final File segmentsFile = new File(outputDir, outputPrefix + fileSuffix);
        logger.info(String.format("Writing segments to %s...", segmentsFile.getAbsolutePath()));
        segments.write(segmentsFile);
    }
}
