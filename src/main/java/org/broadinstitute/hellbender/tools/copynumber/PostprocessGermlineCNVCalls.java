package org.broadinstitute.hellbender.tools.copynumber;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.BaselineCopyNumberCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.CopyNumberPosteriorDistributionCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.IntegerCopyNumberSegmentCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.SimpleIntervalCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.LocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CopyNumberPosteriorDistribution;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.LocatableCopyNumberPosteriorDistribution;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.LocatableIntegerCopyNumber;
import org.broadinstitute.hellbender.tools.copynumber.gcnv.*;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.io.Resource;
import org.broadinstitute.hellbender.utils.python.PythonScriptExecutor;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Postprocesses the output of {@link GermlineCNVCaller} and generates VCF files.
 *
 * <p>Depending on the arguments, this tool either generates a single "intervals" VCF or additionally, performs
 * segmentation on intervals (i.e. bins or targets) and generates genotyped constant copy-number segments.</p>
 *
 * <p>The intervals VCF provides a detailed listing of the most likely copy-number call for each genomic
 * interval included in the call-set, along with call quality, call genotype, and the phred-scaled posterior
 * probability vector of all integer copy-number states.</p>
 *
 * <p>CNV events often span several consecutive intervals such that it may be desirable to coalesce consecutive
 * calls into constant copy-number segments. This tool further performs segmentation and genotyping by calling a
 * dedicated python script in `gcnvkernel`. The segmentation algorithm also provides various quality metrics.
 * Please refer to {@link GermlineCNVSegmentVariantComposer} for a description of the segmentation quality metrics.
 *
 * <p>For both VCF outputs, the CNV genotype is determined as follows: the alternative allele for a CNV call is
 * either set to <code>&lt;DEL&gt;</code> or <code>&lt;DUP&gt;</code>, depending on whether the most likely
 * copy-number call is below or above the reference copy-number of the associated contig. The user may specify the
 * reference copy-number state on autosomal contigs using the argument <code>autosomal-ref-copy-number</code>.
 * The list of allosomal contigs may also be specified via the argument <code>allosomal-contig</code>. All
 * undeclared contigs are assumed to be autosomal. The reference copy-number on an allosomal contig reflects
 * the sex karyotype of the sample and is set to the pre-determined contig ploidy state fetched from the output
 * of {@link DetermineGermlineContigPloidy}.</p>
 *
 * <h3>Required inputs for intervals VCF generation:</h3>
 * <ul>
 *     <li>A list of paths to {@link GermlineCNVCaller} calls shards</li>
 *     <li>Index of the sample in the call-set (which is expected to be the same across all shards)</li>
 *     <li>Output path for writing the intervals VCF</li>
 * </ul>
 *
 * <h3>Required inputs for intervals and segments VCF generation:</h3>
 * <ul>
 *     <li>A list of paths to {@link GermlineCNVCaller} calls shards</li>
 *     <li>A list of paths to {@link GermlineCNVCaller} model shards</li>
 *     <li>Path to the output calls of {@link DetermineGermlineContigPloidy}</li>
 *     <li>Index of the sample in the call-set (which is expected to be the same across all shards)</li>
 *     <li>Output path for writing the intervals VCF</li>
 *     <li>Output path for writing the segments VCF</li>
 * </ul>
 *
 * <p>The order of specifying {@link GermlineCNVCaller} calls or model shards is immaterial and will be put
 * in the ascending order according to the call-set SAM sequence dictionary.</p>
 *
 * <h3>Usage example</h3>
 *
 * <p>To generate intervals VCF only:</p>
 * <pre>
 *   gatk PostprocessGermlineCNVCalls \
 *     --calls-shard-path path/to/shard_1-calls
 *     --calls-shard-path path/to/shard_2-calls
 *     --sample-index 0
 *     --autosomal-ref-copy-number 2
 *     --allosomal-contig X
 *     --allosomal-contig Y
 *     --output-genotyped-intervals sample_0_genotyped_intervals.vcf
 * </pre>
 *
 * <p>To generate boths intervals and segments VCF:</p>
 * <pre>
 *   gatk PostprocessGermlineCNVCalls \
 *     --calls-shard-path path/to/shard_1-calls
 *     --calls-shard-path path/to/shard_2-calls
 *     --model-shard-path path/to/shard_1-model
 *     --model-shard-path path/to/shard_2-model
 *     --sample-index 0
 *     --autosomal-ref-copy-number 2
 *     --allosomal-contig X
 *     --allosomal-contig Y
 *     --output-genotyped-intervals sample_0_genotyped_intervals.vcf
 *     --output-genotyped-segments sample_0_genotyped_segments.vcf
 * </pre>
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 * @author Andrey Smirnov &lt;asmirnov@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Postprocesses the output of GermlineCNVCaller and generates VCF files.",
        oneLineSummary = "Postprocesses the output of GermlineCNVCaller and generates VCF files.",
        programGroup = CopyNumberProgramGroup.class
)
@DocumentedFeature
@BetaFeature
public final class PostprocessGermlineCNVCalls extends GATKTool {
    private static final Logger logger = LogManager.getLogger(PostprocessGermlineCNVCalls.class);

    public static final String SEGMENT_GERMLINE_CNV_CALLS_PYTHON_SCRIPT = "segment_gcnv_calls.py";

    public static final String CALLS_SHARD_PATH_LONG_NAME = "calls-shard-path";
    public static final String MODEL_SHARD_PATH_LONG_NAME = "model-shard-path";
    public static final String CONTIG_PLOIDY_CALLS_LONG_NAME = "contig-ploidy-calls";
    public static final String SAMPLE_INDEX_LONG_NAME = "sample-index";
    public static final String OUTPUT_INTERVALS_VCF_LONG_NAME = "output-genotyped-intervals";
    public static final String OUTPUT_SEGMENTS_VCF_LONG_NAME = "output-genotyped-segments";
    public static final String DRY_RUN_LONG_NAME = "dry-run";
    public static final String AUTOSOMAL_REF_COPY_NUMBER_LONG_NAME = "autosomal-ref-copy-number";
    public static final String ALLOSOMAL_CONTIG_LONG_NAME = "allosomal-contig";

    @Argument(
            doc = "List of paths to GermlineCNVCaller call directories.",
            fullName = CALLS_SHARD_PATH_LONG_NAME,
            minElements = 1
    )
    private List<File> unsortedCallsShardPaths;

    @Argument(
            doc = "List of paths to GermlineCNVCaller model directories.",
            fullName = MODEL_SHARD_PATH_LONG_NAME,
            minElements = 1,
            optional = true
    )
    private List<File> unsortedModelShardPaths;

    @Argument(
            doc = "Path to contig-ploidy calls directory (output of DetermineGermlineContigPloidy).",
            fullName = CONTIG_PLOIDY_CALLS_LONG_NAME,
            optional = true
    )
    private File contigPloidyCallsPath = null;

    @Argument(
            doc = "Sample index in the call-set (must be contained in all shards).",
            fullName = SAMPLE_INDEX_LONG_NAME,
            minValue = 0
    )
    private int sampleIndex;

    @Argument(
            doc = "Reference copy-number on autosomal intervals.",
            fullName = AUTOSOMAL_REF_COPY_NUMBER_LONG_NAME,
            minValue = 0
    )
    private int refAutosomalCopyNumber = 2;

    @Argument(
            doc = "Contigs to treat as allosomal (hence, set their reference copy-number allele according to " +
                    "their neutral baseline ploidy state).",
            fullName = ALLOSOMAL_CONTIG_LONG_NAME,
            optional = true
    )
    private List<String> allosomalContigList;

    @Argument(
            doc = "Output per-interval VCF file.",
            fullName = OUTPUT_INTERVALS_VCF_LONG_NAME
    )
    private File outputIntervalsVCFFile;

    @Argument(
            doc = "Output segments VCF file.",
            fullName = OUTPUT_SEGMENTS_VCF_LONG_NAME,
            optional = true
    )
    private File outputSegmentsVCFFile = null;

    @Argument(
            doc = "Dry run (input validation only).",
            fullName = DRY_RUN_LONG_NAME,
            optional = true
    )
    private boolean dryRun = false;

    /**
     * A list of {@link SimpleIntervalCollection} for each shard
     */
    private List<SimpleIntervalCollection> sortedIntervalCollections;

    /**
     * The sample name corresponding to the provided sample index
     */
    private String sampleName;

    /**
     * Number of shards
     */
    private int numShards;

    /**
     * Reference integer copy-number state for autosomal intervals
     */
    private IntegerCopyNumberState refAutosomalIntegerCopyNumberState;

    /**
     * The collection of all integer copy-number states in the call-set
     */
    private IntegerCopyNumberStateCollection integerCopyNumberStateCollection;

    /**
     * Set of allosomal contigs
     */
    private Set<String> allosomalContigSet;

    /**
     * Call shard directories put in correct order
     */
    private List<File> sortedCallsShardPaths;

    /**
     * Model shard directories put in correct order
     */
    private List<File> sortedModelShardPaths;

    @Override
    public void onStartup() {
        super.onStartup();
        if (!dryRun && outputSegmentsVCFFile != null) {
            /* check for successful import of gcnvkernel */
            PythonScriptExecutor.checkPythonEnvironmentForPackage("gcnvkernel");
        }
    }

    /**
     * Performs various validations on input arguments. Since many of these validations requires loading and parsing
     * reusable data, we store them as global variables (shard interval lists, sample name, etc).
     */
    @Override
    public void onTraversalStart() {
        numShards = unsortedCallsShardPaths.size();

        /* get intervals from each call and model shard in the provided (potentially arbitrary) order */
        final List<SimpleIntervalCollection> unsortedIntervalCollectionsFromCalls =
                getIntervalCollectionsFromPaths(unsortedCallsShardPaths);

        /* assert that all shards have the same SAM sequence dictionary */
        final SAMSequenceDictionary samSequenceDictionary = unsortedIntervalCollectionsFromCalls.get(0)
                .getMetadata().getSequenceDictionary();
        Utils.validateArg(unsortedIntervalCollectionsFromCalls.stream()
                        .map(SimpleIntervalCollection::getMetadata)
                        .map(LocatableMetadata::getSequenceDictionary)
                        .allMatch(shardSAMSequenceDictionary ->
                                shardSAMSequenceDictionary.equals(samSequenceDictionary)),
                "The SAM sequence dictionary is not the same for all of the call shards.");

        /* get the correct shard sort order and sort all collections */
        final int[] sortedCallShardsOrder = SimpleIntervalCollection.getSimpleIntervalCollectionSortedOrder(
                unsortedIntervalCollectionsFromCalls, samSequenceDictionary, true);
        sortedCallsShardPaths = Arrays.stream(sortedCallShardsOrder)
                .mapToObj(unsortedCallsShardPaths::get)
                .collect(Collectors.toList());
        final List<SimpleIntervalCollection> sortedIntervalCollectionsFromCalls = Arrays.stream(sortedCallShardsOrder)
                .mapToObj(unsortedIntervalCollectionsFromCalls::get)
                .collect(Collectors.toList());
        sortedIntervalCollections = sortedIntervalCollectionsFromCalls;

        /* assert that allosomal contigs are contained in the SAM sequence dictionary */
        final Set<String> allContigs = samSequenceDictionary.getSequences().stream()
                .map(SAMSequenceRecord::getSequenceName)
                .collect(Collectors.toSet());
        allosomalContigSet = new HashSet<>(allosomalContigList);
        if (allosomalContigSet.isEmpty()) {
            logger.warn(String.format("Allosomal contigs were not specified; setting ref copy-number allele " +
                    "to (%d) for all intervals.", refAutosomalCopyNumber));
        } else {
            Utils.validateArg(allContigs.containsAll(allosomalContigSet), String.format(
                    "The specified allosomal contigs must be contained in the SAM sequence dictionary of the " +
                            "call-set (specified allosomal contigs: %s, all contigs: %s)",
                    allosomalContigSet.stream().collect(Collectors.joining(", ", "[", "]")),
                    allContigs.stream().collect(Collectors.joining(", ", "[", "]"))));
        }

        /* get sample name from the first shard and assert that all shards have the same sample name */
        sampleName = getShardSampleName(0);
        Utils.validate(IntStream.range(1, numShards)
                        .mapToObj(this::getShardSampleName)
                        .allMatch(shardSampleName -> shardSampleName.equals(sampleName)),
                "The sample name is not the same for all of the shards.");

        /* get a sorted list of contained integer copy-number states */
        integerCopyNumberStateCollection = getCopyNumberStateCollectionFromCopyNumberPosteriorsHeader(0);
        Utils.validate(integerCopyNumberStateCollection.size() >= 3, "There should be at least 3 copy-number states.");
        Utils.validate(IntStream.range(1, numShards)
                        .mapToObj(this::getCopyNumberStateCollectionFromCopyNumberPosteriorsHeader)
                        .allMatch(shardIntegerCopyNumberStateList ->
                                shardIntegerCopyNumberStateList.equals(integerCopyNumberStateCollection)),
                "The copy-number states list is not the same for all of the shards.");

        /* assert that the ref autosomal copy-number is contained in the integer copy-number collection of the call-set */
        refAutosomalIntegerCopyNumberState = new IntegerCopyNumberState(refAutosomalCopyNumber);
        Utils.validateArg(integerCopyNumberStateCollection.getCopyNumberStates().contains(refAutosomalIntegerCopyNumberState),
                String.format("The provided ref copy-number state of autosomal intervals (%s) is not " +
                                "contained by the integer copy-number state collection of the call-set (%s).",
                        refAutosomalIntegerCopyNumberState, integerCopyNumberStateCollection));

        /* if segments VCF is required, further validations are necessary on the inputs of the segmenter python script */
        if (outputSegmentsVCFFile != null) {
            Utils.nonNull(contigPloidyCallsPath, "Contig ploidy calls are required for generating segments " +
                    "VCF file.");
            Utils.validateArg(unsortedModelShardPaths.size() == numShards,
                    "The number of input model shards must match the number of input call shards.");

            /* basic assertion that the model shards are valid and are compatible with call shards */
            final List<SimpleIntervalCollection> unsortedIntervalCollectionsFromModels =
                    getIntervalCollectionsFromPaths(unsortedModelShardPaths);

            Utils.validateArg(unsortedIntervalCollectionsFromModels.stream()
                            .map(SimpleIntervalCollection::getMetadata)
                            .map(LocatableMetadata::getSequenceDictionary)
                            .allMatch(shardSAMSequenceDictionary ->
                                    shardSAMSequenceDictionary.equals(samSequenceDictionary)),
                    "The SAM sequence dictionary is either not the same for all of the model shards, " +
                            "or is different from the SAM sequence dictionary of calls shards.");

            final int[] sortedModelShardsOrder = SimpleIntervalCollection.getSimpleIntervalCollectionSortedOrder(
                    unsortedIntervalCollectionsFromModels, samSequenceDictionary, true);
            final List<SimpleIntervalCollection> sortedIntervalCollectionsFromModels = Arrays.stream(sortedModelShardsOrder)
                    .mapToObj(unsortedIntervalCollectionsFromModels::get)
                    .collect(Collectors.toList());
            Utils.validateArg(sortedIntervalCollectionsFromCalls.equals(sortedIntervalCollectionsFromModels),
                    "The interval lists found in model and call shards do not match. Make sure that the calls and model " +
                            "paths are provided in matching order.");

            /* sort model shards */
            sortedModelShardPaths = Arrays.stream(sortedModelShardsOrder)
                    .mapToObj(unsortedModelShardPaths::get)
                    .collect(Collectors.toList());
        }
    }

    @Override
    public void traverse() {
        if (!dryRun) {
            generateIntervalsVCFFileFromAllShards();
            generateSegmentsVCFFileFromAllShards();
        } else {
            logger.info("Dry run -- only input validations were performed.");
        }
    }

    private void generateIntervalsVCFFileFromAllShards() {
        logger.info("Generating intervals VCF file...");
        final VariantContextWriter intervalsVCFWriter = createVCFWriter(outputIntervalsVCFFile);

        final GermlineCNVIntervalVariantComposer germlineCNVIntervalVariantComposer =
                new GermlineCNVIntervalVariantComposer(intervalsVCFWriter, sampleName,
                        integerCopyNumberStateCollection, refAutosomalIntegerCopyNumberState, allosomalContigSet);
        germlineCNVIntervalVariantComposer.composeVariantContextHeader(getDefaultToolVCFHeaderLines());

        for (int shardIndex = 0; shardIndex < numShards; shardIndex++) {
            logger.info(String.format("Analyzing shard %d...", shardIndex));
            germlineCNVIntervalVariantComposer.writeVariantContext(
                    getShardCopyNumberPosteriorDistributions(shardIndex),
                    getShardBaselineCopyNumbers(shardIndex));
        }
        intervalsVCFWriter.close();
    }

    /**
     * Retrieves intervals either from a `gcnvkernel` model output path or the root directory of calls path.
     */
    private static List<SimpleIntervalCollection> getIntervalCollectionsFromPaths(final List<File> shardsPathList) {
        return shardsPathList.stream()
                .map(shardPath -> new SimpleIntervalCollection(
                        getIntervalFileFromShardDirectory(shardPath)))
                .collect(Collectors.toList());
    }

    /**
     * Extracts sample name from a designated text file in the sample calls directory.
     */
    private String getShardSampleName(final int shardIndex) {
        final File shardSampleNameTextFile = getSampleNameTextFile(sortedCallsShardPaths.get(shardIndex),
                sampleIndex);
        try {
            final BufferedReader reader = new BufferedReader(new FileReader(shardSampleNameTextFile));
            return reader.readLine();
        } catch (final IOException ex) {
            throw new UserException.BadInput(String.format("Could not read the sample name text file at %s.",
                    shardSampleNameTextFile.getAbsolutePath()));
        }
    }

    /**
     * Returns a list of {@link LocatableCopyNumberPosteriorDistribution} for {@link #sampleIndex} from a
     * single calls shard.
     */
    private List<LocatableCopyNumberPosteriorDistribution> getShardCopyNumberPosteriorDistributions(final int shardIndex) {
        /* read copy-number posteriors for the shard */
        final File shardRootDirectory = sortedCallsShardPaths.get(shardIndex);
        final File sampleCopyNumberPosteriorFile = getSampleCopyNumberPosteriorFile(shardRootDirectory, sampleIndex);
        final CopyNumberPosteriorDistributionCollection copyNumberPosteriorDistributionCollection =
                new CopyNumberPosteriorDistributionCollection(sampleCopyNumberPosteriorFile, integerCopyNumberStateCollection);
        final String sampleNameFromCopyNumberPosteriorFile = copyNumberPosteriorDistributionCollection
                .getMetadata().getSampleName();
        Utils.validate(sampleNameFromCopyNumberPosteriorFile.equals(sampleName),
                String.format("Sample name found in the header of copy-number posterior file for shard %d " +
                        "different from the expected sample name (found: %s, expected: %s).",
                        shardIndex, sampleNameFromCopyNumberPosteriorFile, sampleName));

        /* attach the intervals to make locatable posteriors */
        final List<SimpleInterval> shardIntervals = sortedIntervalCollections.get(shardIndex).getIntervals();
        final List<CopyNumberPosteriorDistribution> copyNumberPosteriorDistributionList =
                copyNumberPosteriorDistributionCollection.getRecords();
        Utils.validate(shardIntervals.size() == copyNumberPosteriorDistributionList.size(),
                String.format("The number of entries in the copy-number posterior file for shard %d does " +
                        "not match the number of entries in the shard interval list (posterior list size: %d, " +
                        "interval list size: %d)", shardIndex, copyNumberPosteriorDistributionList.size(),
                        shardIntervals.size()));

        return IntStream.range(0, copyNumberPosteriorDistributionCollection.size())
                        .mapToObj(intervalIndex -> new LocatableCopyNumberPosteriorDistribution(
                                shardIntervals.get(intervalIndex),
                                copyNumberPosteriorDistributionList.get(intervalIndex)))
                        .collect(Collectors.toList());
    }

    /**
     * Returns the baseline copy-number states as a list of {@link LocatableIntegerCopyNumber} for
     * {@link #sampleIndex} from a single calls shard.
     */
    private List<LocatableIntegerCopyNumber> getShardBaselineCopyNumbers(final int shardIndex) {
        /* read copy-number baselines for the shard */
        final File shardRootDirectory = sortedCallsShardPaths.get(shardIndex);
        final File sampleBaselineCopyNumberFile = getSampleBaselineCopyNumberFile(shardRootDirectory, sampleIndex);
        final BaselineCopyNumberCollection baselineCopyNumberCollection = new BaselineCopyNumberCollection(sampleBaselineCopyNumberFile);
        final String sampleNameFromBaselineCopyNumberFile = baselineCopyNumberCollection
                .getMetadata().getSampleName();
        Utils.validate(sampleNameFromBaselineCopyNumberFile.equals(sampleName),
                String.format("Sample name found in the header of baseline copy-number file for shard %d " +
                                "different from the expected sample name (found: %s, expected: %s).",
                        shardIndex, sampleNameFromBaselineCopyNumberFile, sampleName));

        /* attach the intervals to make locatable posteriors */
        final List<SimpleInterval> shardIntervals = sortedIntervalCollections.get(shardIndex).getIntervals();
        final List<IntegerCopyNumberState> baselineCopyNumberList = baselineCopyNumberCollection.getRecords();
        return IntStream.range(0, baselineCopyNumberCollection.size())
                .mapToObj(intervalIndex -> new LocatableIntegerCopyNumber(
                        shardIntervals.get(intervalIndex),
                        baselineCopyNumberList.get(intervalIndex)))
                .collect(Collectors.toList());
    }

    /**
     * Returns a list of {@link IntegerCopyNumberState} from the table column names of the copy-number
     * posterior TSV file.
     */
    private IntegerCopyNumberStateCollection getCopyNumberStateCollectionFromCopyNumberPosteriorsHeader(final int shardIndex) {
        final File shardRootDirectory = sortedCallsShardPaths.get(shardIndex);
        final File sampleCopyNumberPosteriorFile = getSampleCopyNumberPosteriorFile(shardRootDirectory, sampleIndex);
        final List<String> copyNumberStateColumns = CopyNumberPosteriorDistributionCollection
                .extractCopyNumberPosteriorFileColumns(sampleCopyNumberPosteriorFile);
        return new IntegerCopyNumberStateCollection(copyNumberStateColumns);
    }

    /**
     * Gets the copy-number posterior file from the shard directory for a given sample index.
     */
    private static File getSampleCopyNumberPosteriorFile(final File callShardPath, final int sampleIndex) {
        return Paths.get(
                callShardPath.getAbsolutePath(),
                GermlineCNVNamingConstants.SAMPLE_PREFIX + sampleIndex,
                GermlineCNVNamingConstants.COPY_NUMBER_POSTERIOR_FILE_NAME).toFile();
    }

    /**
     * Gets the copy-number posterior file from the shard directory for a given sample index.
     */
    private static File getSampleBaselineCopyNumberFile(final File callShardPath, final int sampleIndex) {
        return Paths.get(
                callShardPath.getAbsolutePath(),
                GermlineCNVNamingConstants.SAMPLE_PREFIX + sampleIndex,
                GermlineCNVNamingConstants.BASELINE_COPY_NUMBER_FILE_NAME).toFile();
    }

    /**
     * Gets the sample name text file from the shard directory for a given sample index.
     */
    private static File getSampleNameTextFile(final File callShardRootPath, final int sampleIndex) {
        return Paths.get(
                callShardRootPath.getAbsolutePath(),
                GermlineCNVNamingConstants.SAMPLE_PREFIX + sampleIndex,
                GermlineCNVNamingConstants.SAMPLE_NAME_TXT_FILE).toFile();
    }

    /**
     * Gets the interval list file from the shard directory
     */
    private static File getIntervalFileFromShardDirectory(final File shardPath) {
        return new File(shardPath, GermlineCNVNamingConstants.INTERVAL_LIST_FILE_NAME);
    }

    private void generateSegmentsVCFFileFromAllShards() {
        if (outputSegmentsVCFFile != null) {
            logger.info("Generating segments VCF file...");

            /* perform segmentation */
            final File pythonScriptOutputPath = IOUtils.tempDir("gcnv-segmented-calls", "");
            final boolean pythonScriptSucceeded = executeSegmentGermlineCNVCallsPythonScript(
                    sampleIndex, contigPloidyCallsPath, sortedCallsShardPaths, sortedModelShardPaths,
                    pythonScriptOutputPath);
            if (!pythonScriptSucceeded) {
                throw new UserException("Python return code was non-zero.");
            }

            /* parse segments */
            final File copyNumberSegmentsFile = getCopyNumberSegmentsFile(pythonScriptOutputPath, sampleIndex);
            final IntegerCopyNumberSegmentCollection integerCopyNumberSegmentCollection =
                    new IntegerCopyNumberSegmentCollection(copyNumberSegmentsFile);
            final String sampleNameFromSegmentCollection = integerCopyNumberSegmentCollection
                    .getMetadata().getSampleName();
            Utils.validate(sampleNameFromSegmentCollection.equals(sampleName),
                    String.format("Sample name found in the header of copy-number segments file is " +
                                    "different from the expected sample name (found: %s, expected: %s).",
                            sampleNameFromSegmentCollection, sampleName));

            /* write variants */
            final VariantContextWriter segmentsVCFWriter = createVCFWriter(outputSegmentsVCFFile);
            final GermlineCNVSegmentVariantComposer germlineCNVSegmentVariantComposer =
                    new GermlineCNVSegmentVariantComposer(segmentsVCFWriter, sampleName,
                            refAutosomalIntegerCopyNumberState, allosomalContigSet);
            germlineCNVSegmentVariantComposer.composeVariantContextHeader(getDefaultToolVCFHeaderLines());
            germlineCNVSegmentVariantComposer.writeAll(integerCopyNumberSegmentCollection);
            segmentsVCFWriter.close();
        } else {
            logger.info("No segments output VCF file was provided -- skipping segmentation.");
        }
    }

    /**
     * Runs the segmenter python scripts and returns the exit code (true for a successful python call).
     */
    private static boolean executeSegmentGermlineCNVCallsPythonScript(final int sampleIndex,
                                                                      final File contigPloidyCallsPath,
                                                                      final List<File> sortedCallDirectories,
                                                                      final List<File> sortedModelDirectories,
                                                                      final File pythonScriptOutputPath) {
        /* the inputs to this method are expected to be previously validated */
        try {
            Utils.nonNull(contigPloidyCallsPath);
            Utils.nonNull(sortedCallDirectories);
            Utils.nonNull(sortedModelDirectories);
            Utils.nonNull(pythonScriptOutputPath);
            sortedCallDirectories.forEach(Utils::nonNull);
            sortedModelDirectories.forEach(Utils::nonNull);
        } catch (final IllegalArgumentException ex) {
            throw new GATKException.ShouldNeverReachHereException(ex);
        }

        final PythonScriptExecutor executor = new PythonScriptExecutor(true);
        final List<String> arguments = new ArrayList<>();
        arguments.add("--ploidy_calls_path");
        arguments.add(contigPloidyCallsPath.getAbsolutePath());
        arguments.add("--model_shards");
        arguments.addAll(sortedModelDirectories.stream().map(File::getAbsolutePath).collect(Collectors.toList()));
        arguments.add("--calls_shards");
        arguments.addAll(sortedCallDirectories.stream().map(File::getAbsolutePath).collect(Collectors.toList()));
        arguments.add("--output_path");
        arguments.add(pythonScriptOutputPath.getAbsolutePath());
        arguments.add("--sample_index");
        arguments.add(String.valueOf(sampleIndex));

        return executor.executeScript(
                new Resource(SEGMENT_GERMLINE_CNV_CALLS_PYTHON_SCRIPT, PostprocessGermlineCNVCalls.class),
                null,
                arguments);
    }

    /**
     * Gets the copy-number segments file from the output of `segment_gcnv_calls.py` for a given sample index.
     */
    private static File getCopyNumberSegmentsFile(final File pythonSegmenterOutputPath, final int sampleIndex) {
        return Paths.get(
                pythonSegmenterOutputPath.getAbsolutePath(),
                GermlineCNVNamingConstants.SAMPLE_PREFIX + sampleIndex,
                GermlineCNVNamingConstants.COPY_NUMBER_SEGMENTS_FILE_NAME).toFile();
    }
}
