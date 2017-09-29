package org.broadinstitute.hellbender.engine.spark;

import com.google.cloud.genomics.dataflow.utils.GCSOptions;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLinePluginDescriptor;
import org.broadinstitute.hellbender.cmdline.GATKPlugin.GATKReadFilterPluginDescriptor;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.*;
import org.broadinstitute.hellbender.engine.TraversalParameters;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.engine.datasources.ReferenceWindowFunctions;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.engine.FeatureManager;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.WellformedReadFilter;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSink;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SequenceDictionaryUtils;
import org.broadinstitute.hellbender.utils.SerializableFunction;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadsWriteFormat;

import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * Base class for GATK spark tools that accept standard kinds of inputs (reads, reference, and/or intervals).
 * Centralizes handling of tool inputs to enforce consistency, reduce duplicated boilerplate code
 * in tools, and apply standardized validation (such as sequence dictionary validation).
 *
 * Spark tools that do not fit into this pattern should extend SparkCommandLineProgram directly instead of
 * this class.
 *
 * USAGE:
 *
 * -Tools must implement {@link #runTool}.
 *
 * -Tools should override {@link #requiresReference}, {@link #requiresReads}, and/or {@link #requiresIntervals}
 *  as appropriate to indicate required inputs.
 *
 * -Tools can query whether certain inputs are present via {@link #hasReference}, {@link #hasReads}, and
 *  {@link #hasIntervals}.
 *
 * -Tools can load the reads via {@link #getReads}, access the reference via {@link #getReference}, and
 *  access the intervals via {@link #getIntervals}. Any intervals specified are automatically applied
 *  to the reads. Input metadata is available via {@link #getHeaderForReads}, {@link #getReferenceSequenceDictionary},
 *  and {@link #getBestAvailableSequenceDictionary}.
 *
 * -Tools that require a custom reference window function (extra bases of reference context around each read)
 *  may override {@link #getReferenceWindowFunction} to supply one. This function will be propagated to the
 *  reference source returned by {@link #getReference}.
 */
public abstract class GATKSparkTool extends SparkCommandLineProgram {
    private static final long serialVersionUID = 1L;

    @ArgumentCollection
    public final ReferenceInputArgumentCollection referenceArguments = requiresReference() ? new RequiredReferenceInputArgumentCollection() :  new OptionalReferenceInputArgumentCollection();

    @ArgumentCollection
    public final ReadInputArgumentCollection readArguments = requiresReads() ? new RequiredReadInputArgumentCollection() : new OptionalReadInputArgumentCollection();

    @ArgumentCollection
    protected IntervalArgumentCollection intervalArgumentCollection = requiresIntervals() ? new RequiredIntervalArgumentCollection() : new OptionalIntervalArgumentCollection();

    @Argument(doc = "maximum number of bytes to read from a file into each partition of reads. " +
            "Setting this higher will result in fewer partitions. Note that this will not be equal to the size of the partition in memory. " +
            "Defaults to 0, which uses the default split size (determined by the Hadoop input format, typically the size of one HDFS block).",
            fullName = "bamPartitionSize", shortName = "bps", optional = true)
    protected long bamPartitionSplitSize = 0;

    @Argument(fullName = StandardArgumentDefinitions.DISABLE_SEQUENCE_DICT_VALIDATION_NAME, shortName = StandardArgumentDefinitions.DISABLE_SEQUENCE_DICT_VALIDATION_NAME, doc = "If specified, do not check the sequence dictionaries from our inputs for compatibility. Use at your own risk!", optional = true)
    private boolean disableSequenceDictionaryValidation = false;

    @Argument(doc = "For tools that write an output, write the output in multiple pieces (shards)", shortName = "shardedOutput", fullName = "shardedOutput", optional = true)
    protected boolean shardedOutput = false;

    @Argument(doc="For tools that shuffle data or write an output, sets the number of reducers. Defaults to 0, which gives one partition per 10MB of input.",
            shortName = "numReducers", fullName = "numReducers", optional = true)
    protected int numReducers = 0;

    private ReadsSparkSource readsSource;
    private SAMFileHeader readsHeader;
    private String readInput;
    private ReferenceMultiSource referenceSource;
    private SAMSequenceDictionary referenceDictionary;
    private List<SimpleInterval> intervals;
    protected FeatureManager features;

    /**
     * Return the list of GATKCommandLinePluginDescriptor objects to be used for this CLP.
     * Use the read filter plugin.
     */
    @Override
    public List<? extends CommandLinePluginDescriptor<?>> getPluginDescriptors() {
        return Collections.singletonList(new GATKReadFilterPluginDescriptor(getDefaultReadFilters()));
    }

    /**
     * Does this tool require reference data? Tools that do should override to return true.
     *
     * @return true if this tool requires a reference, otherwise false
     */
    public boolean requiresReference() {
        return false;
    }

    /**
     * Does this tool require reads? Tools that do should override to return true.
     *
     * @return true if this tool requires reads, otherwise false
     */
    public boolean requiresReads() {
        return false;
    }

    /**
     * Does this tool require intervals? Tools that do should override to return true.
     *
     * @return true if this tool requires intervals, otherwise false
     */
    public boolean requiresIntervals() {
        return false;
    }

    /**
     * Is a source of reference data available?
     *
     * @return true if a reference is available, otherwise false
     */
    public final boolean hasReference() {
        return referenceSource != null;
    }

    /**
     * Are sources of reads available?
     *
     * @return true if reads are available, otherwise false
     */
    public final boolean hasReads() {
        return readsSource != null;
    }

    /**
     * Are sources of intervals available?
     *
     * @return true if intervals are available, otherwise false
     */
    public final boolean hasIntervals() {
        return intervals != null;
    }

    /**
     * Window function that controls how much reference context to return for each read when
     * using the reference source returned by {@link #getReference}. Tools should override
     * as appropriate. The default function is the identity function (ie., return exactly
     * the reference bases that span each read).
     *
     * @return reference window function used to initialize the reference source
     */
    public SerializableFunction<GATKRead, SimpleInterval> getReferenceWindowFunction() {
        return ReferenceWindowFunctions.IDENTITY_FUNCTION;
    }

    /**
     * Returns the "best available" sequence dictionary. This will be the reference sequence dictionary if
     * there is a reference, otherwise it will be the sequence dictionary constructed from the reads if
     * there are reads, otherwise it will be null.
     *
     * TODO: check interval file(s) as well for a sequence dictionary
     *
     * @return best available sequence dictionary given our inputs
     */
    public SAMSequenceDictionary getBestAvailableSequenceDictionary() {
        return hasReference() ? referenceDictionary : (hasReads() ? readsHeader.getSequenceDictionary() : null);
    }

    /**
     * @return sequence dictionary for the reference, or null if there is no reference available
     */
    public SAMSequenceDictionary getReferenceSequenceDictionary() {
        return referenceDictionary;
    }

    /**
     * @return header for the reads, or null if there are no reads available
     */
    public SAMFileHeader getHeaderForReads() {
        return readsHeader;
    }

    /**
     * Loads the reads into a {@link JavaRDD} using the intervals specified, and filters them using
     * the filter returned by {@link #makeReadFilter}.
     *
     * If no intervals were specified, returns all the reads (both mapped and unmapped).
     *
     * @return all reads from our reads input(s) as a {@link JavaRDD}, bounded by intervals if specified,
     *         and filtered using the filter from {@link #makeReadFilter}.
     */
    public JavaRDD<GATKRead> getReads() {
        final ReadFilter filter = makeReadFilter();
        return getUnfilteredReads().filter(read -> filter.test(read));
    }

    /**
     * Loads the reads into a {@link JavaRDD} using the intervals specified, and returns them
     * without applying any filtering.
     *
     * If no intervals were specified, returns all the reads (both mapped and unmapped).
     *
     * @return all reads from our reads input(s) as a {@link JavaRDD}, bounded by intervals if specified, and unfiltered.
     */
    public JavaRDD<GATKRead> getUnfilteredReads() {
        TraversalParameters traversalParameters;
        if ( intervalArgumentCollection.intervalsSpecified() ) {
            traversalParameters = intervalArgumentCollection.getTraversalParameters(getHeaderForReads().getSequenceDictionary());
        } else if ( hasIntervals() ) { // intervals may have been supplied by editIntervals
            traversalParameters = new TraversalParameters(getIntervals(), false);
        } else {
            traversalParameters = null; // no intervals were specified so return all reads (mapped and unmapped)
        }

        // TODO: This if statement is a temporary hack until #959 gets resolved.
        if (readInput.endsWith(".adam")) {
            try {
                return readsSource.getADAMReads(readInput, traversalParameters, getHeaderForReads());
            } catch (IOException e) {
                throw new UserException("Failed to read ADAM file " + readInput, e);
            }

        } else {
            if (hasCramInput() && !hasReference()){
                throw new UserException.MissingReference("A reference file is required when using CRAM files.");
            }
            final String refPath = hasReference() ?  referenceArguments.getReferenceFile().getAbsolutePath() : null;
            return readsSource.getParallelReads(readInput, refPath, traversalParameters, bamPartitionSplitSize);
        }
    }

    /**
     * Writes the reads from a {@link JavaRDD} to an output file.
     * @param ctx the JavaSparkContext to write.
     * @param outputFile path to the output bam/cram.
     * @param reads reads to write.
     */
    public void writeReads(final JavaSparkContext ctx, final String outputFile, JavaRDD<GATKRead> reads) {
        try {
            ReadsSparkSink.writeReads(ctx, outputFile,
                    hasReference() ? referenceArguments.getReferenceFile().getAbsolutePath() : null,
                    reads, readsHeader, shardedOutput ? ReadsWriteFormat.SHARDED : ReadsWriteFormat.SINGLE,
                    getRecommendedNumReducers());
        } catch (IOException e) {
            throw new UserException.CouldNotCreateOutputFile(outputFile,"writing failed", e);
        }
    }

    /**
     * Return the recommended number of reducers for a pipeline processing the reads. The number is
     * calculated by finding the total size (in bytes) of all the files in the input path, then
     * dividing by the target split size (determined by {@link #getTargetPartitionSize()}.
     * Subclasses that want to control the recommended number of reducers should typically override
     * {@link #getTargetPartitionSize()} rather than this method.
     * @return the recommended number of reducers
     */
    public int getRecommendedNumReducers() {
        if (numReducers != 0) {
            return numReducers;
        }
        return 1 + (int) (BucketUtils.dirSize(getReadSourceName()) / getTargetPartitionSize());
    }

    /**
     * Returns the size of each input partition (in bytes) that is used to determine the recommended number of reducers
     * for running a processing pipeline. The larger the number of reducers used, the smaller the amount of memory
     * each one needs.
     *
     * Defaults to 10MB, but subclasses can override to change the value. Memory intensive pipelines should decrease
     * the partition size, while pipelines with lighter memory requirements may increase the partition size.
     */
    public int getTargetPartitionSize() {
        return 10 * 1024 * 1024; // 10MB
    }

    /**
     * Helper method that simply returns a boolean regarding whether the input has CRAM files or not.
     */
    private boolean hasCramInput() {
        return readArguments.getReadFiles().stream().anyMatch(IOUtils::isCramFile);
    }

    /**
     * Returns a read filter (simple or composite) that can be applied to the reads returned from {@link #getReads}.
     * This implementation combines the default read filters for this tool (returned by {@link #getDefaultReadFilters}
     * along with any read filter command line directives specified by the user (such as enabling other filters or
     * disabling default filters); and returns a single composite filter resulting from the list by and'ing them together.
     *
     * NOTE: Most tools will not need to override the method, and should only do so in order to provide custom
     * behavior or processing of the final merged read filter. To change the default read filters used by the tool,
     * override {@link #getDefaultReadFilters} instead.
     *
     * Multiple filters can be composed by using {@link org.broadinstitute.hellbender.engine.filters.ReadFilter}
     * composition methods.
     */
    public ReadFilter makeReadFilter() {
        final GATKReadFilterPluginDescriptor readFilterPlugin =
                getCommandLineParser().getPluginDescriptor(GATKReadFilterPluginDescriptor.class);
        return readFilterPlugin.getMergedReadFilter(getHeaderForReads());
    }

    /**
     * Returns the default list of ReadFilters that are used for this tool. The filters returned
     * by this method are subject to selective enabling/disabling by the user via the command line. The
     * default implementation uses the {@link WellformedReadFilter} filter with all default options. Subclasses
     * can override to provide alternative filters.
     *
     * Note: this method is called before command line parsing begins, and thus before a SAMFileHeader is
     * available through {@link #getHeaderForReads}. The actual SAMFileHeader is propagated to the read filters
     * by {@link #makeReadFilter} after the filters have been merged with command line arguments.
     *
     * @return List of individual filters to be applied for this tool.
     */
    public List<ReadFilter> getDefaultReadFilters() {
        return Arrays.asList(new WellformedReadFilter());
    }

    /**
     * Returns the name of the source of reads data. It can be a file name or URL.
     */
    protected String getReadSourceName(){
        return readInput;
    }

    /**
     * @return our reference source, or null if no reference is present
     */
    public ReferenceMultiSource getReference() {
        return referenceSource;
    }

    /**
     * @return our intervals, or null if no intervals were specified
     */
    public List<SimpleInterval> getIntervals() {
        return intervals;
    }

    @Override
    protected void runPipeline( JavaSparkContext sparkContext ) {
        initializeToolInputs(sparkContext);
        validateToolInputs();
        runTool(sparkContext);
    }

    /**
     * Initialize standard tool inputs.
     */
    private void initializeToolInputs(final JavaSparkContext sparkContext) {
        initializeReference();
        initializeReads(sparkContext); // reference must be initialized before reads
        initializeFeatures();
        initializeIntervals();
    }

    /**
     * Initializes our reads source (but does not yet load the reads into a {@link JavaRDD}).
     * Does nothing if no reads inputs are present.
     */
    private void initializeReads(final JavaSparkContext sparkContext) {
        if ( readArguments.getReadFilesNames().isEmpty() ) {
            return;
        }

        if ( readArguments.getReadFilesNames().size() != 1 ) {
            throw new UserException("Sorry, we only support a single reads input for spark tools for now.");
        }

        readInput = readArguments.getReadFilesNames().get(0);
        readsSource = new ReadsSparkSource(sparkContext, readArguments.getReadValidationStringency());
        readsHeader = readsSource.getHeader(
                readInput,
                hasReference() ?  referenceArguments.getReferenceFile().getAbsolutePath() : null);
    }

    /**
     * Initializes our reference source. Does nothing if no reference was specified.
     */
    private void initializeReference() {
        final GCSOptions gcsOptions = getAuthenticatedGCSOptions(); // null if we have no api key
        final String referenceURL = referenceArguments.getReferenceFileName();
        if ( referenceURL != null ) {
            referenceSource = new ReferenceMultiSource(gcsOptions, referenceURL, getReferenceWindowFunction());
            referenceDictionary = referenceSource.getReferenceSequenceDictionary(readsHeader != null ? readsHeader.getSequenceDictionary() : null);
            if (referenceDictionary == null) {
                throw new UserException.MissingReferenceDictFile(referenceURL);
            }
        }
    }

    /**
     * Initialize our source of Feature data (or set it to null if no Feature argument(s) were provided).
     *
     * Package-private so that engine classes can access it, but concrete tool child classes cannot.
     * May be overridden by traversals that require custom initialization of Feature data sources.
     *
     * By default, this method initializes the FeatureManager to use the lookahead cache of {@link FeatureDataSource#DEFAULT_QUERY_LOOKAHEAD_BASES} bases.
     */
    void initializeFeatures() {
        features = new FeatureManager(this);
        if ( features.isEmpty() ) {  // No available sources of Features discovered for this tool
            features = null;
        }
    }

    /**
     * Loads our intervals using the best available sequence dictionary (as returned by {@link #getBestAvailableSequenceDictionary})
     * to parse/verify them. Does nothing if no intervals were specified.
     */
    private void initializeIntervals() {
        if ( intervalArgumentCollection.intervalsSpecified() ) {
            final SAMSequenceDictionary intervalDictionary = getBestAvailableSequenceDictionary();
            if ( intervalDictionary == null ) {
                throw new UserException("We require at least one input source that has a sequence dictionary (reference or reads) when intervals are specified");
            }

            intervals = intervalArgumentCollection.getIntervals(intervalDictionary);
        }
        intervals = editIntervals(intervals);
    }

    /**
     * Transform the intervals during loading.
     *
     * Developers can override this method to do custom interval handling during initialization of their GATKSparkTool
     *
     * @param rawIntervals Intervals specified on command line by user (-L).  Can be {@code null}
     * @return Transformed set of intervals.  Allowed to return non-null, if null was specified in the input.
     */
    protected List<SimpleInterval> editIntervals(final List<SimpleInterval> rawIntervals) {
        return rawIntervals;
    }

    /**
     * Validates standard tool inputs against each other.
     */
    private void validateToolInputs() {
        if ( ! disableSequenceDictionaryValidation ) {
            // Validate the reference sequence dictionary against the reads sequence dictionary, if both are present,
            // using standard GATK validation settings (requiring a common subset of equivalent contigs without respect
            // to ordering).
            // Check the reference dictionary against the reads dictionary
            if ( hasReference() && hasReads() ) {
                if ( hasCramInput() ) {
                    // Use stricter validation for CRAM vs. the reference
                    SequenceDictionaryUtils.validateCRAMDictionaryAgainstReference(referenceDictionary, readsHeader.getSequenceDictionary());
                }
                else {
                    // Use standard validation settings for non-CRAM reads input vs. the reference
                    SequenceDictionaryUtils.validateDictionaries("reference", referenceDictionary, "reads", readsHeader.getSequenceDictionary());
                }
            }
        }
    }

    /**
     * Runs the tool itself after initializing and validating inputs. Must be implemented by subclasses.
     *
     * @param ctx our Spark context
     */
    protected abstract void runTool( JavaSparkContext ctx );
}
