package org.broadinstitute.hellbender.engine.spark;

import com.google.cloud.dataflow.sdk.transforms.SerializableFunction;
import com.google.cloud.genomics.dataflow.utils.GCSOptions;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.IntervalArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.OptionalIntervalArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.OptionalReadInputArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.OptionalReferenceInputArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.ReadInputArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.ReferenceInputArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.RequiredIntervalArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.RequiredReadInputArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.RequiredReferenceInputArgumentCollection;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.engine.datasources.ReferenceWindowFunctions;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SequenceDictionaryUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.IOException;
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
    private static final long serialVersionUID = 1l;

    @ArgumentCollection
    public final ReferenceInputArgumentCollection referenceArguments = requiresReference() ? new RequiredReferenceInputArgumentCollection() :  new OptionalReferenceInputArgumentCollection();

    @ArgumentCollection
    public final ReadInputArgumentCollection readArguments = requiresReads() ? new RequiredReadInputArgumentCollection() : new OptionalReadInputArgumentCollection();

    @ArgumentCollection
    protected IntervalArgumentCollection intervalArgumentCollection = requiresIntervals() ? new RequiredIntervalArgumentCollection() : new OptionalIntervalArgumentCollection();

    @Argument(doc = "maximum number of bytes to read from a file into each partition of reads. " +
            "Setting this higher will result in fewer partitions. Note that this will not be equal to the size of the partition in memory",
            fullName = "bamPartitionSize", shortName = "bps", optional = true)
    protected long bamPartitionSplitSize = ReadsSparkSource.DEFAULT_SPLIT_SIZE;

    @Argument(fullName = "disableSequenceDictionaryValidation", shortName = "disableSequenceDictionaryValidation", doc = "If specified, do not check the sequence dictionaries from our inputs for compatibility. Use at your own risk!", optional = true)
    private boolean disableSequenceDictionaryValidation = false;

    private ReadsSparkSource readsSource;
    private SAMFileHeader readsHeader;
    private String readInput;
    private ReferenceMultiSource referenceSource;
    private SAMSequenceDictionary referenceDictionary;
    private List<SimpleInterval> intervals;

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
     * Loads the reads into a {@link JavaRDD} using the intervals specified. If no intervals
     * were specified, returns all the reads (both mapped and unmapped).
     *
     * @return all reads from our reads input(s) as a {@link JavaRDD}, bounded by intervals if specified.
     */
    public JavaRDD<GATKRead> getReads() {
        // TODO: This if statement is a temporary hack until #959 gets resolved.
        if (readInput.endsWith(".adam")) {
            try {
                return readsSource.getADAMReads(readInput, intervals, getHeaderForReads());
            } catch (IOException e) {
                throw new UserException("Failed to read ADAM file " + readInput);
            }

        } else {
            // If no intervals were specified (intervals == null), this will return all reads (mapped and unmapped)
            return readsSource.getParallelReads(readInput, intervals, bamPartitionSplitSize);
        }
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
    private void initializeToolInputs(JavaSparkContext sparkContext) {
        initializeReads(sparkContext);
        initializeReference();
        initializeIntervals();
    }

    /**
     * Initializes our reads source (but does not yet load the reads into a {@link JavaRDD}).
     * Does nothing if no reads inputs are present.
     */
    private void initializeReads(JavaSparkContext sparkContext) {
        if ( readArguments.getReadFilesNames().isEmpty() ) {
            return;
        }

        if ( readArguments.getReadFilesNames().size() != 1 ) {
            throw new UserException("Sorry, we only support a single reads input for spark tools for now.");
        }

        readInput = readArguments.getReadFilesNames().get(0);
        readsSource = new ReadsSparkSource(sparkContext);
        readsHeader = ReadsSparkSource.getHeader(sparkContext, readInput, getAuthHolder());
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
    }

    /**
     * Validates standard tool inputs against each other.
     */
    private void validateToolInputs() {
        if ( ! disableSequenceDictionaryValidation ) {
            if ( hasReference() && hasReads() ) {
                SequenceDictionaryUtils.validateDictionaries("reference", referenceDictionary, "reads", readsHeader.getSequenceDictionary());
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
