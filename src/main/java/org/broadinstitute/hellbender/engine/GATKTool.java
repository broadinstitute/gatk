package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.Locatable;
import htsjdk.tribble.Feature;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFSimpleHeaderLine;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLinePluginDescriptor;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.GATKPlugin.GATKReadFilterPluginDescriptor;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.*;
import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.filters.WellformedReadFilter;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.transformers.ReadTransformer;
import org.broadinstitute.hellbender.utils.SequenceDictionaryUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.read.SAMFileGATKReadWriter;
import org.broadinstitute.hellbender.utils.reference.ReferenceUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.io.File;
import java.time.ZonedDateTime;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Stream;

/**
 * Base class for all GATK tools. Tool authors that wish to write a "GATK" tool but not use one of
 * the pre-packaged Walker traversals should feel free to extend this class directly. All other
 * GATK tools should extend one of the Walker classes instead.
 */
public abstract class GATKTool extends CommandLineProgram {

    @ArgumentCollection
    protected IntervalArgumentCollection intervalArgumentCollection = requiresIntervals() ? new RequiredIntervalArgumentCollection() : new OptionalIntervalArgumentCollection();

    @ArgumentCollection
    protected final ReadInputArgumentCollection readArguments = requiresReads() ? new RequiredReadInputArgumentCollection() : new OptionalReadInputArgumentCollection();

    @ArgumentCollection
    protected final ReferenceInputArgumentCollection referenceArguments = requiresReference() ? new RequiredReferenceInputArgumentCollection() :  new OptionalReferenceInputArgumentCollection();

    public static final String SEQUENCE_DICTIONARY_NAME = "sequenceDictionary";

    @Argument(fullName = SEQUENCE_DICTIONARY_NAME,
            shortName = SEQUENCE_DICTIONARY_NAME,
            doc = "Use the the given sequence dictionary for processing.", optional = true, common = true)
    private String masterSequenceDictionary = null;

    public static final String SECONDS_BETWEEN_PROGRESS_UPDATES_NAME = "secondsBetweenProgressUpdates";
    @Argument(fullName = SECONDS_BETWEEN_PROGRESS_UPDATES_NAME, shortName = SECONDS_BETWEEN_PROGRESS_UPDATES_NAME, doc = "Output traversal statistics every time this many seconds elapse", optional = true, common = true)
    private double secondsBetweenProgressUpdates = ProgressMeter.DEFAULT_SECONDS_BETWEEN_UPDATES;

    @Argument(fullName = StandardArgumentDefinitions.DISABLE_SEQUENCE_DICT_VALIDATION_NAME, shortName = StandardArgumentDefinitions.DISABLE_SEQUENCE_DICT_VALIDATION_NAME, doc = "If specified, do not check the sequence dictionaries from our inputs for compatibility. Use at your own risk!", optional = true, common = true)
    private boolean disableSequenceDictionaryValidation = false;

    @Argument(fullName=StandardArgumentDefinitions.CREATE_OUTPUT_BAM_INDEX_LONG_NAME,
            shortName=StandardArgumentDefinitions.CREATE_OUTPUT_BAM_INDEX_SHORT_NAME,
            doc = "If true, create a BAM/CRAM index when writing a coordinate-sorted BAM/CRAM file.", optional=true, common = true)
    public boolean createOutputBamIndex = true;

    @Argument(fullName=StandardArgumentDefinitions.CREATE_OUTPUT_BAM_MD5_LONG_NAME,
            shortName=StandardArgumentDefinitions.CREATE_OUTPUT_BAM_MD5_SHORT_NAME,
            doc = "If true, create a MD5 digest for any BAM/SAM/CRAM file created", optional=true, common = true)
    public boolean createOutputBamMD5 = false;

    @Argument(fullName=StandardArgumentDefinitions.CREATE_OUTPUT_VARIANT_INDEX_LONG_NAME,
            shortName=StandardArgumentDefinitions.CREATE_OUTPUT_VARIANT_INDEX_SHORT_NAME,
            doc = "If true, create a VCF index when writing a coordinate-sorted VCF file.", optional=true, common = true)
    public boolean createOutputVariantIndex = true;

    @Argument(fullName=StandardArgumentDefinitions.CREATE_OUTPUT_VARIANT_MD5_LONG_NAME,
            shortName=StandardArgumentDefinitions.CREATE_OUTPUT_VARIANT_MD5_SHORT_NAME,
            doc = "If true, create a a MD5 digest any VCF file created.", optional=true, common = true)
    public boolean createOutputVariantMD5 = false;

    @Argument(fullName= StandardArgumentDefinitions.LENIENT_LONG_NAME,
            shortName = StandardArgumentDefinitions.LENIENT_SHORT_NAME,
            doc = "Lenient processing of VCF files", common = true, optional = true)
    protected boolean lenientVCFProcessing = false;

    @Argument(fullName = StandardArgumentDefinitions.ADD_OUTPUT_SAM_PROGRAM_RECORD, shortName = StandardArgumentDefinitions.ADD_OUTPUT_SAM_PROGRAM_RECORD, doc = "If true, adds a PG tag to created SAM/BAM/CRAM files.", optional=true, common = true)
    public boolean addOutputSAMProgramRecord = true;

    @Argument(fullName = StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, shortName = StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, doc = "If true, adds a command line header line to created VCF files.", optional=true, common = true)
    public boolean addOutputVCFCommandLine = true;

    // default value of 40MB based on a test with CountReads (it's 5x faster than no prefetching)
    @Argument(fullName = StandardArgumentDefinitions.CLOUD_PREFETCH_BUFFER_LONG_NAME, shortName = StandardArgumentDefinitions.CLOUD_PREFETCH_BUFFER_SHORT_NAME, doc = "Size of the cloud-only prefetch buffer (in MB; 0 to disable).", optional=true)
    public int cloudPrefetchBuffer = getDefaultCloudPrefetchBufferSize();

    @Argument(fullName = StandardArgumentDefinitions.CLOUD_INDEX_PREFETCH_BUFFER_LONG_NAME, shortName = StandardArgumentDefinitions.CLOUD_INDEX_PREFETCH_BUFFER_SHORT_NAME, doc = "Size of the cloud-only prefetch buffer (in MB; 0 to disable). Defaults to cloudPrefetchBuffer if unset.", optional=true)
    public int cloudIndexPrefetchBuffer = getDefaultCloudIndexPrefetchBufferSize();

    @Argument(fullName = StandardArgumentDefinitions.DISABLE_BAM_INDEX_CACHING_LONG_NAME,
            shortName = StandardArgumentDefinitions.DISABLE_BAM_INDEX_CACHING_SHORT_NAME,
            doc = "If true, don't cache bam indexes, this will reduce memory requirements but may harm performance if many intervals are specified.  Caching is automatically disabled if there are no intervals specified.",
            optional = true)
    public boolean disableBamIndexCaching = false;

    /*
     * TODO: Feature arguments for the current tool are currently discovered through reflection via FeatureManager.
     * TODO: Perhaps we should eventually do the same auto-discovery for all input arguments (reads, reference, etc.)
     */

    /*
     * Engine-wide data structures. Package-private so that the various *Walker classes in the engine can access
     * them, but concrete tool implementations that extend from us cannot.
     */

    /**
     * Our source of reference data (null if no reference was provided)
     */
    ReferenceDataSource reference;

    /**
     * Our source of reads data (null if no source of reads was provided)
     */
    ReadsDataSource reads;

    /**
     * Our source of Feature data (null if no source of Features was provided)
     */
    FeatureManager features;

    /**
     * Intervals to be used for traversal (null if no intervals were provided).
     *
     * Walker base classes (ReadWalker, etc.) are responsible for hooking these intervals up to
     * their particular driving data source.
     */
    List<SimpleInterval> intervalsForTraversal;

    /**
     * Progress meter to print out traversal statistics. Subclasses must invoke
     * {@link ProgressMeter#update(Locatable)} after each record processed from
     * the primary input in their {@link #traverse} method.
     */
    protected ProgressMeter progressMeter;

    /**
     * Return the list of GATKCommandLinePluginDescriptors to be used for this tool.
     * Uses the read filter plugin.
     */
    @Override
    public List<? extends CommandLinePluginDescriptor<?>> getPluginDescriptors() {
        return Collections.singletonList(new GATKReadFilterPluginDescriptor(getDefaultReadFilters()));
    }

    /**
     * Returns the default list of ReadFilters that are used for this tool. The filters returned
     * by this method are subject to selective enabling/disabling by the user via the command line. The
     * default implementation returns an empty list. Subclasses can override to provide alternative filters.
     *
     * Note: this method is called before command line parsing begins, and thus before a SAMFileHeader is
     * available through {@link #getHeaderForReads}. The actual SAMFileHeader is propagated to the read filters
     * by {@link #makeReadFilter} after the filters have been merged with command line arguments.
     *
     * @return List of individual filters to be applied for this tool.
     */
    public List<ReadFilter> getDefaultReadFilters() {
        return Collections.singletonList(new WellformedReadFilter());
    }

    /**
     * Returns a read filter (simple or composite) that can be applied to reads. This implementation combines
     * the default read filters for this tool (returned by {@link #getDefaultReadFilters} along with any read filter
     * command line directives specified by the user (such as enabling other filters or disabling default filters);
     * wraps each filter in the resulting list with a CountingReadFilter; and returns a single composite filter
     * resulting from the list by and'ing them together.
     *
     * NOTE: Most tools will not need to override the method, and should only do so in order to provide custom
     * behavior or processing of the final merged read filter. To change the default read filters used by the tool,
     * override {@link #getDefaultReadFilters} instead.
     *
     * Implementations of {@link #traverse()} should call this method once before iterating over the reads, in order to
     * unnecessary avoid object allocation. Nevertheless, keeping state in filter objects is strongly discouraged.
     *
     * Multiple filters can be composed by using {@link org.broadinstitute.hellbender.engine.filters.ReadFilter}
     * composition methods.
     */
     public CountingReadFilter makeReadFilter(){
        final GATKReadFilterPluginDescriptor readFilterPlugin =
                getCommandLineParser().getPluginDescriptor(GATKReadFilterPluginDescriptor.class);
        return hasReads() ?
                readFilterPlugin.getMergedCountingReadFilter(getHeaderForReads()) :
                new CountingReadFilter(ReadFilterLibrary.ALLOW_ALL_READS);
    }

    /**
     * Returns the pre-filter read transformer (simple or composite) that will be applied to the reads before filtering.
     * The default implementation uses the {@link ReadTransformer#identity()}.
     * Default implementation of {@link #traverse()} calls this method once before iterating over the reads and reuses
     * the transformer object to avoid object allocation.
     *
     * Subclasses can extend to provide own transformers (ie override and call super).
     * Multiple transformers can be composed by using {@link ReadTransformer} composition methods.
     */
    public ReadTransformer makePreReadFilterTransformer() {
        return ReadTransformer.identity();
    }

    /**
     * Returns the post-filter read transformer (simple or composite) that will be applied to the reads after filtering.
     * The default implementation uses the {@link ReadTransformer#identity()}.
     * Default implementation of {@link #traverse()} calls this method once before iterating over the reads and reuses
     * the transformer object to avoid object allocation.
     *
     * Subclasses can extend to provide own transformers (ie override and call super).
     * Multiple transformers can be composed by using {@link ReadTransformer} composition methods.
     */
    public ReadTransformer makePostReadFilterTransformer(){
        return ReadTransformer.identity();
    }


    /**
     * Returns a stream over the reads, which are:
     *
     * 1. Transformed with {@link #makePreReadFilterTransformer()}.
     * 2. Filtered with {@code filter}.
     * 3. Transformed with {@link #makePostReadFilterTransformer()}.
     *
     * Note: the filter is passed to keep the state of {@link CountingReadFilter}, obtained with {@link #makeReadFilter()}.
     */
    protected Stream<GATKRead> getTransformedReadStream(final ReadFilter filter) {
        // if has reads, return an transformed/filtered/transformed stream
        if (hasReads()) {
            final ReadTransformer preTransformer = makePreReadFilterTransformer();
            final ReadTransformer postTransformer = makePostReadFilterTransformer();
            return Utils.stream(reads)
                    .map(preTransformer)
                    .filter(filter)
                    .map(postTransformer);
        }
        // returns an empty Stream if there are no reads
        return Stream.empty();
    }

    /**
     * @return Default size in MB of the cloud prefetch buffer. May be overridden by individual tools.
     *         The default implementation returns a value (40 MB) that is suitable for tools with a small
     *         number of large cloud inputs. Tools with large numbers of cloud inputs will likely want to
     *         override to specify a smaller size.
     */
    public int getDefaultCloudPrefetchBufferSize() {
        return 40;
    }

    /**
     * @return Default size in MB of the cloud index prefetch buffer. May be overridden by individual tools.
     *         A return value of -1 means to use the same value as returned by {@link #getDefaultCloudPrefetchBufferSize()}.
     *         The default implementation returns -1.
     */
    public int getDefaultCloudIndexPrefetchBufferSize() {
        return -1;
    }

    /**
     * @return String label to use for records in progress meter output. Defaults to {@link ProgressMeter#DEFAULT_RECORD_LABEL},
     *         but tools may override to provide a more appropriate label (like "reads" or "regions")
     */
    public String getProgressMeterRecordLabel() { return ProgressMeter.DEFAULT_RECORD_LABEL; }

    /**
     * Initialize our source of reference data (or set it to null if no reference argument was provided).
     *
     * Package-private so that engine classes can access it, but concrete tool child classes cannot.
     * May be overridden by traversals that require custom initialization of the reference data source.
     */
    void initializeReference() {
        reference = referenceArguments.getReferenceFile() != null ? ReferenceDataSource.of(referenceArguments.getReferenceFile()) : null;
    }

    /**
     * Initialize our source of reads data (or set it to null if no reads argument(s) were provided).
     *
     * Package-private so that engine classes can access it, but concrete tool child classes cannot.
     * May be overridden by traversals that require custom initialization of the reads data source.
     */
    void initializeReads() {
        if (! readArguments.getReadFiles().isEmpty()) {
            SamReaderFactory factory = SamReaderFactory.makeDefault().validationStringency(readArguments.getReadValidationStringency());
            if (hasReference()) { // pass in reference if available, because CRAM files need it
                factory = factory.referenceSequence(referenceArguments.getReferenceFile());
            }
            else if (hasCramInput()) {
                throw new UserException.MissingReference("A reference file is required when using CRAM files.");
            }

            if(bamIndexCachingShouldBeEnabled()) {
                factory = factory.enable(SamReaderFactory.Option.CACHE_FILE_BASED_INDEXES);
            }

            reads = new ReadsDataSource(readArguments.getReadPaths(), readArguments.getReadIndexPaths(), factory, cloudPrefetchBuffer,
                (cloudIndexPrefetchBuffer < 0 ? cloudPrefetchBuffer : cloudIndexPrefetchBuffer));
        }
        else {
            reads = null;
        }
    }


    private boolean bamIndexCachingShouldBeEnabled() {
        return intervalArgumentCollection.intervalsSpecified() && !disableBamIndexCaching;
    }

    /**
     * Helper method that simply returns a boolean regarding whether the input has CRAM files or not.
     */
    private boolean hasCramInput() {
        return readArguments.getReadFiles().stream().anyMatch(IOUtils::isCramFile);
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
        features = new FeatureManager(this, FeatureDataSource.DEFAULT_QUERY_LOOKAHEAD_BASES, cloudPrefetchBuffer, cloudIndexPrefetchBuffer,
                                      referenceArguments.getReferencePath());
        if ( features.isEmpty() ) {  // No available sources of Features discovered for this tool
            features = null;
        }
    }

    /**
     * Initialize our intervals for traversal.
     *
     * If intervals were specified, requires that another data source (reads or reference) be present
     * to supply a sequence dictionary. This requirement may be removed in the future.
     *
     * Must be called after other data sources have been initialized because of the sequence dictionary
     * requirement.
     *
     * Package-private so that engine classes can access it, but concrete tool child classes cannot.
     * May be overridden by traversals that require custom initialization of intervals.
     */
    void initializeIntervals() {
        if ( intervalArgumentCollection.intervalsSpecified() ) {
            final SAMSequenceDictionary sequenceDictionary = getBestAvailableSequenceDictionary();
            if ( sequenceDictionary == null ) {
                throw new UserException("We require a sequence dictionary from a reference, a source of reads, or a source of variants to process intervals.  " +
                        "Since reference and reads files generally contain sequence dictionaries, this error most commonly occurs " +
                        "for VariantWalkers that do not require a reference or reads.  You can fix the problem by passing a reference file with a sequence dictionary " +
                        "via the -R argument or you can run the tool UpdateVCFSequenceDictionary on your vcf.");
            }

            intervalsForTraversal = intervalArgumentCollection.getIntervals(sequenceDictionary);
        }
    }

    /**
     * Is a source of reference data available?
     *
     * @return true if a reference is available, otherwise false
     */
    public final boolean hasReference() {
        return reference != null;
    }

    /**
     * Are sources of reads available?
     *
     * @return true if reads are available, otherwise false
     */
    public final boolean hasReads() {
        return reads != null;
    }

    /**
     * Are sources of Features available?
     *
     * @return true if sources of Features are available, otherwise false
     */
    public final boolean hasFeatures() {
        return features != null;
    }

    /**
     * Are sources of intervals available?
     *
     * @return true if intervals are available, otherwise false
     */
    public final boolean hasIntervals() {
        return intervalsForTraversal != null;
    }

    /**
     * Does this tool require reference data? Traversals types and/or tools that do should override to return true.
     *
     * @return true if this tool requires a reference, otherwise false
     */
    public boolean requiresReference() {
        return false;
    }

    /**
     * Does this tool require features? Traversals types and/or tools that do should override to return true.
     *
     * @return true if this tool requires reads, otherwise false
     */
    public boolean requiresFeatures() {
        return false;
    }

    /**
     * Does this tool require reads? Traversals types and/or tools that do should override to return true.
     *
     * @return true if this tool requires reads, otherwise false
     */
    public boolean requiresReads() {
        return false;
    }

    /**
     * Does this tool require intervals? Traversals types and/or tools that do should override to return true.
     *
     * @return true if this tool requires intervals, otherwise false
     */
    public boolean requiresIntervals() {
        return false;
    }

    /**
     * Returns the reference sequence dictionary if there is a reference (hasReference() == true), otherwise null.
     * @return reference sequence dictionary if any, or null
     */
    public final SAMSequenceDictionary getReferenceDictionary() {
        return reference != null ? reference.getSequenceDictionary() : null;
    }

    /**
     * Returns the "best available" sequence dictionary or {@code null} if there is no single best dictionary.
     *
     * The algorithm for selecting the best dictionary is as follows:
     * 1) If a master sequence dictionary was specified, use that dictionary
     * 2) if there is a reference, then the best dictionary is the reference sequence dictionary
     * 3) Otherwise, if there are reads, then the best dictionary is the sequence dictionary constructed from the reads.
     * 4) Otherwise, if there are features and the feature data source has only one dictionary, then that one is the best dictionary.
     * 5) Otherwise, the result is {@code null}.
     *
     * TODO: check interval file(s) as well for a sequence dictionary
     *
     * Subclasses may override if they prefer a different algorithm.
     *
     * @return best available sequence dictionary given our inputs or {@code null} if no one dictionary is the best one.
     */
    public SAMSequenceDictionary getBestAvailableSequenceDictionary() {
        if (masterSequenceDictionary != null) {
            return ReferenceUtils.loadFastaDictionary(new File(masterSequenceDictionary));
        } else if (hasReference()){
            return reference.getSequenceDictionary();
        } else if (hasReads()){
            return reads.getSequenceDictionary();
        } else if (hasFeatures()){
            final List<SAMSequenceDictionary> dictionaries = features.getVariantSequenceDictionaries();
            //If there is just one, it clearly is the best. Otherwise, noone is best.
            if (dictionaries.size() == 1){
                return dictionaries.get(0);
            }
        }
        return null;
    }

    /**
     * Returns the SAM header for the reads data source. Will be a merged header if there are multiple inputs for
     * the reads. If there is only a single input, returns its header directly. Null if there are no reads.
     *
     * @return SAM header for our source of reads (null if there are no reads)
     */
    public final SAMFileHeader getHeaderForReads() {
        return hasReads() ? reads.getHeader() : null;
    }

    /**
     * Returns the header for the specified source of Features
     * @param featureDescriptor FeatureInput whose header to retrieve
     * @param <T> type of Feature in our FeatureInput
     * @return header for the provided FeatureInput (null if we have no sources of Features)
     */
    public final <T extends Feature> Object getHeaderForFeatures( final FeatureInput<T> featureDescriptor ) {
        return hasFeatures() ? features.getHeader(featureDescriptor) : null;
    }

    /**
     * Initialize our data sources, make sure that all tool requirements for input data have been satisfied
     * and start the progress meter.
     *
     * The data sources are initialized in the following order:
     *   initializeReference
     *   initializeReads (must be initialized after reference for CRAM files)
     *   initializeFeatures
     *   initializeIntervals (must be initialized after all other sources, because they require a sequence dictionary)
     *
     * Then, data sources are checked by calls to:
     *    validateSequenceDictionaries (unless disabled by disableSequenceDictionaryValidation)
     *    checkToolRequirements
     *
     * Subclasses may extend (ie, they should call super and add functionality).
     */
    @Override
    protected void onStartup() {
        super.onStartup();

        initializeReference();

        initializeReads(); // Must be initialized after reference, in case we are dealing with CRAM and a reference is required

        initializeFeatures();

        initializeIntervals(); // Must be initialized after reference, reads and features, since intervals currently require a sequence dictionary from another data source

        if ( ! disableSequenceDictionaryValidation ) {
            validateSequenceDictionaries();
        }

        checkToolRequirements();

        progressMeter = new ProgressMeter(secondsBetweenProgressUpdates);
        progressMeter.setRecordLabel(getProgressMeterRecordLabel());
    }

    /**
     * Validates all sequence dictionaries by checking them against each other.
     *
     * Currently, the reads dictionary is checked against the reference dictionary,
     * and both are checked against any variant dictionaries. No interval dictionaries are
     * checked at this time.
     *
     * The reference dictionary is required to be a superset of the reads dictionary if there is at least
     * one cram input; otherwise, all dictionaries are required to have a common subset of equivalent contigs,
     * without respect to relative or absolute ordering (GATKTools use contig names rather than indices,
     * and so are not sensitive to contig ordering issues).
     */
    private void validateSequenceDictionaries() {
        final SAMSequenceDictionary masterSequenceDict = (masterSequenceDictionary != null) ?
                ReferenceUtils.loadFastaDictionary(new File(masterSequenceDictionary)) : null;

        final SAMSequenceDictionary refDict = hasReference() ? reference.getSequenceDictionary() : null;
        final SAMSequenceDictionary readDict = hasReads() ? reads.getSequenceDictionary() : null;
        final List<SAMSequenceDictionary> featureDicts = hasFeatures() ? features.getAllSequenceDictionaries() : Collections.emptyList();

        // Check the reference dictionary against the reads dictionary
        if ( hasReference() && hasReads() ) {
            if ( hasCramInput() ) {
                // Use stricter validation for CRAM vs. the reference
                SequenceDictionaryUtils.validateCRAMDictionaryAgainstReference(refDict, readDict);
            }
            else {
                // Use standard validation settings for non-CRAM reads input vs. the reference
                SequenceDictionaryUtils.validateDictionaries("reference", refDict, "reads", readDict);
            }
        }

        // Check all Feature dictionaries against the reference and/or reads dictionaries
        // TODO: pass file names associated with each sequence dictionary into validateDictionaries(); issue #660
        for ( final SAMSequenceDictionary featureDict : featureDicts ) {
            if (hasReference()){
                SequenceDictionaryUtils.validateDictionaries("reference", refDict, "features", featureDict);
            }
            if (hasReads()) {
                SequenceDictionaryUtils.validateDictionaries("reads", readDict, "features", featureDict);
            }
        }

        // Check the master dictionary against the reference / reads / features
        if (masterSequenceDictionary != null) {

            // Check against the reads
            if (hasReads()) {
                SequenceDictionaryUtils.validateDictionaries("sequence", masterSequenceDict, "reads", readDict);
            }

            // Check against the reference
            if (hasReference()) {
                SequenceDictionaryUtils.validateDictionaries("sequence", masterSequenceDict, "reads", refDict);
            }

            // Check against the features
            for (final SAMSequenceDictionary featureDict : featureDicts) {
                SequenceDictionaryUtils.validateDictionaries("sequence", masterSequenceDict, "features", featureDict);
            }

            // When overriding the master sequence dictionary and working with CRAM files, we need to make sure that the
            // reference we're using is completely contained by the sequence dictionary (i.e. the sequence dictionary is
            // equal to the reference OR the sequence dictionary is a superset of the reference).
            // This is accomplished by the call to validateCRAMDictionaryAgainstReference.
            if ( hasCramInput() ) {
                SequenceDictionaryUtils.validateCRAMDictionaryAgainstReference(masterSequenceDict, refDict);
            }

        }

        // we don't currently have access to seqdicts from intervals
        //if (hasIntervals()) {}
    }

    /**
     * Check that all tool requirements for input data have been satisfied.
     *
     * Must be called after data source initialization.
     */
    private void checkToolRequirements() {

        if ( requiresFeatures() && ! hasFeatures() ) {
            throw new UserException("Tool " + getClass().getSimpleName() + " requires features, but none were provided");
        }

    }

    /*
     * Create a common SAMFileWriter using the reference and read header for this tool.
     *
     * @param outputFile    - if this file has a .cram extension then a reference is required. Can not be null.
     * @param preSorted     - if true then the records must already be sorted to match the header sort order
     *
     * @throws UserException if outputFile ends with ".cram" and no reference is provided
     * @return SAMFileWriter
     */
    public final SAMFileGATKReadWriter createSAMWriter(final File outputFile, final boolean preSorted) {
        if (!hasReference() && IOUtils.isCramFile(outputFile)) {
            throw new UserException.MissingReference("A reference file is required for writing CRAM files");
        }

        return new SAMFileGATKReadWriter(
                        ReadUtils.createCommonSAMWriter(
                                outputFile,
                                referenceArguments.getReferenceFile(),
                                getHeaderForSAMWriter(),
                                preSorted,
                                createOutputBamIndex,
                                createOutputBamMD5
                        )
        );
    }

    /**
     * Creates a VariantContextWriter whose outputFile type is determined by
     * the outFile's extension, using the best available sequence dictionary for
     * this tool, and default index, leniency and md5 generation settings.
     *
     * @param outFile output File for this writer. May not be null.
     * @returns VariantContextWriter must be closed by the caller
     */
    public VariantContextWriter createVCFWriter(final File outFile) {
        Utils.nonNull(outFile);

        final SAMSequenceDictionary sequenceDictionary = getBestAvailableSequenceDictionary();

        List<Options> options = new ArrayList<>();
        if (lenientVCFProcessing) {
            options.add(Options.ALLOW_MISSING_FIELDS_IN_HEADER);
        }

        if (createOutputVariantIndex) {
            if (null == sequenceDictionary) {
                logger.warn("An variant index will not be created - a sequence dictionary is required to create an output index");
                // fall through and create without an index
            } else {
                options.add(Options.INDEX_ON_THE_FLY);
            }
        }

        return GATKVariantContextUtils.createVCFWriter(
                outFile,
                sequenceDictionary,
                createOutputVariantMD5,
                options.toArray(new Options[options.size()]));
    }

    /**
     * Returns the SAM header suitable for writing SAM/BAM/CRAM files produced by this tool.
     *
     * The default implementation calls {@link #getHeaderForReads} (and makes an empty header if that call returns null)
     * and optionally adds program tag to the header with a program version {@link #getVersion()}, program name {@link #getToolName()}
     * and command line {@link #getCommandLine()}.
     *
     * Subclasses may override.
     *
     * @return SAM header for the SAM writer with (optionally, if {@link #addOutputSAMProgramRecord} is true) program record appropriately.
     */
    protected SAMFileHeader getHeaderForSAMWriter(){
        final SAMFileHeader header = getHeaderForReads() == null ? new SAMFileHeader(): getHeaderForReads();
        if (addOutputSAMProgramRecord) {
            final SAMProgramRecord programRecord = new SAMProgramRecord(createProgramGroupID(header));
            programRecord.setProgramVersion(getVersion());
            programRecord.setCommandLine(getCommandLine());
            programRecord.setProgramName(getToolName());
            header.addProgramRecord(programRecord);
        }
        return header;
    }

    /**
     * @return If addOutputVCFCommandLine is true, a set of VCF header lines containing the tool name, version,
     * date and command line, otherwise an empty set.
     */
    protected Set<VCFHeaderLine> getDefaultToolVCFHeaderLines() {
        final Set<VCFHeaderLine> gatkToolHeaderLines = new HashSet<>();
        if (addOutputVCFCommandLine) {
            final Map<String, String> simpleHeaderLineMap = new HashMap<>(4);
            simpleHeaderLineMap.put("ID", this.getClass().getSimpleName());
            simpleHeaderLineMap.put("Version", getVersion());
            simpleHeaderLineMap.put("Date", Utils.getDateTimeForDisplay((ZonedDateTime.now())));
            simpleHeaderLineMap.put("CommandLine", getCommandLine());
            gatkToolHeaderLines.add(new VCFHeaderLine("source", this.getClass().getSimpleName()));
            gatkToolHeaderLines.add(new VCFSimpleHeaderLine(String.format("%sCommandLine", getToolkitName()), simpleHeaderLineMap));
        }
        return gatkToolHeaderLines;
    }

    /**
     * Returns the program group ID that will be used in the SAM writer.
     * Starts with {@link #getToolName} and looks for the first available ID by appending consecutive integers.
     */
    private String createProgramGroupID(final SAMFileHeader header) {
        final String toolName = getToolName();

        String pgID = toolName;
        SAMProgramRecord record = header.getProgramRecord(pgID);
        int count = 1;
        while (record != null){
            pgID = toolName + "." + String.valueOf(count++);
            record = header.getProgramRecord(pgID);
        }
        return pgID;
    }

    /**
     * @return The name of toolkit for this tool. Subclasses may override to provide a custom toolkit name.
     */
    protected String getToolkitName() { return "GATK"; }

    /**
     * Returns the name of this tool.
     * The default implementation returns the result of calling {@link# getToolkitName} followed by the simple
     * name of the class. Subclasses may override.
     */
    public String getToolName() {
        return String.format("%s %s", getToolkitName(), getClass().getSimpleName());
    }

    /**
     * Close all data sources on shutdown.
     * Subclasses may extend (ie they must call super).
     */
    @Override
    protected void onShutdown() {
        super.onShutdown();

        if ( hasReference() ) {
            reference.close();
        }

        if ( hasReads() ) {
            reads.close();
        }

        if ( hasFeatures() ) {
            features.close();
        }
    }

    /**
     * Operations performed just prior to the start of traversal. Should be overridden by tool authors
     * who need to process arguments local to their tool or perform other kinds of local initialization.
     *
     * Default implementation does nothing.
     */
    public void onTraversalStart() {}

    /**
     * A complete traversal from start to finish. Tool authors who wish to "roll their own" traversal
     * from scratch can extend this class directly and implement this method. Walker authors should
     * instead extend a Walker class and implement the Walker-appropriate apply() method, since the
     * Walker base classes implement the various kinds of traversals for you.
     */
    public abstract void traverse();

    /**
     * Operations performed immediately after a successful traversal (ie when no uncaught exceptions were thrown during the traversal).
     * Should be overridden by tool authors who need to close local resources, etc., after traversal.
     * Also allows tools to return a value representing the traversal result, which is printed by the engine.
     *
     * Default implementation does nothing and returns null.
     *
     * @return Object representing the traversal result, or null if a tool does not return a value
     */
    public Object onTraversalSuccess() { return null; }

    @Override
    protected final Object doWork() {
        try {
            onTraversalStart();
            progressMeter.start();
            traverse();
            progressMeter.stop();
            return onTraversalSuccess();
        } finally {
            closeTool();
        }
    }

    /**
     * This method is called by the GATK framework at the end of the {@link #doWork} template method.
     * It is called regardless of whether the {@link #traverse} has succeeded or not.
     * It is called <em>after</em> the {@link #onTraversalSuccess} has completed (successfully or not)
     * but before the {@link #doWork} method returns.
     *
     * In other words, on successful runs both {@link #onTraversalSuccess} and {@link #closeTool} will be called (in this order) while
     * on failed runs (when {@link #traverse} causes an exception), only {@link #closeTool} will be called.
     *
     * The default implementation does nothing.
     * Subclasses should override this method to close any resources that must be closed regardless of the success of traversal.
     */
    public void closeTool(){
    }
}
