package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.Locatable;
import htsjdk.tribble.Feature;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.argumentcollections.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SequenceDictionaryUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.read.SAMFileGATKReadWriter;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

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

    @Argument(fullName = "secondsBetweenProgressUpdates", shortName = "secondsBetweenProgressUpdates", doc = "Output traversal statistics every time this many seconds elapse", optional = true)
    private double secondsBetweenProgressUpdates = ProgressMeter.DEFAULT_SECONDS_BETWEEN_UPDATES;

    @Argument(fullName = "disableSequenceDictionaryValidation", shortName = "disableSequenceDictionaryValidation", doc = "If specified, do not check the sequence dictionaries from our inputs for compatibility. Use at your own risk!", optional = true)
    private boolean disableSequenceDictionaryValidation = false;

    @Argument(fullName="createOutputBamIndex", shortName="createOutputBamIndex", doc = "If true, create a BAM/CRAM index when writing a coordinate-sorted BAM/CRAM file", optional=true)
    public boolean createOutputBamIndex = true;

    @Argument(fullName="createOutputBamMD5", shortName="createOutputBamMD5", doc = "If true, create a MD5 digest for any BAM/SAM/CRAM file created", optional=true)
    public boolean createOutputBamMD5 = false;

    @Argument(fullName="addOutputSAMProgramRecord", shortName="addOutputSAMProgramRecord", doc = "If true, adds a PG tag to created SAM/BAM/CRAM files.", optional=true)
    public boolean addOutputSAMProgramRecord = true;

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
    ProgressMeter progressMeter;

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
            reads = new ReadsDataSource(readArguments.getReadFiles(), factory);
        }
        else {
            reads = null;
        }
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
        features = new FeatureManager(this);
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
                throw new UserException("We currently require a sequence dictionary (from a reference or source of reads) " +
                                        "to process intervals. This restriction may be removed in the future.");
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
     * Returns the "best available" sequence dictionary. This will be the reference sequence dictionary if
     * there is a reference, otherwise it will be the sequence dictionary constructed from the reads if
     * there are reads, otherwise it will be null.
     *
     * TODO: check interval file(s) as well for a sequence dictionary
     *
     * @return best available sequence dictionary given our inputs
     */
    public final SAMSequenceDictionary getBestAvailableSequenceDictionary() {
        return reference != null ? reference.getSequenceDictionary() : (reads != null ? reads.getSequenceDictionary() : null);
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
     *   initializeReads
     *   initializeIntervals
     *   initializeFeatures
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

        initializeIntervals(); // Must be initialized after reference and reads, since intervals currently require a sequence dictionary from another data source

        initializeFeatures();

        if ( ! disableSequenceDictionaryValidation ) {
            validateSequenceDictionaries();
        }

        checkToolRequirements();

        progressMeter = new ProgressMeter(secondsBetweenProgressUpdates);
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
        final SAMSequenceDictionary refDict = hasReference() ? reference.getSequenceDictionary() : null;
        final SAMSequenceDictionary readDict = hasReads() ? reads.getSequenceDictionary() : null;
        List<SAMSequenceDictionary> variantDicts = new ArrayList<>();
        if (hasFeatures()){
            List<VCFHeader> variantHeaders = features.getAllVariantHeaders();
            for (VCFHeader header : variantHeaders) {
                SAMSequenceDictionary headerDict = header.getSequenceDictionary();
                if (headerDict != null) {
                    variantDicts.add(headerDict);
                }
            }
        }

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

        // Check all variants dictionaries against the reference and/or reads dictionaries
        // TODO: pass file names associated with each sequence dictionary into validateDictionaries(); issue #660
        for ( SAMSequenceDictionary variantsDict : variantDicts ) {
            if (hasReference()){
                SequenceDictionaryUtils.validateDictionaries("reference", refDict, "variants", variantsDict);
            }
            if (hasReads()) {
                SequenceDictionaryUtils.validateDictionaries("reads", readDict, "variants", variantsDict);
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
     * Returns the name of this GATK tool.
     * The default implementation return the the string "GATK " followed by the simple name of the class.
     * Subclasses may override.
     */
    public String getToolName() {
        return "GATK " + getClass().getSimpleName();
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
        onTraversalStart();
        progressMeter.start();
        traverse();
        progressMeter.stop();
        return onTraversalSuccess();
    }
}
