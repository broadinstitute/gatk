package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.utils.GenomeLoc;
import org.broadinstitute.hellbender.utils.GenomeLocParser;

import java.io.File;
import java.util.List;
import java.util.Optional;
import java.util.stream.StreamSupport;

import static org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary.*;

/**
 * A ReadWalker is a tool that processes a single read at a time from one or multiple sources of reads, with
 * optional contextual information from a reference and/or sets of variants.
 *
 * If multiple sources of reads are specified, they are merged together into a single sorted stream of reads.
 *
 * ReadWalker authors must implement the apply() method to process each read, and may optionally implement
 * onTraversalStart() and/or onTraversalDone(). See the PrintReadsWithReference walker for an example.
 */
public abstract class ReadWalker extends GATKTool {

    @Argument(fullName = StandardArgumentDefinitions.INPUT_LONG_NAME, shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME, doc = "One or more BAM/SAM/CRAM files containing reads", common = false, optional = false, minElements = 1)
    public List<File> READS_FILES;

    @Argument(fullName = StandardArgumentDefinitions.REFERENCE_LONG_NAME, shortName = StandardArgumentDefinitions.REFERENCE_SHORT_NAME, doc = "Reference sequence", common = false, optional = true)
    public File REFERENCE_FILE;

    private ReadsDataSource reads = null;
    private ReferenceDataSource reference = null;
    private FeatureManager features = null;

    /**
     * Create the reads and reference data sources.
     *
     * Marked final so that tool authors don't override it. Tool authors should override onTraversalStart() instead.
     */
    @Override
    protected final void onStartup() {
        super.onStartup();

        // Need to delay initialization of members until after the argument-parsing system has injected argument values
        reads = new ReadsDataSource(READS_FILES);
        reference = REFERENCE_FILE != null ? new ReferenceDataSource(REFERENCE_FILE) : null;

        features = new FeatureManager(this);
        if ( features.isEmpty() ) // No available sources of Features for this tool
            features = null;

        if(intervalArgumentCollection.intervalsSpecified()){
            reads.setIntervalsForTraversal(intervalArgumentCollection.getIntervals(getBestAvailableSequenceDictionary()));
        }
    }

    /**
     * Returns true if a reference was provided for this traversal, otherwise false
     *
     * @return true if a reference was provided for this traversal, otherwise false
     */
    public boolean referenceIsPresent() {
        return REFERENCE_FILE != null;
    }

    /**
     * Returns the SAM header for this reads traversal. Will be a merged header if there are multiple inputs for
     * the reads. If there is only a single input, returns its header directly.
     *
     * @return SAM header for this reads traversal
     */
    public SAMFileHeader getHeaderForReads() {
        return reads.getHeader();
    }

    /**
     * Returns the "best available" sequence dictionary. This will be the reference sequence dictionary if
     * there is a reference, otherwise it will be the sequence dictionary constructed from the reads.
     *
     * @return best available sequence dictionary given our inputs (never null)
     */
    public SAMSequenceDictionary getBestAvailableSequenceDictionary() {
        return reference != null ? reference.getSequenceDictionary() : reads.getSequenceDictionary();
    }

    /**
     * Implementation of read-based traversal.
     * Subclasses can override to provide own behavior but default implementation should be suitable for most uses.
     *
     * The default implementation creates filters using {@link #makeReadFilter}
     * and then iterates over all reads, applies the filter and hands the resulting reads to the {@link #apply}
     * function of the walker (long with additional contextual information, if present, such as reference bases).
     */
    @Override
    public void traverse() {
        final GenomeLocParser genomeLocParser = new GenomeLocParser(getBestAvailableSequenceDictionary());

        // Process each read in the input stream.
        // Supply reference bases spanning each read, if a reference is available.
        ReadFilter filter = makeReadFilter();
        StreamSupport.stream(reads.spliterator(), false)
                .filter(filter)
                .forEach(read -> {
                    final GenomeLoc readInterval = genomeLocParser.createGenomeLoc(read);
                    apply(read,
                          reference == null ? Optional.empty() : Optional.of(new ReferenceContext(reference, readInterval)),
                          features == null  ? Optional.empty() : Optional.of(new FeatureContext(features, readInterval)));
                });
    }

    /**
     * Returns the read filter (simple or composite) that will be applied to the reads before calling {@link #apply}.
     * The default implementation uses the {#link MalformedReadFilter} filter with all default options.
     * Default implementation of {@link #traverse()} calls this method once before iterating
     * over the reads and reuses the filter object to avoid object allocation. Nevertheless, keeping state in filter objects is strongly discouraged.
     *
     * Subclasses can extend to provide own filters (ie override and call super).
     * Multiple filters can be composed by using {@link org.broadinstitute.hellbender.engine.filters.ReadFilter} composition methods.
     */
    public ReadFilter makeReadFilter(){
          return ALLOW_ALL_READS;//TODO reanable MalformedReadFilter https://github.com/broadinstitute/hellbender/issues/180
    }

    /**
     * Process an individual read (with optional contextual information). Must be implemented by tool authors.
     * In general, tool authors should simply stream their output from apply(), and maintain as little internal state
     * as possible.
     *
     * TODO: Determine whether and to what degree the GATK engine should provide a reduce operation
     * TODO: to complement this operation. At a minimum, we should make apply() return a value to
     * TODO: discourage statefulness in walkers, but how this value should be handled is TBD.
     * @param read current read
     * @param referenceContext Reference bases spanning the current read (Optional.empty() if no reference was specified).
     *                         Can request extra bases of context around the current read's interval by invoking
     *                         {@link ReferenceContext#setWindow} on this object before calling {@link ReferenceContext#getBases}
     * @param featureContext Features spanning the current read (Optional.empty() if no Feature inputs were specified for this tool)
     */
    public abstract void apply( SAMRecord read, Optional<ReferenceContext> referenceContext, Optional<FeatureContext> featureContext );

    /**
     * Close the reads and reference data sources.
     *
     * Marked final so that tool authors don't override it. Tool authors should override onTraversalDone() instead.
     */
    @Override
    protected final void onShutdown() {
        super.onShutdown();

        if ( reads != null )
            reads.close();

        if ( reference != null )
            reference.close();

        if ( features != null )
            features.close();
    }
}
