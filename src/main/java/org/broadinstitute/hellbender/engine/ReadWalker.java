package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.filters.MalformedReadFilter;
import org.broadinstitute.hellbender.utils.GenomeLocParser;

import java.io.File;
import java.util.List;
import java.util.function.Predicate;
import java.util.stream.StreamSupport;

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

        if(intervalArgumentCollection.intervalsSpecified()){
            SAMSequenceDictionary sequenceDict = reference != null ? reference.getSequenceDictionary() : reads.getSequenceDictionary();
            reads.setIntervalsForTraversal(intervalArgumentCollection.getIntervals(sequenceDict));
        }
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
     * Returns the sequence dictionary for the reference, if a reference is present. If there is no
     * reference, returns null.
     *
     * @return sequence dictionary for the reference, or null if there is no reference
     */
    public SAMSequenceDictionary getReferenceDictionary() {
        return reference != null ? reference.getSequenceDictionary() : null;
    }

    /**
     * Implementation of read-based traversal.
     * Subclasses can override to provide own behavior but default implementation should be suitable for most uses.
     *
     * The default implementation creates filters using {@link #makeReadFilter}
     * and then iterates over all reads, applies the filter and hands the resulting reads to the {@link #apply}
     * function of the walker.
     */
    @Override
    public void traverse() {
        final GenomeLocParser genomeLocParser = reference != null ? new GenomeLocParser(reference.getSequenceDictionary()) : null;

        // Process each read in the input stream.
        // Supply reference bases spanning each read, if a reference is available.
        Predicate<SAMRecord> filter = makeReadFilter();
        StreamSupport.stream(reads.spliterator(), false)
                .filter(filter)
                .forEach(read -> {
                    final ReferenceContext refContext = reference == null ? null : new ReferenceContext(reference, genomeLocParser.createGenomeLoc(read));
                    apply(read, refContext);
                });
    }

    /**
     * Returns the read filter (simple or composite) that will be applied to the reads before calling {@link #apply}.
     * The default implementation uses the {#link MalformedReadFilter} filter with all default options.
     * Default implementation of {@link #traverse()} calls this method once before iterating
     * over the reads and reuses the filter object to avoid object allocation. Nevertheless, keeping state in filter objects is strongly discouraged.
     *
     * Subclasses can extend to provide own filters (ie override and call super).
     * Multiple filters can be composed by using {@link Predicate} composition methods.
     */
    public Predicate<SAMRecord> makeReadFilter(){
          return read -> true;//HACK to check
        //return new MalformedReadFilter();
    }

    /**
     * Process an individual read (with optional contextual information). Must be implemented by tool authors.
     * In general, tool authors should simply stream their output from apply(), and maintain as little internal state
     * as possible.
     *
     * TODO: Determine whether and to what degree the GATK engine should provide a reduce operation
     * TODO: to complement this operation. At a minimum, we should make apply() return a value to
     * TODO: discourage statefulness in walkers, but how this value should be handled is TBD.
     *
     * @param read current read
     * @param referenceContext Reference bases spanning the current read (null if no reference was specified).
     *                         Can request extra bases of context around the current read's interval by invoking
     *                         {@link ReferenceContext#setWindow} on this object before calling {@link ReferenceContext#getBases}
     */
    public abstract void apply(final SAMRecord read, final ReferenceContext referenceContext);

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
    }
}
