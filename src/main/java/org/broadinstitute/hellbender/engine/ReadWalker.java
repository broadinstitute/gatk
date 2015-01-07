package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.cmdline.Option;
import org.broadinstitute.hellbender.cmdline.StandardOptionDefinitions;
import org.broadinstitute.hellbender.utils.GenomeLocParser;

import java.io.File;
import java.util.List;

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

    @Option(fullName = StandardOptionDefinitions.INPUT_LONG_NAME, shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "One or more BAM/SAM/CRAM files containing reads", common = false, optional = false, overridable = true, minElements = 1)
    public List<File> READS_FILES;

    @Option(fullName = StandardOptionDefinitions.REFERENCE_LONG_NAME, shortName = StandardOptionDefinitions.REFERENCE_SHORT_NAME, doc = "Reference sequence", common = false, optional = true, overridable = true)
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
    }

    /**
     * Implementation of read-based traversal.
     */
    @Override
    public void traverse() {
        final GenomeLocParser genomeLocParser = reference != null ? new GenomeLocParser(reference.getSequenceDictionary()) : null;

        // Invoke apply() on each read in the input stream.
        // Supply reference bases spanning each read, if a reference is available.
        for ( SAMRecord read : reads ) {
            final ReferenceContext refContext = reference != null ? new ReferenceContext(reference, genomeLocParser.createGenomeLoc(read)) : null;
            apply(read, refContext);
        }
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
     * @param referenceContext reference bases spanning the current read (null if no reference was specified)
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
