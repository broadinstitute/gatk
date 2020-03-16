package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMTag;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.filters.WellformedReadFilter;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.consensus.MoleculeID;
import org.broadinstitute.hellbender.tools.walkers.consensus.ReadsWithSameUMI;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.ArrayList;
import java.util.List;

/**
 * A walker that processes duplicate reads that share the same Unique molecule Identifier (UMI) as a single unit.
 *
 * This tool assumes that the input bam has been sorted by UMI (the {@link SAMTag.MI} tag to be specific) with FGBio GroupReadsByUmi:
 * http://fulcrumgenomics.github.io/fgbio/tools/latest/GroupReadsByUmi.html
 */
public abstract class DuplicateSetWalker extends ReadWalker {
    public static final String MIN_REQUIRED_READS_NAME = "min-reads";
    public static final String MIN_REQUIRED_READS_PER_STRAND_NAME = "min-per-strand-reads";

    private static final int DEFAULT_MINIMUM_READS_PER_SET = 1;
    private static final int DEFAULT_MINIMUM_READS_PER_STRAND = 0;

    @Argument(fullName = MIN_REQUIRED_READS_NAME, doc = "The mininum total number of reads required in the set", optional = true, minValue = 0)
    private int minimumRequiredReadsPerUMI = DEFAULT_MINIMUM_READS_PER_SET;

    // The user may choose to only keep read sets containing both strands (i.d. duplex evidence) by setting this argument to a positive number
    @Argument(fullName = MIN_REQUIRED_READS_PER_STRAND_NAME, doc = "The mininum total number of reads in each strand", optional = true, minValue = 0)
    private int minimumRequiredReadsPerStrand = DEFAULT_MINIMUM_READS_PER_STRAND;

    protected ReadsWithSameUMI currentReadsWithSameUMI = null;

    @Override
    public final void traverse(){
        super.traverse();
        processLastReadSet();
    }

    /***
     * FGBio GroupByUMI returns reads sorted by molecule ID: For example, the input bam may look like
     * read1: ... MI:Z:0/A ...
     * read2: ... MI:Z:0/A ...
     * read3: ... MI:Z:0/B ...
     * read4: ... MI:Z:0/B ...
     * read5: ... MI:Z:1/A ...
     * read6: ... MI:Z:1/B ...
     * read7: ... MI:Z:1/B ...
     *
     * Thus it's sufficient to go through the reads in order and collect them in a list until
     * we encounter the next molecule ID, at which point we pass the list to the {@code apply} method,
     * process the set based on the child class's implementation of the method, and clear the {@code currentDuplicateSet} variable and start collecting reads again.
     *
     * Notice there are two apply() methods in this class:
     * This apply() inherited from ReadWalker is marked final to discourage subclassing.
     * A subclass must override the other apply() method that takes in the DuplicateSet.
     */
    @Override
    public final void apply(GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext) {
        if (currentReadsWithSameUMI == null){ // evaluates to true for the very first read
            currentReadsWithSameUMI = new ReadsWithSameUMI(read);
            return;
        }

        final int readMoleculeNumber = MoleculeID.getMoleculeNumberOfRead(read);
        final int duplicateSetMoleculeNumber = currentReadsWithSameUMI.getMoleculeNumber();

        // If the incoming read has the molecule id less than that of the currentDuplicateSet,
        // the input bam is not sorted properly by the MI tag
        if (duplicateSetMoleculeNumber > readMoleculeNumber){
            throw new UserException(String.format("The input bam must be sorted by the molecule ID (%s) tag.", SAMTag.MI.name()));
        }

        if (duplicateSetMoleculeNumber < readMoleculeNumber) {
            // The incoming read's molecule ID is greater than that of the current set, meaning we've reached the end of the current set.
            if (rejectSet(currentReadsWithSameUMI)){
                currentReadsWithSameUMI = new ReadsWithSameUMI(read);
                return;
            }

            // Call the apply() method to process the current set and start a new set.
            apply(currentReadsWithSameUMI,
                    new ReferenceContext(reference, currentReadsWithSameUMI.getInterval()), // Will create an empty ReferenceContext if reference or readInterval == null
                    new FeatureContext(features, currentReadsWithSameUMI.getInterval()));
            currentReadsWithSameUMI = new ReadsWithSameUMI(read);
            return;
        }

        // The incoming read has the same UMI as the current set; simply add to the current set
        currentReadsWithSameUMI.addRead(read);
    }

    /**
     * A subclass must specify how to process the duplicate sets by overriding this method.
     *
     * @param readsWithSameUMI A set of reads with the matching UMIs with the same fragment start and end
     * @param referenceContext A reference context object over the intervals determined by the duplicate set.
     * @param featureContext Entries from a secondary feature file (e.g. vcf) if provided
     *
     */
    public abstract void apply(ReadsWithSameUMI readsWithSameUMI, ReferenceContext referenceContext, FeatureContext featureContext );

    private void processLastReadSet(){
        if (currentReadsWithSameUMI.getReads().size() > 0){
            apply(currentReadsWithSameUMI,
                    new ReferenceContext(reference, currentReadsWithSameUMI.getInterval()),
                    new FeatureContext(features, currentReadsWithSameUMI.getInterval()));
        }
    }

    /**
     * Returns true for duplicate sets that does not meet required criteria for further processing.
     * We encourage the user to override this method to meet their needs.
     */
    protected boolean rejectSet(final ReadsWithSameUMI readsWithSameUMI){
        // Check that the set contains the minimum required number of reads in each strand
        final Pair<Integer, Integer> strandCounts = MoleculeID.countStrands(readsWithSameUMI.getReads());
        if (Math.min(strandCounts.getLeft(), strandCounts.getRight()) < minimumRequiredReadsPerStrand){
            return true;
        }

        // Check that the read set is paired (this check may not reject some sets that contain unpaired reads)
        if (readsWithSameUMI.getReads().size() % 2 == 1){
            return true;
        }

        // Check that the total number of reads (from both strands) exceeds a specified threshold
        return readsWithSameUMI.getReads().size() < minimumRequiredReadsPerUMI;
    }

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        final List<ReadFilter> readFilters = new ArrayList<>(6);
        readFilters.add(new WellformedReadFilter());
        readFilters.addAll(getDuplicateSetWalkerSpecificReadFilterList());
        return readFilters;
    }

    private static List<ReadFilter> getDuplicateSetWalkerSpecificReadFilterList() {
        List<ReadFilter> filters = new ArrayList<>(5);
        filters.add(ReadFilterLibrary.MAPPING_QUALITY_AVAILABLE);
        filters.add(ReadFilterLibrary.MAPPED);
        filters.add(ReadFilterLibrary.NOT_SECONDARY_ALIGNMENT);
        filters.add(ReadFilterLibrary.PASSES_VENDOR_QUALITY_CHECK);
        filters.add(ReadFilterLibrary.NON_ZERO_REFERENCE_LENGTH_ALIGNMENT);
        return filters;
    }
}
