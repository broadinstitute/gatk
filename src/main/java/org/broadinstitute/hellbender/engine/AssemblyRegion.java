package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.mutect.AlignmentData;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.clipping.ReadClipper;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadCoordinateComparator;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Region of the genome that gets assembled by the local assembly engine.
 *
 * As AssemblyRegion is defined by two intervals -- a primary interval containing a territory for variant calling and a second,
 * padded, interval for assembly -- as well as the reads overlapping the padded interval.  Although we do not call variants in the padded interval,
 * assembling over a larger territory improves calls in the primary territory.
 *
 * This concept is complicated somewhat by the fact that these intervals are mutable and the fact that the AssemblyRegion onject lives on after
 * assembly during local realignment during PairHMM.  Here is an example of the life cycle of an AssemblyRegion:
 *
 * Suppose that the HaplotypeCaller engine finds an evidence for a het in a pileup at locus 400 -- that is, it produces
 * an {@code ActivityProfileState} with non-zero probability at site 400 and passes it to its {@code ActivityProfile}.
 * The {@code ActivityProfile} eventually produces an AssemblyRegion based on the {@code AssemblyRegionArgumentCollection} parameters.
 * Let's suppose that this initial region has primary span 350-450 and padded span 100 - 700.
 *
 * Next, the assembly engine assembles all reads that overlap the padded interval to find variant haplotypes and the variants
 * they contain.  The AssemblyRegion is then trimmed down to a new primary interval bound by all assembled variants within the original primary interval
 * and a new padded interval.  The amount of padding of the new padded interval around the variants depends on the needs of local realignment
 * and as such need not equal the original padding that was used for assembly.
 */
public final class AssemblyRegion implements Locatable {

    private final SAMFileHeader header;

    /**
     * The reads included in this assembly region.  May be empty upon creation, and expand / contract
     * as reads are added or removed from this region.
     */
    private final List<GATKRead> reads;

    /**
     * The reads are specifically used for haplotype generation to kmerize reads to match with haplotype kmers.
     */
    private final List<GATKRead> hardClippedPileupReads;

    /**
     * The active span in which this AssemblyRegion is responsible for calling variants
     */
    private final SimpleInterval activeSpan;

    /**
     * The padded span in which we perform assembly etc in order to call variants within the active span
     */
    private final SimpleInterval paddedSpan;

    /**
     * Does this region represent an active region (all isActiveProbs above threshold) or
     * an inactive region (all isActiveProbs below threshold)?
     */
    private boolean isActive;

    /**
     * Indicates whether the region has been finalized
     */
    private boolean hasBeenFinalized;

    private List<AlignmentData> alignmentData = new ArrayList<>();

    /**
     * Create a new AssemblyRegion containing no reads
     *  @param activeSpan the span of this active region
     * @param isActive indicates whether this is an active region, or an inactive one
     * @param padding the active region padding to use for this active region
     */
    public AssemblyRegion(final SimpleInterval activeSpan, final boolean isActive, final int padding, final SAMFileHeader header) {
        this(activeSpan, makePaddedSpan(activeSpan, padding, header), isActive, header);
    }

    private static SimpleInterval makePaddedSpan(final SimpleInterval activeSpan, final int padding, final SAMFileHeader header) {
        final String contig = activeSpan.getContig();
        final SAMSequenceRecord sequence = header.getSequence(contig);
        if( sequence == null) {
            throw new UserException.MissingContigInSequenceDictionary(contig, header.getSequenceDictionary());
        }
        return IntervalUtils.trimIntervalToContig(contig, activeSpan.getStart() - padding, activeSpan.getEnd() + padding, sequence.getSequenceLength());
    }

    /**
     * Create a new AssemblyRegion containing no reads
     * @param activeSpan the span of this active region
     * @param paddedSpan    the padded span of this active region
     * @param isActive indicates whether this is an active region, or an inactive one
     */
    public AssemblyRegion(final SimpleInterval activeSpan, final SimpleInterval paddedSpan, final boolean isActive, final SAMFileHeader header) {
        this.header = Utils.nonNull(header);
        this.activeSpan = Utils.nonNull(activeSpan);
        this.paddedSpan = Utils.nonNull(paddedSpan);

        Utils.validateArg( activeSpan.size() > 0, () -> "Active region cannot be of zero size, but got " + activeSpan);
        Utils.validate(paddedSpan.contains(activeSpan), "Padded span must contain active span.");

        reads = new ArrayList<>();
        hardClippedPileupReads = new ArrayList<>();
        this.isActive = isActive;
    }

    /**
     * Simple interface to create an assembly region that isActive without any profile state
     */
    public AssemblyRegion(final SimpleInterval activeSpan, final int padding, final SAMFileHeader header) {
        this(activeSpan, true, padding, header);
    }

    /**
     * Method for obtaining the alignment data which is attached to the assembly region.
     *
     * @return The list of AlignmentData objects associated with ActiveRegion.
     */
    public List<AlignmentData> getAlignmentData() {
        return alignmentData;
    }

    /**
     * Method for adding alignment data to the collection of AlignmentData associated with
     * the ActiveRegion.
     */
    public void addAllAlignmentData(List<AlignmentData> alignmentData) {
        this.alignmentData.addAll(alignmentData);
    }

    @Override
    public String getContig() {
        return activeSpan.getContig();
    }

    @Override
    public int getStart() {
        return activeSpan.getStart();
    }

    @Override
    public int getEnd() {
        return activeSpan.getEnd();
    }

    @Override
    public String toString() {
        return "AssemblyRegion "  + activeSpan.toString() + " active?=" + isActive + " nReads=" + reads.size();
    }

    /**
     * Does this region represent an active region (all isActiveProbs above threshold) or
     * an inactive region (all isActiveProbs below threshold)?
     */
    public boolean isActive() {
        return isActive;
    }

    /**
     * Override activity state of the region
     *
     * Note: Changing the isActive state after construction is a debug-level operation that only engine classes
     * like AssemblyRegionWalker should be able to do
     *
     * @param value new activity state of this region
     */
    void setIsActive(final boolean value) {
        isActive = value;
    }

    /**
     * Get the span of this assembly region including the padding value
     * @return a non-null SimpleInterval
     */
    public SimpleInterval getPaddedSpan() { return paddedSpan; }

    /**
     * Get the raw span of this assembly region (excluding the padding)
     * @return a non-null SimpleInterval
     */
    public SimpleInterval getSpan() { return activeSpan; }

    /**
     * Get an unmodifiable copy of the list of reads currently in this assembly region.
     *
     * The reads are sorted by their coordinate position.
     * @return an unmodifiable and inmutable copy of the reads in the assembly region.
    */
    public List<GATKRead> getReads(){
        return Collections.unmodifiableList(new ArrayList<>(reads));
    }

    /**
     * Get an unmodifiable copy of the list of reads currently in this assembly region.
     *
     * The reads are sorted by their coordinate position.
     * @return an unmodifiable and inmutable copy of the reads in the assembly region.
     */
    public List<GATKRead> getHardClippedPileupReads(){
        return Collections.unmodifiableList(new ArrayList<>(hardClippedPileupReads));
    }

    /**
     * Returns the header for the reads in this region.
     */
    public SAMFileHeader getHeader(){
        return header;
    }

    /**
     * Trim this region to just the span, producing a new assembly region without any reads that has only
     * the extent of newExtend intersected with the current extent
     * @param span the new extend of the active region we want
     * @param padding the padding size we want for the newly trimmed active region
     * @return a non-null, empty assembly region
     */
    public AssemblyRegion trim(final SimpleInterval span, final int padding) {
        Utils.nonNull(span, "Active region extent cannot be null");
        Utils.validateArg( padding >= 0, "the padding size must be 0 or greater");
        final SimpleInterval paddedSpan = span.expandWithinContig(padding, header.getSequenceDictionary());
        return trim(span, paddedSpan);
    }

    /**
     * Trim this region to no more than the span, producing a new assembly region with properly trimmed reads that
     * attempts to provide the best possible representation of this region covering the span.
     *
     * The challenge here is that span may (1) be larger than can be represented by this assembly region
     * + its original padding and (2) the padding must be symmetric on both sides.  This algorithm
     * therefore determines how best to represent span as a subset of the span of this
     * region with a padding value that captures as much of the span as possible.
     *
     * For example, suppose this active region is
     *
     * Active:    100-200 with padding of 50, so that the true span is 50-250
     * NewExtent: 150-225 saying that we'd ideally like to just have bases 150-225
     *
     * Here we represent the assembly region as a region from 150-200 with 25 bp of padding.
     *
     * The overall constraint is that the region can never exceed the original region, and
     * the padding is chosen to maximize overlap with the desired region
     *
     * @param span the new extend of the active region we want
     * @return a non-null, empty active region
     */
    public AssemblyRegion trim(final SimpleInterval span, final SimpleInterval paddedSpan) {
        Utils.nonNull(span, "Active region extent cannot be null");
        Utils.nonNull(paddedSpan, "Active region padded span cannot be null");
        Utils.validateArg(paddedSpan.contains(span), "The requested padded span must fully contain the requested span");

        final SimpleInterval newActiveSpan = getSpan().intersect(span);
        final SimpleInterval newPaddedSpan = getPaddedSpan().intersect(paddedSpan);

        final AssemblyRegion result = new AssemblyRegion( newActiveSpan, newPaddedSpan, isActive, header );

        final List<GATKRead> trimmedReads = reads.stream()
                .map(read -> {GATKRead clipped = ReadClipper.hardClipToRegion(read, newPaddedSpan.getStart(), newPaddedSpan.getEnd());
                              return clipped;})
                .filter(read -> !read.isEmpty() && read.overlaps(result.paddedSpan))
                .sorted(new ReadCoordinateComparator(header))
                .collect(Collectors.toList());

        result.clearReads();
        result.addAll(trimmedReads);
        return result;
    }


    /**
     * Add read to this region
     *
     * Read must have alignment start >= than the last read currently in this active region.
     *
     * @throws IllegalArgumentException if read doesn't overlap the padded region of this active region
     *
     * @param read a non-null GATKRead
     */
    public void add( final GATKRead read ) {
        addToReadCollectionAndValidate(read, reads);
    }

    private void addToReadCollectionAndValidate(final GATKRead read, final List<GATKRead> collection) {
        Utils.nonNull(read, "Read cannot be null");
        final SimpleInterval readLoc = new SimpleInterval( read );
        Utils.validateArg(paddedSpan.overlaps(read), () ->
                "Read location " + readLoc + " doesn't overlap with active region padded span " + paddedSpan);

        if ( ! collection.isEmpty() ) {
            final GATKRead lastRead = collection.get(collection.size() - 1);
            Utils.validateArg(Objects.equals(lastRead.getContig(), read.getContig()), () ->
                    "Attempting to add a read to ActiveRegion not on the same contig as other reads: lastRead " + lastRead + " attempting to add " + read);
            Utils.validateArg( read.getStart() >= lastRead.getStart(), () ->
                    "Attempting to add a read to ActiveRegion out of order w.r.t. other reads: lastRead " + lastRead + " at " + lastRead.getStart() + " attempting to add " + read + " at " + read.getStart());
        }

        collection.add( read );
    }

    /**
     * Get the number of reads currently in this region
     * @return an integer >= 0
     */
    public int size() { return reads.size(); }

    /**
     * Clear all of the reads currently in this region
     */
    public void clearReads() {
        reads.clear();
        hardClippedPileupReads.clear();
    }

    /**
     * Remove all of the reads in readsToRemove from this region
     * @param readsToRemove the set of reads we want to remove
     */
    public void removeAll( final Collection<GATKRead> readsToRemove ) {
        Utils.nonNull(readsToRemove);
        reads.removeAll(readsToRemove);
    }

    /**
     * Add all readsToAdd to this region
     * @param readsToAdd a collection of readsToAdd to add to this active region
     */
    public void addAll(final Collection<GATKRead> readsToAdd){
        Utils.nonNull(readsToAdd).forEach(r -> add(r));
    }

    public void addHardClippedPileupReads(final Collection<GATKRead> readsToAdd) {
        Utils.nonNull(readsToAdd).forEach(r -> addToReadCollectionAndValidate(r, hardClippedPileupReads));
    }

    /**
     * Get the reference bases from referenceReader spanned by the padded span of this region,
     * including additional padding bp on either side.  If this expanded region would exceed the boundaries
     * of the active region's contig, the returned result will be truncated to only include on-genome reference
     * bases.
     *
     * @param referenceReader the source of the reference genome bases
     * @param padding the padding, in BP, we want to add to either side of this active region padded span
     * @param genomeLoc a non-null genome loc indicating the base span of the bp we'd like to get the reference for
     * @return a non-null array of bytes holding the reference bases in referenceReader
     */
    private static byte[] getReference(final ReferenceSequenceFile referenceReader, final int padding, final SimpleInterval genomeLoc) {
        Utils.nonNull(referenceReader, "referenceReader cannot be null");
        Utils.nonNull(genomeLoc, "genomeLoc cannot be null");
        Utils.validateArg( padding >= 0, () -> "padding must be a positive integer but got " + padding);
        Utils.validateArg( genomeLoc.size() > 0, () -> "GenomeLoc must have size > 0 but got " + genomeLoc);

        final String contig = genomeLoc.getContig();
        final SAMSequenceDictionary sequenceDictionary = referenceReader.getSequenceDictionary();
        final SAMSequenceRecord sequence = sequenceDictionary.getSequence(contig);
        if ( sequence == null ) {
            throw new UserException.MissingContigInSequenceDictionary("Contig: " + contig + " not found in reference dictionary." +
                    "\nPlease check that you are using a compatible reference for your data." +
                    "\nReference Contigs: " + ReadUtils.prettyPrintSequenceRecords(sequenceDictionary));
        }
        return referenceReader.getSubsequenceAt(contig,
                Math.max(1, genomeLoc.getStart() - padding),
                Math.min(sequence.getSequenceLength(), genomeLoc.getEnd() + padding)).getBases();
    }

    /**
     * See {@link #getAssemblyRegionReference} with padding == 0
     */
    public byte[] getAssemblyRegionReference( final ReferenceSequenceFile referenceReader ) {
        return getAssemblyRegionReference(referenceReader, 0);
    }

    /**
     * Get the reference bases from referenceReader spanned by the padded span of this active region,
     * including additional padding bp on either side.  If this expanded region would exceed the boundaries
     * of the active region's contig, the returned result will be truncated to only include on-genome reference
     * bases
     *
     * @param referenceReader the source of the reference genome bases
     * @param padding the padding, in BP, we want to add to either side of this active region padded region
     * @return a non-null array of bytes holding the reference bases in referenceReader
     */
    public byte[] getAssemblyRegionReference(final ReferenceSequenceFile referenceReader, final int padding ) {
        return getReference(referenceReader, padding, paddedSpan);
    }

    public void setFinalized(final boolean value) {
        hasBeenFinalized = value;
    }

    public boolean isFinalized() {
        return hasBeenFinalized;
    }

}
