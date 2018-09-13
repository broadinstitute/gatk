package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.activityprofile.ActivityProfileState;
import org.broadinstitute.hellbender.utils.clipping.ReadClipper;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadCoordinateComparator;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Region of the genome that gets assembled by the local assembly engine.
 */
public final class AssemblyRegion implements Locatable {

    private final SAMFileHeader header;

    /**
     * The reads included in this assembly region.  May be empty upon creation, and expand / contract
     * as reads are added or removed from this region.
     */
    private final List<GATKRead> reads;

    /**
     * An ordered list (by genomic coordinate) of the ActivityProfileStates that went
     * into this assembly region.  May be empty, which says that no supporting states were
     * provided when this region was created.
     */
    private final List<ActivityProfileState> supportingStates;

    /**
     * The raw span of this assembly region, not including the region extension
     */
    private final SimpleInterval activeRegionLoc;

    /**
     * The span of this assembly region on the genome, including the region extension
     */
    private final SimpleInterval extendedLoc;

    /**
     * The extension, in bp, of this region. The extension is >= 0 bp in size, and indicates how much padding was
     * requested for the region.
     */
    private final int extension;

    /**
     * Does this region represent an active region (all isActiveProbs above threshold) or
     * an inactive region (all isActiveProbs below threshold)?
     */
    private final boolean isActive;

    /**
     * The span of this assembly region, including the bp covered by all reads in this
     * region.  This union of extensionLoc and the loc of all reads in this region.
     *
     * Must be at least as large as extendedLoc, but may be larger when reads
     * partially overlap this region.
     */
    private SimpleInterval spanIncludingReads;

    /**
     * Indicates whether the region has been finalized
     */
    private boolean hasBeenFinalized;

    /**
     * Create a new AssemblyRegion containing no reads
     *
     * @param activeRegionLoc the span of this active region
     * @param supportingStates the states that went into creating this region, or null / empty if none are available.
     *                         If not empty, must have exactly one state for each bp in activeRegionLoc
     * @param isActive indicates whether this is an active region, or an inactive one
     * @param extension the active region extension to use for this active region
     */
    public AssemblyRegion( final SimpleInterval activeRegionLoc, final List<ActivityProfileState> supportingStates, final boolean isActive, final int extension , final SAMFileHeader header ) {
        Utils.nonNull(activeRegionLoc, "activeRegionLoc cannot be null");
        Utils.nonNull(header, "header cannot be null");
        Utils.validateArg( activeRegionLoc.size() > 0, () -> "Active region cannot be of zero size, but got " + activeRegionLoc);
        Utils.validateArg( extension >= 0, () -> "extension cannot be < 0 but got " + extension);

        this.header = header;
        this.reads = new ArrayList<>();
        this.activeRegionLoc = activeRegionLoc;
        this.supportingStates = supportingStates == null ? Collections.emptyList() : Collections.unmodifiableList(new ArrayList<>(supportingStates));
        this.isActive = isActive;
        this.extension = extension;
        this.extendedLoc = trimIntervalToContig(activeRegionLoc.getContig(), activeRegionLoc.getStart() - extension, activeRegionLoc.getEnd() + extension);
        this.spanIncludingReads = extendedLoc;

        checkStates(activeRegionLoc);
    }

    /**
     * Create a new SimpleInterval, bounding start and stop by the start and end of contig
     *
     * This function will return null if start and stop cannot be adjusted in any reasonable way
     * to be on the contig.  For example, if start and stop are both past the end of the contig,
     * there's no way to fix this, and null will be returned.
     *
     * @param contig our contig
     * @param start our start as an arbitrary integer (may be negative, etc)
     * @param stop our stop as an arbitrary integer (may be negative, etc)
     * @return a valid genome loc over contig, or null if a meaningful genome loc cannot be created
     */
    private SimpleInterval trimIntervalToContig(final String contig, final int start, final int stop) {
        final int contigLength = header.getSequence(contig).getSequenceLength();
        return IntervalUtils.trimIntervalToContig(contig, start, stop, contigLength);
    }

    private void checkStates(final SimpleInterval activeRegionLoc) {
        if ( ! this.supportingStates.isEmpty() ) {
            Utils.validateArg( this.supportingStates.size() == activeRegionLoc.size(), () ->
                    "Supporting states wasn't empty but it doesn't have exactly one state per bp in the active region: states " + this.supportingStates.size() + " vs. bp in region = " + activeRegionLoc.size());
            SimpleInterval lastStateLoc = null;
            for ( final ActivityProfileState state : this.supportingStates ) {
                if ( lastStateLoc != null ) {
                    if ( state.getLoc().getStart() != lastStateLoc.getStart() + 1 || !state.getLoc().getContig().equals(lastStateLoc.getContig())) {
                        throw new IllegalArgumentException("Supporting state has an invalid sequence: last state was " + lastStateLoc + " but next state was " + state);
                    }
                }
                lastStateLoc = state.getLoc();
            }
        }
    }

    /**
     * Simple interface to create an assembly region that isActive without any profile state
     */
    public AssemblyRegion(final SimpleInterval activeRegionLoc, final int extension, final SAMFileHeader header) {
        this(activeRegionLoc, Collections.<ActivityProfileState>emptyList(), true, extension, header);
    }

    @Override
    public String getContig() {
        return activeRegionLoc.getContig();
    }

    @Override
    public int getStart() {
        return activeRegionLoc.getStart();
    }

    @Override
    public int getEnd() {
        return activeRegionLoc.getEnd();
    }

    @Override
    public String toString() {
        return "AssemblyRegion "  + activeRegionLoc.toString() + " active?=" + isActive + " nReads=" + reads.size();
    }

    /**
     * Does this region represent an active region (all isActiveProbs above threshold) or
     * an inactive region (all isActiveProbs below threshold)?
     */
    public boolean isActive() {
        return isActive;
    }

    /**
     * Get the span of this assembly region including the extension value
     * @return a non-null SimpleInterval
     */
    public SimpleInterval getExtendedSpan() { return extendedLoc; }

    /**
     * Get the raw span of this assembly region (excluding the extension)
     * @return a non-null SimpleInterval
     */
    public SimpleInterval getSpan() { return activeRegionLoc; }

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
     * Returns the header for the reads in this region.
     */
    public SAMFileHeader getHeader(){
        return header;
    }

    /**
     * Intersect this assembly region with the allowed intervals, returning a list of active regions
     * that only contain locations present in intervals
     *
     * Note: modifications to the returned list have no effect on this region object.
     *
     * Note that the returned list may be empty, if this active region doesn't overlap the set at all
     *
     * Note that the resulting regions are all empty, regardless of whether the current active region has reads
     *
     * @param intervals a non-null set of intervals that are allowed
     * @return an ordered list of active region where each interval is contained within intervals
     */
    public List<AssemblyRegion> splitAndTrimToIntervals(final Set<SimpleInterval> intervals) {
        return intervals.stream()
                .filter(loc -> loc.overlaps(activeRegionLoc))
                .map(gl -> trim(gl, extension))
                .collect(Collectors.toList());
    }

    /**
     * Trim this region to just the span, producing a new assembly region without any reads that has only
     * the extent of newExtend intersected with the current extent
     * @param span the new extend of the active region we want
     * @param extensionSize the extensionSize size we want for the newly trimmed active region
     * @return a non-null, empty assembly region
     */
    public AssemblyRegion trim(final SimpleInterval span, final int extensionSize) {
        Utils.nonNull(span, "Active region extent cannot be null");
        Utils.validateArg( extensionSize >= 0, "the extensionSize size must be 0 or greater");
        final int extendStart = Math.max(1,span.getStart() - extensionSize);
        final int maxStop = header.getSequence(span.getContig()).getSequenceLength();
        final int extendStop = Math.min(span.getEnd() + extensionSize, maxStop);
        final SimpleInterval extendedSpan = new SimpleInterval(span.getContig(), extendStart, extendStop);
        return trim(span, extendedSpan);

//TODO - Inconsistent support of substates trimming. Check lack of consistency!!!!
//        final GenomeLoc subLoc = getLocation().intersect(span);
//        final int subStart = subLoc.getStart() - getLocation().getStart();
//        final int subEnd = subStart + subLoc.size();
//        final List<ActivityProfileState> subStates = supportingStates.isEmpty() ? supportingStates : supportingStates.subList(subStart, subEnd);
//        return new ActiveRegion( subLoc, subStates, isActive, genomeLocParser, extensionSize );

    }

    /**
     * Equivalent to trim(span,span).
     */
    public AssemblyRegion trim(final SimpleInterval span) {
        return trim(span, span);
    }

    /**
     * Trim this region to no more than the span, producing a new assembly region with properly trimmed reads that
     * attempts to provide the best possible representation of this region covering the span.
     *
     * The challenge here is that span may (1) be larger than can be represented by this assembly region
     * + its original extension and (2) the extension must be symmetric on both sides.  This algorithm
     * therefore determines how best to represent span as a subset of the span of this
     * region with a padding value that captures as much of the span as possible.
     *
     * For example, suppose this active region is
     *
     * Active:    100-200 with extension of 50, so that the true span is 50-250
     * NewExtent: 150-225 saying that we'd ideally like to just have bases 150-225
     *
     * Here we represent the assembly region as a region from 150-200 with 25 bp of padding.
     *
     * The overall constraint is that the region can never exceed the original region, and
     * the extension is chosen to maximize overlap with the desired region
     *
     * @param span the new extend of the active region we want
     * @return a non-null, empty active region
     */
    public AssemblyRegion trim(final SimpleInterval span, final SimpleInterval extendedSpan) {
        Utils.nonNull(span, "Active region extent cannot be null");
        Utils.nonNull(extendedSpan, "Active region extended span cannot be null");
        Utils.validateArg(extendedSpan.contains(span), "The requested extended span must fully contain the requested span");

        final SimpleInterval subActive = getSpan().intersect(span);
        final int requiredOnRight = Math.max(extendedSpan.getEnd() - subActive.getEnd(), 0);
        final int requiredOnLeft = Math.max(subActive.getStart() - extendedSpan.getStart(), 0);
        final int requiredExtension = Math.min(Math.max(requiredOnLeft, requiredOnRight), getExtension());

        final AssemblyRegion result = new AssemblyRegion( subActive, Collections.<ActivityProfileState>emptyList(), isActive, requiredExtension, header );

        final List<GATKRead> myReads = getReads();
        final SimpleInterval resultExtendedLoc = result.getExtendedSpan();
        final int resultExtendedLocStart = resultExtendedLoc.getStart();
        final int resultExtendedLocStop = resultExtendedLoc.getEnd();

        final List<GATKRead> trimmedReads = new ArrayList<>(myReads.size());
        for( final GATKRead read : myReads ) {
            final GATKRead clippedRead = ReadClipper.hardClipToRegion(read, resultExtendedLocStart, resultExtendedLocStop);
            if( result.readOverlapsRegion(clippedRead) && !clippedRead.isEmpty() ) {
                trimmedReads.add(clippedRead);
            }
        }
        result.clearReads();

        trimmedReads.sort(new ReadCoordinateComparator(header));
        result.addAll(trimmedReads);
        return result;
    }

    /**
     * Returns true if read would overlap the extended extent of this region
     * @param read the read we want to test
     * @return true if read can be added to this region, false otherwise
     */
    public boolean readOverlapsRegion(final GATKRead read) {
        if ( read.isEmpty() || read.getStart() > read.getEnd() ) {
            return false;
        }

        final SimpleInterval readLoc = new SimpleInterval( read );
        return readLoc.overlaps(extendedLoc);
    }

    /**
     * Add read to this region
     *
     * Read must have alignment start >= than the last read currently in this active region.
     *
     * @throws IllegalArgumentException if read doesn't overlap the extended region of this active region
     *
     * @param read a non-null GATKRead
     */
    public void add( final GATKRead read ) {
        Utils.nonNull(read, "Read cannot be null");
        final SimpleInterval readLoc = new SimpleInterval( read );
        Utils.validateArg(readOverlapsRegion(read), () ->
                "Read location " + readLoc + " doesn't overlap with active region extended span " + extendedLoc);

        spanIncludingReads = spanIncludingReads.mergeWithContiguous( readLoc );

        if ( ! reads.isEmpty() ) {
            final GATKRead lastRead = reads.get(size() - 1);
            Utils.validateArg(Objects.equals(lastRead.getContig(), read.getContig()), () ->
                    "Attempting to add a read to ActiveRegion not on the same contig as other reads: lastRead " + lastRead + " attempting to add " + read);
            Utils.validateArg( read.getStart() >= lastRead.getStart(), () ->
                    "Attempting to add a read to ActiveRegion out of order w.r.t. other reads: lastRead " + lastRead + " at " + lastRead.getStart() + " attempting to add " + read + " at " + read.getStart());
        }

        reads.add( read );
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
        spanIncludingReads = extendedLoc;
        reads.clear();
    }

    /**
     * Remove all of the reads in readsToRemove from this region
     * @param readsToRemove the set of reads we want to remove
     */
    public void removeAll( final Collection<GATKRead> readsToRemove ) {
        Utils.nonNull(readsToRemove);
        reads.removeAll(readsToRemove);
        spanIncludingReads = extendedLoc;
        for (final GATKRead read : reads) {
            spanIncludingReads = spanIncludingReads.mergeWithContiguous(read);
        }
    }

    /**
     * Add all readsToAdd to this region
     * @param readsToAdd a collection of readsToAdd to add to this active region
     */
    public void addAll(final Collection<GATKRead> readsToAdd){
        Utils.nonNull(readsToAdd).forEach(r -> add(r));
    }

    /**
     * Get the extension applied to this region
     *
     * The extension is >= 0 bp in size, and indicates how much padding was requested for the region
     *
     * @return the size in bp of the region extension
     */
    public int getExtension() { return extension; }

    /**
     * The span of this assembly region, including the bp covered by all reads in this
     * region.  This union of extensionLoc and the loc of all reads in this region.
     *
     * Must be at least as large as extendedLoc, but may be larger when reads
     * partially overlap this region.
     */
    public SimpleInterval getReadSpanLoc() {
        return spanIncludingReads;
    }

    /**
     * An ordered list (by genomic coordinate) of the ActivityProfileStates that went
     * into this active region.  May be empty, which says that no supporting states were
     * provided when this region was created.
     * The returned list is unmodifiable.
     */
    public List<ActivityProfileState> getSupportingStates() {
        return supportingStates;
    }

    /**
     * See #getActiveRegionReference but using the span including regions not the extended span
     */
    public byte[] getFullReference( final IndexedFastaSequenceFile referenceReader ) {
        return getFullReference(referenceReader, 0);
    }

    /**
     * See #getActiveRegionReference but using the span including regions not the extended span
     */
    public byte[] getFullReference( final IndexedFastaSequenceFile referenceReader, final int padding ) {
        return getReference(referenceReader, padding, spanIncludingReads);
    }

    /**
     * Get the reference bases from referenceReader spanned by the extended location of this region,
     * including additional padding bp on either side.  If this expanded region would exceed the boundaries
     * of the active region's contig, the returned result will be truncated to only include on-genome reference
     * bases.
     *
     * @param referenceReader the source of the reference genome bases
     * @param padding the padding, in BP, we want to add to either side of this active region extended region
     * @param genomeLoc a non-null genome loc indicating the base span of the bp we'd like to get the reference for
     * @return a non-null array of bytes holding the reference bases in referenceReader
     */
    private static byte[] getReference(final ReferenceSequenceFile referenceReader, final int padding, final SimpleInterval genomeLoc) {
        Utils.nonNull(referenceReader, "referenceReader cannot be null");
        Utils.nonNull(genomeLoc, "genomeLoc cannot be null");
        Utils.validateArg( padding >= 0, () -> "padding must be a positive integer but got " + padding);
        Utils.validateArg( genomeLoc.size() > 0, () -> "GenomeLoc must have size > 0 but got " + genomeLoc);

        return referenceReader.getSubsequenceAt( genomeLoc.getContig(),
                Math.max(1, genomeLoc.getStart() - padding),
                Math.min(referenceReader.getSequenceDictionary().getSequence(genomeLoc.getContig()).getSequenceLength(), genomeLoc.getEnd() + padding) ).getBases();
    }

    /**
     * See {@link #getAssemblyRegionReference} with padding == 0
     */
    public byte[] getAssemblyRegionReference( final ReferenceSequenceFile referenceReader ) {
        return getAssemblyRegionReference(referenceReader, 0);
    }

    /**
     * Get the reference bases from referenceReader spanned by the extended location of this active region,
     * including additional padding bp on either side.  If this expanded region would exceed the boundaries
     * of the active region's contig, the returned result will be truncated to only include on-genome reference
     * bases
     *
     * @param referenceReader the source of the reference genome bases
     * @param padding the padding, in BP, we want to add to either side of this active region extended region
     * @return a non-null array of bytes holding the reference bases in referenceReader
     */
    public byte[] getAssemblyRegionReference(final ReferenceSequenceFile referenceReader, final int padding ) {
        return getReference(referenceReader, padding, extendedLoc);
    }

    /**
     * Is this region equal to other, excluding any reads in either region in the comparison
     * @param other the other active region we want to test
     * @return true if this region is equal, excluding any reads and derived values, to other
     */
    public boolean equalsIgnoreReads(final AssemblyRegion other) {
        if ( other == null ) {
            return false;
        }
        if ( ! activeRegionLoc.equals(other.activeRegionLoc)) {
            return false;
        }
        if ( isActive() != other.isActive()) {
            return false;
        }
        if ( extension != other.extension ) {
            return false;
        }
        return extendedLoc.equals(other.extendedLoc);
    }

    public void setFinalized(final boolean value) {
        hasBeenFinalized = value;
    }

    public boolean isFinalized() {
        return hasBeenFinalized;
    }

}
