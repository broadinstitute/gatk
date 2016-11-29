package org.broadinstitute.hellbender.tools.exome;

import com.google.api.services.genomics.model.Read;
import com.google.cloud.dataflow.sdk.repackaged.com.google.common.collect.Iterators;
import htsjdk.samtools.*;
import htsjdk.samtools.util.Locatable;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.*;

/**
 * Read that is the result of merging two read records.
 * <p>
 *     This class is tailored to be used by {@link CalculateTargetBaseCallCoverage}
 *     when fragment based processing is requested.
 * </p>
 * <p>
 *     Its current implementation has limited functionality and is not intended as
 *     a general read-pair merging utility class.
 * </p>
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
final class MergedGATKReadPair implements GATKRead {

    private final String name;

    private final byte[] bases;

    private final byte[] qualities;

    private final Cigar cigar;

    private final String readGroup;

    private final SimpleInterval interval;

    private final int mappingQuality;

    /**
     * Construct an instance given the reads key properties.
     * @param name the read name.
     * @param readGroup the read-group name.
     * @param interval the read alignment enclosing interval.
     * @param bases the base call sequence.
     * @param qualities the base call qualities sequence.
     * @param cigar the new read alignment cigar.
     */
    private MergedGATKReadPair(final String name, final String readGroup,
                               final SimpleInterval interval, final byte[] bases,
                               final byte[] qualities, final Cigar cigar, final int mappingQuality) {
        this.name = name;
        this.bases = bases;
        this.qualities = qualities;
        this.cigar = cigar;
        this.readGroup = readGroup;
        this.interval = interval;
        this.mappingQuality = mappingQuality;
    }

    /**
     * Compose a read record out of a read-pair.
     *
     * @param first the first read in the pair.
     * @param second the second read in the pair.
     * @return never {@code null}.
     * @throws IllegalArgumentException if any of the input read pair is null or it is unmapped.
     */
    public static GATKRead mergeReadPair(final GATKRead first, final GATKRead second) {
        Utils.nonNull(first, "the first read cannot be null");
        Utils.nonNull(second, "the second read cannot be null");
        Utils.validateArg(!first.isUnmapped() && first.getContig().equals(second.getContig()), "the first and second reads must align to the same contig");
        return first.getStart() <= second.getStart() ? mergeReadSortedPair(first, second, true) : mergeReadSortedPair(second, first, false);
    }

    /**
     * Compose a read record out of a read-pair.
     *
     * @param left the first read in the pair, the one whose alignment starts earlier on the reference.
     * @param right the second read in the pair, the one that aligns later on the reference.
     * @param leftOverRight if {@code true} for overlapped based calls with the same quality but different base
     *                       the left read base is used. If {@code false} the right based is used instead.
     * @return never {@code null}.
     */
    private static MergedGATKReadPair mergeReadSortedPair(final GATKRead left, final GATKRead right, final boolean leftOverRight) {
        final List<AlignmentBlock> leftAlignmentBlocks = extractAlignmentBlocks(left, leftOverRight);
        final List<AlignmentBlock> rightAlignmentBlocks = extractAlignmentBlocks(right, !leftOverRight);
        final List<AlignmentBlock> combinedAlignmentBlocks = combineAlignmentBlocks(leftAlignmentBlocks, rightAlignmentBlocks);
        final SimpleInterval interval = calculateMatchInterval(left.getContig(), combinedAlignmentBlocks);
        final Pair<byte[], byte[]> basesAndQuals = extractBasesAndQualsFromAlignmentBlockList(combinedAlignmentBlocks);
        final Cigar cigar = calculateCigarFromAlignmentBlocks(combinedAlignmentBlocks);
        return new MergedGATKReadPair(left.getName(), left.getReadGroup(), interval, basesAndQuals.getLeft(), basesAndQuals.getRight(),
                cigar, Math.max(left.getMappingQuality(), right.getMappingQuality()));
    }

    /**
     * Compose the cigar for the sequence of alignment blocks provided adding any deletion operations
     * needed to fill the gap between the blocks.
     * <p>At this point we can assume that the gaps are not 0 so there is always going to be
     *   at least one base deletion between consecutive alignment blocks in the input</p>
     */
    private static Cigar calculateCigarFromAlignmentBlocks(final List<AlignmentBlock> combinedAlignmentBlocks) {
        final List<CigarElement> elements = new ArrayList<>(combinedAlignmentBlocks.size() * 2 - 1);
        final Iterator<AlignmentBlock> it = combinedAlignmentBlocks.listIterator();
        AlignmentBlock previousBlock = null;
        while (it.hasNext()) {
            final AlignmentBlock nextBlock = it.next();
            if (previousBlock != null) {
                final int gap = nextBlock.referenceStart - previousBlock.referenceEnd;
                elements.add(new CigarElement(gap, CigarOperator.D));
            }
            elements.add(nextBlock.toCigarElement());
            previousBlock = nextBlock;
        }
        return new Cigar(elements);
    }

    /**
     * Combine two sorted alignment block lists into one by resolving all overlaps appropriately.
     * @param leftAlignmentBlocks the left read alignment blocks.
     * @param rightAlignmentBlocks the right read alignment blocks.
     * @return the new combined read alignment blocks.
     */
    private static List<AlignmentBlock> combineAlignmentBlocks(final List<AlignmentBlock> leftAlignmentBlocks,
                                                               final List<AlignmentBlock> rightAlignmentBlocks) {
        final LinkedList<AlignmentBlock> result = new LinkedList<>();
        final Comparator<AlignmentBlock> comparator =
                Comparator.comparing(block -> block.referenceStart);

        final Iterator<AlignmentBlock> iterator = Iterators.mergeSorted(Arrays.asList(leftAlignmentBlocks.iterator(), rightAlignmentBlocks.iterator()),
                comparator);
        result.add(iterator.next());
        while (iterator.hasNext()) {
            final AlignmentBlock last = result.getLast();
            final AlignmentBlock next = iterator.next();
            if (last.referenceEnd >= next.referenceStart) {
                result.removeLast(); // remove the last block as we are going to merge it with the next.
                result.add(mergeOverlappingAlignmentBlocks(last, next));
            } else {
                result.add(next);
            }
        }
        return result;
    }

    /**
     * Merges the bases and qualities from two overlapping alignment blocks.
     * <p>Assumes that there is some overlap between the two blocks and that the {@code left} block starts
     *  earlier aligned against the reference (i.e. {@code left.referenceStart <= right.referenceStart}).</p>
     * @param left the left alignment block.
     * @param right the right alignment block.
     * @return never null.
     */
    private static AlignmentBlock mergeOverlappingAlignmentBlocks(final AlignmentBlock left, final AlignmentBlock right) {

        // The resulting alignment block is divided into a preceding 'left-overhang' a central 'overlap'
        // and a following 'right-overhang'.
        // The 'left-overhang base and qualities come from the left alignment block,
        // The central overlap base and qualities are the consensus between the overlapping calls from
        // the left and right blocks. This block may have size == 0 if the blocks are contiguous but not overlapping.
        // Finally the right-overhang base and qualities may come from either the left or right block
        // depending on which one has the largest reference-end position.

        // We start by calculating some useful lengths and offsets.
        final int referenceEnd = Math.max(left.referenceEnd, right.referenceEnd);
        final int totalLength = referenceEnd - left.referenceStart;
        final int leftOverhangLength = right.referenceStart - left.referenceStart;
        final int overlapLength = Math.min(left.referenceEnd, right.referenceEnd) - right.referenceStart;
        final int rightOverhangLength = totalLength - leftOverhangLength - overlapLength;
        final int rightOverhangOffset = leftOverhangLength + overlapLength;

        // the resulting base and quals arrays.
        final byte[] bases = new byte[totalLength];
        final byte[] quals = new byte[totalLength];

        mergeOverlappingAlignmentBlockLeftOverhangSection(left, leftOverhangLength, bases, quals);

        mergeOverlappingAlignmentBlockOverlappingSection(left, right, leftOverhangLength, overlapLength, bases, quals);

        mergeOverlappingAlignmentBlocksRightOverhangSection(left, right, overlapLength,
                rightOverhangLength, rightOverhangOffset, bases, quals);

        return new AlignmentBlock(bases, quals, 0, left.referenceStart, totalLength, false);
    }

    private static void mergeOverlappingAlignmentBlockLeftOverhangSection(final AlignmentBlock left,
                                                                          final int leftOverhangLength,
                                                                          final byte[] resultBases,
                                                                          final byte[] resultQuals) {
        System.arraycopy(left.bases, left.readStart, resultBases, 0, leftOverhangLength);
        System.arraycopy(left.quals, left.readStart, resultQuals, 0, leftOverhangLength);
    }

    private static void mergeOverlappingAlignmentBlocksRightOverhangSection(final AlignmentBlock left,
                                                                            final AlignmentBlock right,
                                                                            final int overlapLength,
                                                                            final int rightOverhangLength,
                                                                            final int rightOverhangOffset,
                                                                            final byte[] bases,
                                                                            final byte[] quals) {
        // Determine block where the base/quals come from
        // depends on which one extends the farthest along the reference.
        final boolean rightBlockExtendsFurther = left.referenceEnd < right.referenceEnd;

        // Get the appropriate source of bases and quals (block and offset within the block).
        final AlignmentBlock rightOverhangBlock = rightBlockExtendsFurther ? right : left;
        final int rightOverhangBlockOffset = rightBlockExtendsFurther ? right.readStart + overlapLength : left.readStart + rightOverhangOffset;

        System.arraycopy(rightOverhangBlock.bases, rightOverhangBlockOffset, bases, rightOverhangOffset, rightOverhangLength);
        System.arraycopy(rightOverhangBlock.quals, rightOverhangBlockOffset, quals, rightOverhangOffset, rightOverhangLength);
    }

    /**
     * Merges the base calls and qualities in overlapping alignment blocks' overlap section.
     *
     * <p>
     * The quality of an overlapping call is the maximum of the two qualities.
     * The base is the call of the one associated with the maximum quality (left block's call in case of a tie).
     * </p>
     * <p>
     * Notice that when there are mismatched based calls between calls the quality is
     * kept as the maximum quality which is kind of unintuitive. This is done in order to
     * preserve the behavior in GATK3's DepthOfCoverage tool when COUNT_FRAGMENTS is
     * selected which is the default in XHMM
     * (minBaseQuality is set to 0 so it does not really matter anyway if we would do something different).
     * </p>
     *
     * @param left the left overlapping alignment block.
     * @param right the right overlapping alignment block.
     * @param leftOverhangLength the size of the left-alignment block's left-overhang.
     * @param overlapLength the overlap size.
     * @param bases the output bases array.
     * @param quals the output quals array.
     */
    private static void mergeOverlappingAlignmentBlockOverlappingSection(final AlignmentBlock left,
                                                                         final AlignmentBlock right,
                                                                         final int leftOverhangLength,
                                                                         final int overlapLength,
                                                                         final byte[] bases,
                                                                         final byte[] quals) {
        final int leftBlockOverlapOffset = left.readStart + leftOverhangLength;
        for (int i = 0; i < overlapLength; i++) {
            final int leftIndex = leftBlockOverlapOffset + i;
            final int rightIndex = right.readStart + i;
            final int thisIndex = leftOverhangLength + i;
            final byte leftBase = left.bases[leftIndex];
            final byte leftQual = left.quals[leftIndex];
            final byte rightBase = right.bases[rightIndex];
            final byte rightQual = right.quals[rightIndex];
            bases[thisIndex] = (leftQual > rightQual) ? leftBase
                             : (rightQual > leftQual) ? rightBase
                             : (right.priority) ? rightBase : leftBase;
            quals[thisIndex] = (byte) Math.max(leftQual, rightQual);
        }
    }

    private static Pair<byte[], byte[]> extractBasesAndQualsFromAlignmentBlockList(final List<AlignmentBlock> alignmentBlocks) {
        final int totalLength = alignmentBlocks.stream().mapToInt(c -> c.length).sum();
        final byte[] bases = new byte[totalLength];
        final byte[] quals = new byte[totalLength];
        int nextIndex = 0;
        for (final AlignmentBlock block : alignmentBlocks) {
            System.arraycopy(block.bases, block.readStart, bases, nextIndex, block.length);
            System.arraycopy(block.quals, block.readStart, quals, nextIndex, block.length);
            nextIndex += block.length;
        }
        return new ImmutablePair<>(bases, quals);
    }

    @Override
    public String getName() {
        return name;
    }

    @Override
    public void setName(final String name) {
        throw new UnsupportedOperationException();
    }

    @Override
    public int getLength() {
        return bases.length;
    }

    @Override
    public void setPosition(final String contig, final int start) {
        throw new UnsupportedOperationException();
    }

    @Override
    public void setPosition(final Locatable locatable) {
        throw new UnsupportedOperationException();
    }

    @Override
    public int getUnclippedStart() {
        return getStart();
    }

    @Override
    public int getUnclippedEnd() {
        return getEnd();
    }

    @Override
    public String getMateContig() {
        return null;
    }

    @Override
    public int getMateStart() {
        return 0;
    }

    @Override
    public void setMatePosition(final String contig, final int start) {
        throw new UnsupportedOperationException();
    }

    @Override
    public void setMatePosition(final Locatable locatable) {
        throw new UnsupportedOperationException();
    }

    @Override
    public int getFragmentLength() {
        return 0;
    }

    @Override
    public void setFragmentLength(final int fragmentLength) {
        throw new UnsupportedOperationException();
    }

    @Override
    public int getMappingQuality() {
        return mappingQuality;
    }

    @Override
    public void setMappingQuality(final int mappingQuality) {
        throw new UnsupportedOperationException();
    }

    @Override
    public byte[] getBases() {
        return bases;
    }

    @Override
    public void setBases(final byte[] bases) {
        throw new UnsupportedOperationException();
    }

    @Override
    public byte[] getBaseQualities() {
        return qualities;
    }

    @Override
    public void setBaseQualities(final byte[] baseQualities) {
        throw new UnsupportedOperationException();
    }

    @Override
    public Cigar getCigar() {
        return cigar;
    }

    @Override
    public void setCigar(final Cigar cigar) {
        throw new UnsupportedOperationException();
    }

    @Override
    public void setCigar(final String cigarString) {
        throw new UnsupportedOperationException();
    }

    @Override
    public String getReadGroup() {
        return readGroup;
    }

    @Override
    public void setReadGroup(final String readGroupID) {
        throw new UnsupportedOperationException();
    }

    @Override
    public boolean isPaired() {
        return false;
    }

    @Override
    public void setIsPaired(final boolean isPaired) {
        throw new UnsupportedOperationException();
    }

    @Override
    public boolean isProperlyPaired() {
        return false;
    }

    @Override
    public void setIsProperlyPaired(final boolean isProperlyPaired) {
        throw new UnsupportedOperationException();
    }

    @Override
    public boolean isUnmapped() {
        return false;
    }

    @Override
    public void setIsUnmapped() {
        throw new UnsupportedOperationException();
    }

    @Override
    public boolean mateIsUnmapped() {
        return false;
    }

    @Override
    public void setMateIsUnmapped() {
        throw new UnsupportedOperationException();
    }

    @Override
    public boolean isReverseStrand() {
        return false;
    }

    @Override
    public void setIsReverseStrand(final boolean isReverseStrand) {
        throw new UnsupportedOperationException();
    }

    @Override
    public boolean mateIsReverseStrand() {
        return false;
    }

    @Override
    public void setMateIsReverseStrand(final boolean mateIsReverseStrand) {
        throw new UnsupportedOperationException();
    }

    @Override
    public boolean isFirstOfPair() {
        return false;
    }

    @Override
    public void setIsFirstOfPair() {
        throw new UnsupportedOperationException();
    }

    @Override
    public boolean isSecondOfPair() {
        return false;
    }

    @Override
    public void setIsSecondOfPair() {
        throw new UnsupportedOperationException();
    }

    @Override
    public boolean isSecondaryAlignment() {
        return false;
    }

    @Override
    public void setIsSecondaryAlignment(final boolean isSecondaryAlignment) {
        throw new UnsupportedOperationException();
    }

    @Override
    public boolean isSupplementaryAlignment() {
        return false;
    }

    @Override
    public void setIsSupplementaryAlignment(final boolean isSupplementaryAlignment) {
        throw new UnsupportedOperationException();
    }

    @Override
    public boolean failsVendorQualityCheck() {
        return false;
    }

    @Override
    public void setFailsVendorQualityCheck(final boolean failsVendorQualityCheck) {
        throw new UnsupportedOperationException();
    }

    @Override
    public boolean isDuplicate() {
        return false;
    }

    @Override
    public void setIsDuplicate(final boolean isDuplicate) {
        throw new UnsupportedOperationException();
    }

    @Override
    public boolean hasAttribute(final String attributeName) {
        throw new UnsupportedOperationException();
    }

    @Override
    public Integer getAttributeAsInteger(final String attributeName) {
        throw new UnsupportedOperationException();
    }

    @Override
    public String getAttributeAsString(final String attributeName) {
        throw new UnsupportedOperationException();
    }

    @Override
    public byte[] getAttributeAsByteArray(final String attributeName) {
        throw new UnsupportedOperationException();
    }

    @Override
    public void setAttribute(final String attributeName, final Integer attributeValue) {
        throw new UnsupportedOperationException();
    }

    @Override
    public void setAttribute(final String attributeName, final String attributeValue) {
        throw new UnsupportedOperationException();
    }

    @Override
    public void setAttribute(final String attributeName, final byte[] attributeValue) {
        throw new UnsupportedOperationException();
    }

    @Override
    public void clearAttribute(final String attributeName) {
        throw new UnsupportedOperationException();
    }

    @Override
    public void clearAttributes() {
        throw new UnsupportedOperationException();
    }

    @Override
    public GATKRead copy() {
        return new MergedGATKReadPair(name, readGroup, interval, bases, qualities, cigar, mappingQuality);
    }

    @Override
    public GATKRead deepCopy() {
        return new MergedGATKReadPair(name, readGroup, interval, bases.clone(), qualities.clone(), cigar, mappingQuality);
    }

    @Override
    public SAMRecord convertToSAMRecord(final SAMFileHeader header) {
        throw new UnsupportedOperationException();
    }

    @Override
    public Read convertToGoogleGenomicsRead() {
        throw new UnsupportedOperationException();
    }

    @Override
    public String getSAMString() {
        throw new UnsupportedOperationException();
    }

    @Override
    public String getContig() {
        return interval.getContig();
    }

    @Override
    public String getAssignedContig() {
        return getContig();
    }

    @Override
    public int getStart() {
        return interval.getStart();
    }

    @Override
    public int getAssignedStart() {
        return getStart();
    }

    @Override
    public int getEnd() {
        return interval.getEnd();
    }

    /**
     * An alignment block is a contiguous section of a read that aligns with
     * a contiguous region on the reference.
     */
    private static class AlignmentBlock {

        /**
         * Whether the calls on this block have preference over other overlapping calls
         * with different base but same quality.
         */
        public boolean priority;


        /**
         * bases of the read the alignment block refers to.
         */
        public final byte[] bases;

        /**
         * Qualities scores of the read the alignment block refers to.
         */
        public final byte[] quals;

        /**
         * Offset of the first base in the read this alignment block covers.
         * This is a 0-based index.
         */
        public final int readStart;

        /**
         * Offset of the read base just after the last one covered by this alignment block.
         */
        public final int readEnd;

        /**
         * Offset of the first base in the reference this alignment cover.
         * This is a 1-based index.
         */
        public final int referenceStart;

        /**
         * Offset of the reference base right after the last covered by this alignment block.
         * This is a 1-based index.
         */
        public final int referenceEnd;

        /**
         * Length of the alignment block in base-pairs.
         */
        public final int length;

        /**
         * Creates a new alignment block given the corresponding read,
         * @param read the alignment blocks read.
         * @param readOffset the offset in the read for the first base/qual in the alignment block. This index is 0-based.
         * @param referenceOffset the offset in the reference sequence. This index is 1 based.
         * @param length the length of the alignment block.
         * @param priority whether the block is a priority call block (see {@link #priority}).
         * @throws IllegalArgumentException if {@code read} is {@code null}, the value given to other arguments are
         * inconsistent with the length of the read or {@code referenceOffset} is not a valid 1-based contig coordinate value.
         */
        public AlignmentBlock(final GATKRead read, final int readOffset, final int referenceOffset, final int length, final boolean priority) {
            this(Utils.nonNull(read).getBases(), read.getBaseQualities(), readOffset, referenceOffset, length, priority);
        }

        /**
         * Composes an alignment block given base and quality arrays.
         *
         * @param bases the bases from the read the alignment block makes reference to.
         * @param quals the quals from the read the alignment block makes reference to.
         * @param readOffset the offset in the read for the first base/qual in the alignment block. This index is 0-based.
         * @param referenceOffset the offset in the reference sequence. This index is 1 based.
         * @param length the length of the alignment block.
         * @param priority whether it is a priority block.
         */
        public AlignmentBlock(final byte[] bases, final byte[] quals, final int readOffset, final int referenceOffset, final int length, final boolean priority) {
            this.priority = priority;
            this.length = ParamUtils.isPositive(length, "the input length must be greater than 0");
            this.bases = Utils.nonNull(bases, "input bases array must not be null");
            this.quals = Utils.nonNull(quals, "input quals array must not be null");
            Utils.validateArg(this.bases.length == this.quals.length, "input base and quals array must have the same length");
            readStart = ParamUtils.inRange(readOffset, 0, bases.length - length, "the input read-offset is invalid");
            readEnd = readOffset + length;
            referenceStart = ParamUtils.isPositive(referenceOffset, "the reference offset must be greater than 0");
            referenceEnd = referenceOffset + length;
        }

        /**
         * Returns the cigar element that would represent this alignment block.
         * @return never {@code null}.
         */
        public CigarElement toCigarElement() {
            return new CigarElement(length, CigarOperator.MATCH_OR_MISMATCH);
        }
    }

    private static SimpleInterval calculateMatchInterval(final String contig, final List<AlignmentBlock> relevantCigarElements) {
        final int referenceStart = relevantCigarElements.get(0).referenceStart;
        final AlignmentBlock lastRelevantCigarElement = relevantCigarElements.get(relevantCigarElements.size() - 1);
        final int referenceStop = lastRelevantCigarElement.referenceStart
                + lastRelevantCigarElement.length - 1;
        return new SimpleInterval(contig, referenceStart, referenceStop);
    }

    private static List<AlignmentBlock> extractAlignmentBlocks(final GATKRead read, final boolean priority) {
        int readStart = 0;
        int referenceStart = read.getStart();
        final List<AlignmentBlock> result = new ArrayList<>(read.getCigar().numCigarElements());
        for (final CigarElement element : read.getCigar().getCigarElements()) {
            final boolean consumesReadBases = element.getOperator().consumesReadBases();
            final boolean consumesReferenceBases = element.getOperator().consumesReferenceBases();
            if (consumesReadBases && consumesReferenceBases && element.getLength() > 0) {
                result.add(new AlignmentBlock(read, readStart, referenceStart, element.getLength(), priority));
            }
            readStart += consumesReadBases ? element.getLength() : 0;
            referenceStart += consumesReferenceBases ? element.getLength() : 0;
        }
        return result;
    }
}
