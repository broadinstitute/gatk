package org.broadinstitute.hellbender.utils;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import it.unimi.dsi.fastutil.ints.*;
import org.broadinstitute.hellbender.utils.read.CigarUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;

import java.util.*;
import java.util.stream.Collectors;

/**
 * An interval map where bases and qualities are explicitly stored in matrices
 * for fast look-up.
 */
class ByteMapIntervalPileup implements IntervalPileup {

    private final ReferenceBases referenceBases;
    private final List<GATKRead> reads;
    private final List<IntervalPileup.Element> elements;
    private final Int2ObjectMap<IntervalPileup.Element> elementByIndex;
    private final Map<GATKRead, IntervalPileup.Element> elementByRead;
    private final int width;
    private final int height;
    private byte[][] bases;
    private byte[][] quals;
    private final List<Insert> insertsBuffer = new ArrayList<>(10);
    private final IntList insertsBufferOffsets = new IntArrayList(10);


    ByteMapIntervalPileup(final ReferenceBases referenceBases, final List<GATKRead> reads) {
        this.referenceBases = referenceBases;
        this.width = referenceBases.getInterval().size();
        final SimpleInterval referenceInterval = referenceBases.getInterval();
        final int referenceStart = referenceInterval.getStart();
        final int referenceEnd = referenceInterval.getEnd();
        this.reads = Collections.unmodifiableList(reads.stream()
                .filter(read -> !read.isUnmapped() && read.getStart() <= referenceEnd && read.getEnd() >= referenceStart)
                .sorted(Comparator.comparingInt(GATKRead::getStart).thenComparing(GATKRead::getEnd).thenComparing(GATKRead::getName))
                .collect(Collectors.toList()));
        this.elements = new ArrayList<>(this.reads.size());
        this.elementByIndex = new Int2ObjectOpenHashMap<>(this.reads.size());
        this.elementByRead = new HashMap<>(this.reads.size());
        this.height = this.reads.size();
        bases = new byte[height][width];
        quals = new byte[height][width];
        for (int i = 0; i < height; i++) {
            final IntervalPileup.Element element = new Element(this, this.reads.get(i), i);
            elements.add(element);
            elementByIndex.put(i, element);
            elementByRead.put(this.reads.get(i), element);
        }
    }

    @Override
    public List<GATKRead> reads() {
        return reads;
    }

    @Override
    public ReferenceBases reference() {
        return referenceBases;
    }

    @Override
    public int width() {
        return width;
    }

    @Override
    public int height() {
        return height;
    }

    @Override
    public byte baseAt(final int row, final int column) {
        return bases[row][column];
    }

    @Override
    public byte qualAt(int row, int column) {
        return quals[row][column];
    }

    @Override
    public IntervalPileup.Insert insertAt(final int row, final int column) {
        final IntervalPileup.Element element = elementByIndex.get(row);
        return element != null ? element.insertAt(column) : null;
    }

    @Override
    public GATKRead readAt(int row, int column) {
        final IntervalPileup.Element element = elements.get(row);
        if (element.minColumn() > column || element.maxColumn() < column) {
            return null;
        } else {
            return element.read();
        }
    }

    @Override
    public List<GATKRead> readsAt(int row, int column) {
        final GATKRead read = readAt(row, column);
        return read != null ? Collections.singletonList(read) : Collections.emptyList();
    }

    @Override
    public IntervalPileup.Element element(final GATKRead read) {
        return elementByRead.get(read);
    }

    @Override
    public boolean hasInsertAt(final int row, final int column) {
        final IntervalPileup.Element element = elements.get(row);
        return element.hasInsertAt(column);
    }

    @Override
    public byte[] insertBasesAt(final int row, final int column) {
        final IntervalPileup.Element element = elements.get(row);
        return element.insertBasesAt(column);
    }

    @Override
    public byte[] insertQualsAt(final int row, final int column) {
        final IntervalPileup.Element element = elements.get(row);
        return element.insertQualsAt(column);
    }

    static class Insert implements  IntervalPileup.Insert {
        private final GATKRead enclosingRead;
        private final int offset;
        private final int length;
        private final int column;
        private transient int hashCode = 0;

        Insert(final GATKRead enclosingRead, final int offset, final int column, final int length) {
            this.enclosingRead = enclosingRead;
            this.offset = offset;
            this.length = length;
            this.column = column;
        }

        public int column() {
            return column;
        }

        @Override
        public int length() {
            return length;
        }

        @Override
        public byte[] bases() {
            final byte[] result = new byte[length];
            final int copied = enclosingRead.copyBases(offset, result, 0, length);
            if (copied < length) {
                Arrays.fill(result, copied, result.length, (byte) 'N');
            }
            return result;
        }

        @Override
        public byte[] quals() {
            final byte[] result = new byte[length];
            final int copied = enclosingRead.copyBaseQualities(offset, result, 0, length);
            if (copied < length) {
                Arrays.fill(result, copied, result.length, NO_BQ);
            }
            return result;
        }

        @Override
        public int copyBases(int offset, byte[] dest, int destOffset, int maxLength) {
            final int actualMaxLength = Math.min(length, maxLength);
            final int copied = enclosingRead.copyBases(this.offset + offset, dest, destOffset, actualMaxLength);
            if (copied < actualMaxLength) {
                Arrays.fill(dest, destOffset + copied, destOffset + actualMaxLength, (byte) 'N');
            }
            return actualMaxLength;
        }

        @Override
        public int copyQuals(int offset, byte[] dest, int destOffset, int maxLength) {
            final int actualMaxLength = Math.min(length, maxLength);
            final int copied = enclosingRead.copyBaseQualities(this.offset + offset, dest, destOffset, actualMaxLength);
            if (copied < actualMaxLength) {
                Arrays.fill(dest, destOffset + copied, destOffset + actualMaxLength, (byte) 'N');
            }
            return actualMaxLength;
        }

        @Override
        public int hashCode() {
            if (hashCode == 0) {
                final byte[] readBases = enclosingRead.getBasesNoCopy();
                if (readBases != null && readBases.length > offset) {
                    final int to = offset + length;
                    final int to2 = to <= readBases.length ? to : readBases.length;
                    int i;
                    hashCode = 1;
                    for (i = offset; i < to2; i++) {
                        hashCode = hashCode * 31 + readBases[i];
                    }
                    for (; i < to; i++) {
                        hashCode = hashCode * 31 + 'N';
                    }
                }
            }
            return hashCode;
        }

        @Override
        public boolean equals(final Object other) {
            return other == this || (other instanceof Insert && equals((IntervalPileup.Insert) other));
        }

        private boolean equals(final IntervalPileup.Insert other) {
            return length == other.length() && hashCode() == other.hashCode() &&
                    equalBases(other) && equalQualities(other);
        }

        private boolean equalBases(final IntervalPileup.Insert other) {
            if (other instanceof Insert) {
                equalBases((Insert) other);
            }
            final byte[] otherBases = other.bases();
            final byte[] readBases = enclosingRead.getBasesNoCopy();
            for (int i = 0; i < otherBases.length; i++) {
                if (otherBases[i] != readBases[offset + i]) {
                    return false;
                }
            }
            return true;
        }

        private boolean equalBases(final Insert other) {
            final byte[] otherBases = other.enclosingRead.getBasesNoCopy();
            final byte[] thisBases = enclosingRead.getBasesNoCopy();
            for (int i = 0; i < length; i++) {
                if (otherBases[other.offset + i] != thisBases[offset + i]) {
                    return false;
                }
            }
            return true;
        }

        private boolean equalQualities(final IntervalPileup.Insert other) {
            if (other instanceof Insert) {
                equalQualities((Insert) other);
            }
            final byte[] otherQuals = other.quals();
            final byte[] thisQuals = enclosingRead.getBaseQualitiesNoCopy();
            for (int i = 0; i < otherQuals.length; i++) {
                if (otherQuals[i] != thisQuals[offset + i]) {
                    return false;
                }
            }
            return true;
        }

        private boolean equalQualities(final Insert other) {
            final byte[] otherQuals = other.enclosingRead.getBaseQualitiesNoCopy();
            final byte[] thisQuals = enclosingRead.getBaseQualitiesNoCopy();
            for (int i = 0; i < length; i++) {
                if (otherQuals[other.offset + i] != thisQuals[offset + i]) {
                    return false;
                }
            }
            return true;
        }

        @Override
        public String toString() {
            return new String(bases()) + "/" + new String(quals());
        }
    }

    static class Element implements IntervalPileup.Element {

        private final GATKRead read;
        private final int row;
        private final int minColumn;
        private final int maxColumn;
        private final Int2ObjectMap<Insert> inserts;
        private final ByteMapIntervalPileup pileup;

        private Element(final ByteMapIntervalPileup pileup, final GATKRead read, final int row) {
            this.pileup = pileup;
            this.row = row;
            this.read = read;
            final Cigar cigar = read.getCigar();
            final int referenceWidth = cigar.getReferenceLength();
            int referenceOffset = read.getStart() - pileup.referenceBases.getInterval().getStart();
            int readOffset = 0;
            minColumn = Math.max(0, referenceOffset);
            maxColumn = Math.min(pileup.width - 1, referenceOffset + referenceWidth - 1);
            Arrays.fill(pileup.bases[row], 0, minColumn, NO_BASE);
            Arrays.fill(pileup.bases[row], maxColumn + 1, pileup.width, NO_BASE);
            Arrays.fill(pileup.quals[row], 0, minColumn, NO_BQ);
            Arrays.fill(pileup.quals[row], maxColumn + 1, pileup.width, NO_BQ);
            final List<CigarElement> cigarElements = cigar.getCigarElements();
            pileup.insertsBuffer.clear();
            pileup.insertsBufferOffsets.clear();
            int i;
            for (i = 0; referenceOffset <= maxColumn + 1 &&  i < cigarElements.size(); i++) {
                final CigarElement element = cigarElements.get(i);
                final CigarOperator op = element.getOperator();
                final int length = element.getLength();
                final int newReferenceOffset = referenceOffset + (op.consumesReferenceBases() ? length : 0);
                final int newReadOffset = readOffset + (op.consumesReadBases() ? length : 0);
                if (newReferenceOffset >= minColumn) {
                    if (op.isAlignment()) {
                        final int leftOverhang = referenceOffset < minColumn ? minColumn - referenceOffset : 0;
                        final int from = referenceOffset + leftOverhang;
                        final int len = newReferenceOffset > maxColumn ? maxColumn + 1 - from :
                                  length - leftOverhang;
                        final int readFrom = readOffset + leftOverhang;
                        final int copiedBases = read.copyBases(readFrom, pileup.bases[row], from, len);
                        final int copiedQuals = read.copyBaseQualities(readFrom, pileup.quals[row], from, len);
                        if (copiedBases < len) {
                            Arrays.fill(pileup.bases[row], from + copiedBases, from + len - copiedBases, (byte) 'N');
                        }
                        if (copiedQuals < len) {
                            Arrays.fill(pileup.quals[row], from + copiedQuals, from + len - copiedQuals, NO_BQ);
                        }
                    } else if (op == CigarOperator.I) {
                        pileup.insertsBuffer.add(new Insert(read, readOffset, referenceOffset - 1, length));
                        pileup.insertsBufferOffsets.add(referenceOffset - 1);
                    } else if (op == CigarOperator.D || op == CigarOperator.N) {
                        final int from = referenceOffset < minColumn ? minColumn : referenceOffset;
                        final int len = (newReferenceOffset > maxColumn ? maxColumn + 1 : newReferenceOffset) - from;
                        Arrays.fill(pileup.bases[row], from, from + len, GAP);
                        Arrays.fill(pileup.quals[row], from, from + len, NO_BQ);
                    }
                }
                readOffset = newReadOffset;
                referenceOffset = newReferenceOffset;
            }
            // merge adjacent inserts if any.
            mergeAdjacentInserts(read, pileup.insertsBuffer, pileup.insertsBufferOffsets);
            inserts = consolidateInserts(pileup.insertsBuffer, pileup.insertsBufferOffsets);
        }

        private static Int2ObjectMap<Insert> consolidateInserts(final List<Insert> inserts, final IntList offsets) {
            final int size = inserts.size();
            if (size == 0) {
                return Int2ObjectMaps.emptyMap();
            } else if (size == 1) {
                return Int2ObjectMaps.singleton(offsets.get(0) , inserts.get(0));
            } else {
                final Int2ObjectMap<Insert> result = size < 5 ? new Int2ObjectArrayMap<>(size)
                        : new Int2ObjectOpenHashMap<>(size);
                for (int ins = 0; ins < size; ins++) {
                    result.put(offsets.get(ins) , inserts.get(ins));
                }
                return result;
            }
        }

        private static void mergeAdjacentInserts(final GATKRead read, final List<Insert> inserts, final IntList offsets) {
            int i;
            for (i = inserts.size() - 1; i > 0; i--) {
                if (offsets.get(i) == offsets.get(i-1)) {
                    final Insert left = inserts.get(i-1);
                    final Insert right = inserts.get(i);
                    final Insert merged = new Insert(read, left.offset, left.column,left.length + right.length);
                    inserts.remove(i);
                    inserts.set(i - 1, merged);
                    offsets.remove(i);
                }
            }
        }

        @Override
        public GATKRead read() {
            return read;
        }

        @Override
        public int row() {
            return row;
        }

        @Override
        public int minColumn() {
            return minColumn;
        }

        @Override
        public int maxColumn() {
            return maxColumn;
        }


        @Override
        public boolean hasInsertAt(final int column) {
            return inserts.containsKey(column);
        }

        @Override
        public int insertSize(final int column) {
            final Insert insert = inserts.get(column);
            return insert != null ? insert.length() : 0;
        }

        @Override
        public int copyInsertQuals(final int column, final byte[] dest, final int offset, final int maxLength) {
            Utils.nonNull(dest);
            Utils.validIndex(offset, dest.length);
            Utils.validIndex(offset + maxLength, dest.length);
            final Insert insert = inserts.get(column);
            if (maxLength != 0 && insert != null) {
                return insert.copyQuals(0, dest, offset, maxLength);
            } else {
                return 0;
            }
        }

        @Override
        public int copyInsertBases(final int column, final byte[] dest, final int offset, final int maxLength) {
            Utils.nonNull(dest);
            Utils.validIndex(offset, dest.length);
            Utils.validIndex(offset + maxLength, dest.length);
            final Insert insert = inserts.get(column);
            if (maxLength != 0 && insert != null) {
                return insert.copyBases(0, dest, offset, maxLength);
            } else {
                return 0;
            }
        }

        @Override
        public byte[] insertBasesAt(final int column) {
            final Insert insert = inserts.get(column);
            return insert != null ? insert.bases() : null;
        }

        @Override
        public byte[] insertQualsAt(final int column) {
            final Insert insert = inserts.get(column);
            return insert != null ? insert.quals() : null;
        }

        @Override
        public Insert insertAt(final int column) {
            return inserts.get(column);
        }

        @Override
        public List<IntervalPileup.Insert> inserts() {
            return new ArrayList<>(inserts.values());
        }

        @Override
        public byte baseAt(final int column) {
            return pileup.bases[row][column];
        }

        @Override
        public byte qualAt(final int column) {
            return pileup.quals[row][column];
        }
    }
}
