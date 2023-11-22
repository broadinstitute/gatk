package org.broadinstitute.hellbender.tools.walkers.softcliprealignment;

import com.google.common.collect.Lists;
import htsjdk.samtools.*;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.MergingIterator;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.engine.ReadsPathDataSource;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.bwa.BwaReadAligner;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.codehaus.plexus.util.FileUtils;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.function.Predicate;

/**
 * Multipurpose class for performing read alignment on a subset of reads using BWA. Automatically detects paired and
 * unpaired reads using read names (flags are ignored for this) and uses the appropriate alignment mode for each.
 * It merges the realigned subset ("selected" reads) with the rest of the original reads ("non-selected" reads)
 * and returns them in coordinate sorted order.
 *
 * This class is used as follows. First, each read is passed in with the {@link #addRead(GATKRead, Predicate)}
 * method, which requires one to specify the predicate for whether that read should be realigned. The reads must
 * be provided in coordinate sorted order. Once all reads have been added, a single call to {@link #alignAndMerge()}
 * can be made to perform the alignment and merge with the non-realigned reads.
 *
 * To ensure a consistently small memory footprint, reads are cached to disk using temporary BAM files and buffered
 * for realignment. Also remember to call {@link #close()} or use try-with-resources to ensure that cached data is
 * properly cleaned up.
 *
 * Note that instances of this class are NOT idempotent because calls to {@link #addRead(GATKRead, Predicate)} and
 * {@link #alignAndMerge()} must be made in this specific order.
 */
public final class SubsettingRealignmentEngine implements AutoCloseable {

    // Read buffers
    private GATKRead lastRead = null;
    private List<GATKRead> pairedBuffer;
    private List<GATKRead> unpairedBuffer;
    private final int bufferSize;

    // Read headers
    private final SAMFileHeader inputHeader;  // coordinate sorted (input and final output)
    private final SAMFileHeader alignmentHeader; // same as input header but queryname sorted

    // Aligners
    private final BwaReadAligner pairedAligner;
    private final BwaReadAligner unpairedAligner;

    // Temporary files for the first pass, which divides reads into "selected" and "non-selected" groups
    private File selectedReadsBam;
    private File nonselectedReadsBam;
    private SAMFileWriter selectedReadsWriter;
    private SAMFileWriter nonselectedReadsWriter;

    // Iterators for the final merge of non-selected and aligned reads
    private CloseableIterator<GATKRead> unselectedRecordsIter;
    private CloseableIterator<GATKRead> alignedRecordsIter;

    // Read counters, for logging
    private long selectedReadsCount = 0;
    private long nonselectedReadsCount = 0;
    private long pairedAlignmentReadsCount = 0;
    private long unpairedAlignmentReadsCount = 0;

    // Tag assigned to reads that were realigned
    public static final String REALIGNED_READ_TAG = "RA";
    public static final int REALIGNED_READ_TAG_VALUE = 1;

    /**
     * Constructor
     * @param indexImagePath path to image created with {@link org.broadinstitute.hellbender.tools.BwaMemIndexImageCreator}
     * @param inputHeader    input reads header (coordinate-sorted)
     * @param bufferSize     number of reads for the alignment buffers
     * @param bwaThreads     number of bwa threads
     * @param retainDuplicateFlag  keep duplicate flag for each realigned read
     */
    public SubsettingRealignmentEngine(final String indexImagePath, final SAMFileHeader inputHeader,
                                       final int bufferSize, final int bwaThreads, final boolean retainDuplicateFlag) {
        if (!inputHeader.getSortOrder().equals(SAMFileHeader.SortOrder.coordinate)) {
            throw new UserException.BadInput("Input is not coordinate sorted");
        }
        this.bufferSize = bufferSize;
        pairedBuffer = new ArrayList<>(bufferSize);
        unpairedBuffer = new ArrayList<>(bufferSize);
        this.inputHeader = inputHeader;

        // The input must be queryname sorted
        alignmentHeader = inputHeader.clone();
        alignmentHeader.setSortOrder(SAMFileHeader.SortOrder.queryname);

        pairedAligner = new BwaReadAligner(indexImagePath, alignmentHeader, true, retainDuplicateFlag);
        pairedAligner.setSoftClipSupplementaryAlignments(true);  // as in the Warp pipeline
        pairedAligner.setNThreads(bwaThreads);

        unpairedAligner = new BwaReadAligner(indexImagePath, alignmentHeader, false, retainDuplicateFlag);
        pairedAligner.setSoftClipSupplementaryAlignments(true);
        unpairedAligner.setNThreads(bwaThreads);

        selectedReadsBam = FileUtils.createTempFile("realign_1", FileExtensions.BAM, null);
        nonselectedReadsBam = FileUtils.createTempFile("realign_2", FileExtensions.BAM, null);
        selectedReadsWriter = new SAMFileWriterFactory().makeWriter(alignmentHeader, false, selectedReadsBam, null);
        nonselectedReadsWriter = new SAMFileWriterFactory().makeWriter(inputHeader, true, nonselectedReadsBam, null);
    }

    /**
     * Cleans up temporary files and writers
     */
    @Override
    public void close() {
        closeFirstPassWriters();
        cleanUpFile(selectedReadsBam);
        selectedReadsBam = null;
        cleanUpFile(nonselectedReadsBam);
        nonselectedReadsBam = null;
        if (unselectedRecordsIter != null) {
            unselectedRecordsIter.close();
            unselectedRecordsIter = null;
        }
        if (alignedRecordsIter != null) {
            alignedRecordsIter.close();
            alignedRecordsIter = null;
        }
    }

    private void closeFirstPassWriters() {
        if (selectedReadsWriter != null) {
            selectedReadsWriter.close();
            selectedReadsWriter = null;
        }
        if (nonselectedReadsWriter != null) {
            nonselectedReadsWriter.close();
            nonselectedReadsWriter = null;
        }
    }

    /**
     * For deleting temporary files
     */
    private void cleanUpFile(final File file) {
        if (file == null) {
            return;
        }
        try {
            FileUtils.forceDelete(file);
        } catch (final IOException e) {
            throw new GATKException("Could not clean up temporary file", e);
        }
    }

    /**
     * Add an input read. Note that any {@link #REALIGNED_READ_TAG} tags should be cleared prior. All reads are
     * cached on disk in BAM format. This method may NOT be called after any call to @{link #alignAndMerge} or
     * @{link #close}. Also, it is assumed that mates have the same name (i.e. no "/1" and "/2") to detect pairs.
     * @param read  candidate realignment read
     * @param predicate  returns true iff the read should be realigned
     */
    public void addRead(final GATKRead read, final Predicate<GATKRead> predicate) {
        Utils.nonNull(selectedReadsWriter, "This instance has been closed");
        Utils.nonNull(nonselectedReadsWriter, "This instance has been closed");
        if (predicate.test(read)) {
            if (ReadFilterLibrary.PRIMARY_LINE.test(read)) {
                selectedReadsWriter.addAlignment(read.convertToSAMRecord(selectedReadsWriter.getFileHeader()));
                selectedReadsCount++;
            }
        } else {
            // TODO add back in?
            //nonselectedReadsWriter.addAlignment(read.convertToSAMRecord(nonselectedReadsWriter.getFileHeader()));
            //nonselectedReadsCount++;
        }
    }

    public void addDistantMate(final GATKRead mate) {
        Utils.nonNull(selectedReadsWriter, "This instance has been closed");
        Utils.nonNull(nonselectedReadsWriter, "This instance has been closed");
        if (ReadFilterLibrary.PRIMARY_LINE.test(mate)) {
            selectedReadsWriter.addAlignment(mate.convertToSAMRecord(selectedReadsWriter.getFileHeader()));
            selectedReadsCount++;
        }
    }

    /**
     * Performs realignment on the reads submitted with @{@link #addRead(GATKRead, Predicate)}, and merges
     * the alignments with non-selected reads. Caches realignments to disk before merging. This method may NOT be
     * called after a previous call to this method or @{link #close}.
     * @return merged iterator of realigned and non-selected reads, all coordinate-sorted
     */
    public MergingIterator<GATKRead> alignAndMerge() {
        Utils.nonNull(selectedReadsBam, "This instance has been run already");
        Utils.nonNull(nonselectedReadsBam, "This instance has been closed");
        // Assume we're done with adding reads, so we can close the first set of temp files
        closeFirstPassWriters();
        
        // Align selected reads
        final File alignmentBam = FileUtils.createTempFile("realign_3", FileExtensions.BAM, null);
        try (final ReadsDataSource selectedSource = new ReadsPathDataSource(selectedReadsBam.toPath());
             final SAMFileWriter alignmentWriter = new SAMFileWriterFactory().makeWriter(inputHeader, false, alignmentBam, null)) {
            for (final GATKRead read : selectedSource) {
                writeAlignments(addReadToBuffer(read), alignmentWriter);
            }
            // Flush leftover reads
            writeAlignments(flushPairedBuffer(), alignmentWriter);
            addLastReadToBuffer();
            writeAlignments(flushUnpairedBuffer(), alignmentWriter);
        }

        // Delete unneeded temp file of the original selected reads
        cleanUpFile(selectedReadsBam);
        selectedReadsBam = null;

        // Merge the non-selected reads with the new alignments
        final ReadsDataSource unselectedSource = new ReadsPathDataSource(nonselectedReadsBam.toPath());
        final ReadsDataSource alignmentSource = new ReadsPathDataSource(alignmentBam.toPath());
        unselectedRecordsIter = new ClosableReadSourceIterator(unselectedSource);
        alignedRecordsIter = new ClosableReadSourceIterator(alignmentSource);
        // htsjdk read comparators attempt tie-breaking after comparing coordinates, but this is generally too
        // strict for data in the wild, so we won't enforce it here.
        final Comparator<GATKRead> weakComparator = new Comparator<>() {
            final SAMRecordComparator strongComparator = inputHeader.getSortOrder().getComparatorInstance();
            @Override
            public int compare(GATKRead o1, GATKRead o2) {
                return strongComparator.fileOrderCompare(o1.convertToSAMRecord(inputHeader), o2.convertToSAMRecord(inputHeader));
            }
        };
        return new MergingIterator<>(weakComparator, List.of(unselectedRecordsIter, alignedRecordsIter));
    }

    /**
     * Simple helper method for writing alignment output
     */
    private void writeAlignments(final List<GATKRead> reads, final SAMFileWriter writer) {
        for (final GATKRead read: reads) {
            writer.addAlignment(read.convertToSAMRecord(inputHeader));
        }
    }

    /**
     * Realigns the given reads and properly handles pairedness. Note that the buffers and {@link #lastRead} may
     * not be empty/null at the end and must be flushed afterward.
     */
    private List<GATKRead> addReadToBuffer(final GATKRead read) {
        if (lastRead == null) {
            lastRead = read;
        } else if (lastRead.getName().equals(read.getName())) {
            pairedBuffer.add(lastRead);
            pairedBuffer.add(read);
            lastRead = null;
            if (pairedBuffer.size() == bufferSize || pairedBuffer.size() + 1 == bufferSize) {
                return flushPairedBuffer();
            }
        } else {
            // Names don't match, so the previous read was unpaired
            unpairedBuffer.add(lastRead);
            lastRead = read;
            return flushUnpairedBuffer();
        }
        return Collections.emptyList();
    }
    private List<GATKRead> flushPairedBuffer() {
        pairedAlignmentReadsCount += pairedBuffer.size();
        return flushBuffer(pairedBuffer, pairedAligner);
    }

    private List<GATKRead> flushUnpairedBuffer() {
        unpairedAlignmentReadsCount += unpairedBuffer.size();
        return flushBuffer(unpairedBuffer, unpairedAligner);
    }

    private void addLastReadToBuffer() {
        if (lastRead != null) {
            unpairedBuffer.add(lastRead);
            lastRead = null;
        }
    }

    /**
     * Realign all reads in the given buffer and clear it
     */
    private List<GATKRead> flushBuffer(final List<GATKRead> mutableBuffer, final BwaReadAligner aligner) {
        final List<GATKRead> alignments = Lists.newArrayList(aligner.apply(mutableBuffer));
        alignments.forEach(SubsettingRealignmentEngine::setRealignedTag);
        mutableBuffer.clear();
        return alignments;
    }

    /**
     * Sets tag {@link #REALIGNED_READ_TAG} to mark it as realigned
     */
    private static void setRealignedTag(final GATKRead read) {
        read.setAttribute(REALIGNED_READ_TAG, REALIGNED_READ_TAG_VALUE);
    }

    /**
     * Returns true if the read was realigned
     */
    public static boolean readIsRealigned(final GATKRead read) {
        return read.hasAttribute(REALIGNED_READ_TAG) && read.getAttributeAsInteger(REALIGNED_READ_TAG) == REALIGNED_READ_TAG_VALUE;
    }

    /**
     * Number of reads selected for realignment
     */
    public long getSelectedReadsCount() {
        return selectedReadsCount;
    }

    /**
     * Number of reads not selected for realignment
     */
    public long getNonselectedReadsCount() {
        return nonselectedReadsCount;
    }

    /**
     * Number of realigned paired reads
     */
    public long getPairedAlignmentReadsCount() {
        return pairedAlignmentReadsCount;
    }

    /**
     * Number of realigned unpaired reads
     */
    public long getUnpairedAlignmentReadsCount() {
        return unpairedAlignmentReadsCount;
    }

    /**
     * Helper class to map a {@link ReadsDataSource} iterator to a {@link CloseableIterator} and each input read to a
     * {@link SAMRecord} as well.
     */
    static final class ClosableReadSourceIterator implements CloseableIterator<GATKRead> {

        ReadsDataSource source;
        final Iterator<GATKRead> iterator;

        public ClosableReadSourceIterator(final ReadsDataSource source) {
            this.source = source;
            this.iterator = source.iterator();
        }

        @Override
        public void close() {
            if (source != null) {
                source.close();
                source = null;
            }
        }

        @Override
        public boolean hasNext() {
            if (!iterator.hasNext()) {
                close();
                return false;
            }
            return iterator.hasNext();
        }

        @Override
        public GATKRead next() {
            if (source == null) {
                throw new NoSuchElementException("Source has been closed");
            }
            return iterator.next();
        }
    }

}
