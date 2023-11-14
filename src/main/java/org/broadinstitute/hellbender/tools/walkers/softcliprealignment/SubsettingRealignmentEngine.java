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
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.codehaus.plexus.util.FileUtils;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.function.Predicate;

public final class SubsettingRealignmentEngine implements AutoCloseable {

    private GATKRead lastRead = null;
    private List<GATKRead> pairedBuffer;
    private List<GATKRead> unpairedBuffer;
    private final int bufferSize;

    private final SAMFileHeader inputHeader;
    private final SAMFileHeader alignmentHeader;

    private final BwaReadAligner pairedAligner;
    private final BwaReadAligner unpairedAligner;

    private File selectedReadsBam;
    private File nonselectedReadsBam;
    private SAMFileWriter selectedReadsWriter;
    private SAMFileWriter nonselectedReadsWriter;

    private CloseableIterator<SAMRecord> unselectedRecordsIter;
    private CloseableIterator<SAMRecord> alignedRecordsIter;

    private long selectedReadsCount = 0;
    private long nonselectedReadsCount = 0;
    private long pairedAlignmentReadsCount = 0;
    private long unpairedAlignmentReadsCount = 0;

    public static final String REALIGNED_READ_TAG = "RA";
    public static final int REALIGNED_READ_TAG_VALUE = 1;

    public SubsettingRealignmentEngine(final String indexImagePath, final SAMFileHeader inputHeader,
                                       final int bufferSize, final int bwaThreads) {
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

        pairedAligner = new BwaReadAligner(indexImagePath, alignmentHeader, true, true);
        pairedAligner.setNThreads(bwaThreads);

        unpairedAligner = new BwaReadAligner(indexImagePath, alignmentHeader, false, true);
        unpairedAligner.setNThreads(bwaThreads);

        selectedReadsBam = FileUtils.createTempFile("realign_1", FileExtensions.BAM, null);
        nonselectedReadsBam = FileUtils.createTempFile("realign_2", FileExtensions.BAM, null);
        selectedReadsWriter = new SAMFileWriterFactory().makeWriter(alignmentHeader, false, selectedReadsBam, null);
        nonselectedReadsWriter = new SAMFileWriterFactory().makeWriter(inputHeader, true, nonselectedReadsBam, null);
    }

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

    public void addRead(final GATKRead read, final Predicate<GATKRead> predicate) {
        Utils.nonNull(selectedReadsWriter, "This instance has been closed");
        Utils.nonNull(nonselectedReadsWriter, "This instance has been closed");
        if (predicate.test(read)) {
            if (ReadFilterLibrary.PRIMARY_LINE.test(read)) {
                selectedReadsWriter.addAlignment(read.convertToSAMRecord(selectedReadsWriter.getFileHeader()));
                selectedReadsCount++;
            }
        } else {
            nonselectedReadsWriter.addAlignment(read.convertToSAMRecord(nonselectedReadsWriter.getFileHeader()));
            nonselectedReadsCount++;
        }
    }

    public Iterator<GATKRead> alignAndMerge() {
        Utils.nonNull(selectedReadsBam, "This instance has been run already");
        Utils.nonNull(nonselectedReadsBam, "This instance has been closed");
        // Assume we're done with adding reads, so we can close the first set of temp files
        closeFirstPassWriters();
        
        // Align selected reads
        final File alignmentBam = FileUtils.createTempFile("realign_3", FileExtensions.BAM, null);
        try (final ReadsDataSource selectedSource = new ReadsPathDataSource(selectedReadsBam.toPath());
             final SAMFileWriter alignmentWriter = new SAMFileWriterFactory().makeWriter(inputHeader, false, alignmentBam, null)) {
            for (final GATKRead read : selectedSource) {
                addReadToBuffer(read).forEach(r -> alignmentWriter.addAlignment(r.convertToSAMRecord(inputHeader)));
            }
            // Flush leftover reads
            writeAlignments(flushPairedBuffer(), alignmentWriter);
            writeAlignments(flushUnpairedBuffer(), alignmentWriter);
            writeAlignments(flushLastRead(), alignmentWriter);
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
        final Comparator<SAMRecord> weakComparator = new Comparator<>() {
            final SAMRecordComparator strongComparator = inputHeader.getSortOrder().getComparatorInstance();
            @Override
            public int compare(SAMRecord o1, SAMRecord o2) {
                return strongComparator.fileOrderCompare(o1, o2);
            }
        };
        return new MergingIterator<>(weakComparator, List.of(unselectedRecordsIter, alignedRecordsIter))
                .stream().map(SAMRecordToGATKReadAdapter::new).map(r -> (GATKRead)r).iterator();
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

    private List<GATKRead> flushLastRead() {
        if (lastRead == null) {
            return Collections.emptyList();
        } else {
            final ArrayList<GATKRead> buffer = new ArrayList<>(1);
            buffer.add(lastRead);
            unpairedAlignmentReadsCount++;
            return flushBuffer(buffer, unpairedAligner);
        }
    }

    private List<GATKRead> flushBuffer(final List<GATKRead> mutableBuffer, final BwaReadAligner aligner) {
        final List<GATKRead> alignments = Lists.newArrayList(aligner.apply(mutableBuffer));
        alignments.forEach(SubsettingRealignmentEngine::setRealignedTag);
        mutableBuffer.clear();
        return alignments;
    }

    /**
     * Sets tag {@link #REALIGNED_READ_TAG} to mark it as realigned
     */
    static void setRealignedTag(final GATKRead read) {
        read.setAttribute(REALIGNED_READ_TAG, REALIGNED_READ_TAG_VALUE);
    }

    public static boolean checkRealignedTag(final GATKRead read) {
        return read.hasAttribute(REALIGNED_READ_TAG) && read.getAttributeAsInteger(REALIGNED_READ_TAG) == REALIGNED_READ_TAG_VALUE;
    }

    public long getSelectedReadsCount() {
        return selectedReadsCount;
    }

    public long getNonselectedReadsCount() {
        return nonselectedReadsCount;
    }

    public long getPairedAlignmentReadsCount() {
        return pairedAlignmentReadsCount;
    }

    public long getUnpairedAlignmentReadsCount() {
        return unpairedAlignmentReadsCount;
    }

    /**
     * Helper class to map a {@link ReadsDataSource} iterator to a {@link CloseableIterator} and each input read to a
     * {@link SAMRecord} as well.
     */
    static final class ClosableReadSourceIterator implements CloseableIterator<SAMRecord> {

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
        public SAMRecord next() {
            if (source == null) {
                throw new NoSuchElementException("Source has been closed");
            }
            return iterator.next().convertToSAMRecord(source.getHeader());
        }
    }

}
