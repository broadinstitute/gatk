package org.broadinstitute.hellbender.tools.walkers.softcliprealignment;

import com.google.common.collect.Lists;
import htsjdk.samtools.*;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.MergingIterator;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.engine.ReadsPathDataSource;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.bwa.BwaReadAligner;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.codehaus.plexus.util.FileUtils;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.util.*;

public class SoftClipRealignmentEngine implements Closeable {

    private final Logger logger = LogManager.getLogger(this.getClass());

    private GATKRead lastRead = null;
    private List<GATKRead> readBuffer;
    private final int bufferSize;
    private final BwaReadAligner aligner;
    private final SAMFileHeader inputHeader;
    private final SAMFileHeader alignmentHeader;
    final File clippedReadBam;
    SAMFileWriter clippedWriter;
    final File unclippedReadBam;
    SAMFileWriter unclippedWriter;

    public static final String REALIGNED_READ_TAG = "RA";

    public SoftClipRealignmentEngine(final String indexImagePath, final SAMFileHeader inputHeader,
                                     final int bufferSize, final int bwaThreads) {
        if (!inputHeader.getSortOrder().equals(SAMFileHeader.SortOrder.coordinate)) {
            throw new UserException.BadInput("Input is not coordinate sorted");
        }
        this.bufferSize = bufferSize;
        readBuffer = new ArrayList<>(bufferSize);
        this.inputHeader = inputHeader;

        // The input must be queryname sorted
        alignmentHeader = inputHeader.clone();
        alignmentHeader.setSortOrder(SAMFileHeader.SortOrder.queryname);
        aligner = new BwaReadAligner(indexImagePath, alignmentHeader, true, true);
        aligner.setNThreads(bwaThreads);

        clippedReadBam = FileUtils.createTempFile("realign_1", FileExtensions.BAM, null);
        unclippedReadBam = FileUtils.createTempFile("realign_2", FileExtensions.BAM, null);
        clippedWriter = new SAMFileWriterFactory().makeWriter(alignmentHeader, false, clippedReadBam, null);
        unclippedWriter = new SAMFileWriterFactory().makeWriter(inputHeader, true, unclippedReadBam, null);
    }

    @Override
    public void close() throws IOException {
        closeFirstPassWriters();
    }

    private void closeFirstPassWriters() {
        if (clippedWriter != null) {
            clippedWriter.close();
            clippedWriter = null;
        }
        if (unclippedWriter != null) {
            unclippedWriter.close();
            unclippedWriter = null;
        }
    }

    public void addRead(final GATKRead read, final Set<String> readNames) {
        if (readNames.contains(read.getName())) {
            if (ReadFilterLibrary.PRIMARY_LINE.test(read)) {
                clippedWriter.addAlignment(read.convertToSAMRecord(clippedWriter.getFileHeader()));
            }
        } else {
            unclippedWriter.addAlignment(read.convertToSAMRecord(unclippedWriter.getFileHeader()));
        }
    }

    public Iterator<GATKRead> run() {
        closeFirstPassWriters();
        final File alignmentBam = FileUtils.createTempFile("realign_3", FileExtensions.BAM, null);
        try (final ReadsPathDataSource clippedSource = new ReadsPathDataSource(clippedReadBam.toPath());
            final SAMFileWriter alignmentWriter = new SAMFileWriterFactory().makeWriter(inputHeader, false, alignmentBam, null)) {
            logger.debug("Created temporary file for realigned reads: " + alignmentBam);
            logger.info("Aligning reads...");
            for (final GATKRead read : clippedSource) {
                thirdPass(read).forEach(r -> alignmentWriter.addAlignment(r.convertToSAMRecord(inputHeader)));
            }
            flushBuffer().forEach(r -> alignmentWriter.addAlignment(r.convertToSAMRecord(inputHeader)));
        }

        // Delete unneeded temp file
        try {
            FileUtils.forceDelete(clippedReadBam);
        } catch (final IOException e) {
            throw new GATKException("Could not clean up temporary file", e);
        }

        try (final ReadsDataSource unclippedSource = new ReadsPathDataSource(unclippedReadBam.toPath());
             final ReadsDataSource alignmentSource = new ReadsPathDataSource(alignmentBam.toPath())) {
            final CloseableIterator<SAMRecord> unclippedRecordsIter = new ClosableReadSourceIterator(unclippedSource);
            final CloseableIterator<SAMRecord> alignedRecordsIter = new ClosableReadSourceIterator(alignmentSource);
            final Comparator<SAMRecord> weakComparator = new Comparator<>() {

                final SAMRecordComparator strongComparator = inputHeader.getSortOrder().getComparatorInstance();
                @Override
                public int compare(SAMRecord o1, SAMRecord o2) {
                    return strongComparator.fileOrderCompare(o1, o2);
                }
            };
            return new MergingIterator<>(weakComparator, List.of(unclippedRecordsIter, alignedRecordsIter))
                    .stream().map(SAMRecordToGATKReadAdapter::new).map(r -> (GATKRead)r).iterator();
        }
    }

    /**
     * Realign the clipped reads
     */
    private List<GATKRead> thirdPass(final GATKRead read) {
        if (lastRead == null) {
            lastRead = read;
        } else if (lastRead.getName().equals(read.getName())) {
            reverseComplementIfNecessary(lastRead);
            reverseComplementIfNecessary(read);
            readBuffer.add(lastRead);
            readBuffer.add(read);
            lastRead = null;
            if (readBuffer.size() == bufferSize || readBuffer.size() + 1 == bufferSize) {
                return flushBuffer();
            }
        } else {
            // Drop an unpaired read and move on to this one
            lastRead = read;
        }
        return Collections.emptyList();
    }

    /**
     * Resets read sequence to forward strand
     */
    private static void reverseComplementIfNecessary(final GATKRead read) {
        if (read.isReverseStrand()) {
            read.reverseComplement();
        }
    }

    public List<GATKRead> flushBuffer() {
        final List<GATKRead> alignments = Lists.newArrayList(aligner.apply(readBuffer.iterator()));
        alignments.forEach(SoftClipRealignmentEngine::setRealignedTag);
        readBuffer.clear();
        return alignments;
    }

    /**
     * Sets tag {@link #REALIGNED_READ_TAG} to mark it as realigned
     */
    static void setRealignedTag(final GATKRead read) {
        read.setAttribute(REALIGNED_READ_TAG, 1);
    }

    /**
     * Helper class to map a {@link ReadsDataSource} iterator to a {@link CloseableIterator} and each input read to a
     * {@link SAMRecord} as well.
     */
    static final class ClosableReadSourceIterator implements CloseableIterator<SAMRecord> {

        final ReadsDataSource source;
        final Iterator<GATKRead> iterator;

        public ClosableReadSourceIterator(final ReadsDataSource source) {
            this.source = source;
            this.iterator = source.iterator();
        }

        @Override
        public void close() {
            source.close();
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
            return iterator.next().convertToSAMRecord(source.getHeader());
        }
    }

}
