package org.broadinstitute.hellbender.utils.nio;

import com.google.cloud.storage.contrib.nio.CloudStorageOptions;
import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.seekablestream.ByteArraySeekableStream;
import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.samtools.util.CloseableIterator;
import java.io.IOException;
import java.io.Serializable;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;

/**
 * ReadsIterable gives you all the reads for a given genomic interval.
 *
 * QueryInterval + header --> iterable SAMRecords
 */
public class ReadsIterable implements Iterable<SAMRecord>, Serializable {

    private static final long serialVersionUID = 1L;
    private final String path;
    private final byte[] index;
    private final QueryInterval[] intervals;
    private final boolean removeHeader = true;
    private QueryInterval latestInterval = null;

    class ReadsIterator implements CloseableIterator<SAMRecord> {
        private final static int BUFSIZE = 200 * 1024 * 1024;
        private SamReader bam;
        private SAMRecordIterator query;
        private SAMRecord nextRecord = null;
        private boolean done = false;

        public ReadsIterator() throws IOException {
            Path fpath = IOUtils.getPath(path);
            byte[] indexData = index;
            SeekableStream indexInMemory = new ByteArraySeekableStream(indexData);
            // set high-level retries to deal with servers that might be temporarily overloaded
            // while we're reading a very long file from them.
            SeekableByteChannelPrefetcher chan = new SeekableByteChannelPrefetcher(
                Files.newByteChannel(fpath), BUFSIZE);
            ChannelAsSeekableStream bamOverNIO = new ChannelAsSeekableStream(chan, path);
            bam = SamReaderFactory.makeDefault()
                    .validationStringency(ValidationStringency.LENIENT)
                    .enable(SamReaderFactory.Option.CACHE_FILE_BASED_INDEXES)
                    .open(SamInputResource.of(bamOverNIO).index(indexInMemory));

            query = bam.query(intervals, false);
        }

        /**
         * Returns {@code true} if the iteration has more elements.
         * (In other words, returns {@code true} if {@link #next} would
         * return an element rather than throwing an exception.)
         *
         * @return {@code true} if the iteration has more elements
         */
        @Override
        public boolean hasNext() {
            if (done) return false;

            if (nextRecord!=null) return true;

            nextRecord = fetchRecord();

            boolean ret = (nextRecord != null);
            if (!ret) {
                done = true;
                close();
            }
            return ret;
        }

        /**
         * Returns the next element in the iteration.
         *
         * @return the next element in the iteration
         * @throws NoSuchElementException if the iteration has no more elements
         */
        @Override
        public SAMRecord next() {
            if (!hasNext()) throw new NoSuchElementException();
            SAMRecord ret = nextRecord;
            nextRecord = null;
            return ret;
        }

        private SAMRecord fetchRecord() {
            while (query.hasNext()) {
                SAMRecord sr = query.next();
                int start = sr.getAlignmentStart();
                // only return read if it starts in one of the intervals.
                // First check the previous interval that worked (most likely to work again)
                if (latestInterval != null && start >= latestInterval.start && start <= latestInterval.end) {
                    // read starts in this interval
                    return maybeRemoveHeader(sr);
                }
                for (QueryInterval interval : intervals) {
                    if (start >= interval.start && start <= interval.end) {
                        // read starts in this interval
                        latestInterval = interval;
                        return maybeRemoveHeader(sr);
                    }
                }
                // read doesn't start inside any of our intervals, drop it.
            }
            return null;
        }

        private SAMRecord maybeRemoveHeader(SAMRecord sr) {
            if (removeHeader) {
                sr.setHeader(null);
            }
            return sr;
        }

        @Override
        public void close() {
            if (null==query) return;
            try {
                query.close();
                query = null;
                bam.close();
                bam = null;
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
        }
    }

    // don't modify ins after passing it to us.
    public ReadsIterable(String path, byte[] index, QueryInterval[] ins) {
        this.path = path;
        this.index = index;
        this.intervals = ins;
    }

    @Override
    public Iterator<SAMRecord> iterator() {
        try {
            return new ReadsIterator();
        } catch (IOException x) {
            throw new RuntimeException(x);
        }
    }

}
