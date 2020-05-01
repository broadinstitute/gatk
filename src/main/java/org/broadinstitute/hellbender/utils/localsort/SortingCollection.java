package org.broadinstitute.hellbender.utils.localsort;

import htsjdk.samtools.Defaults;
import htsjdk.samtools.util.*;
import org.apache.avro.generic.GenericRecord;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.tools.variantdb.SchemaUtils;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.lang.reflect.Array;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.TreeSet;

/**
 * Collection to which many records can be added.  After all records are added, the collection can be
 * iterated, and the records will be returned in order defined by the comparator.  Records may be spilled
 * to a temporary directory if there are more records added than will fit in memory.  As a result of this,
 * the objects returned may not be identical to the objects added to the collection, but they should be
 * equal as determined by the codec used to write them to disk and read them back.
 * <p>
 * When iterating over the collection, the number of file handles required is numRecordsInCollection/maxRecordsInRam.
 * If this becomes a limiting factor, a file handle cache could be added.
 * <p>
 * If Snappy DLL is available and snappy.disable system property is not set to true, then Snappy is used
 * to compress temporary files.
 */
public class SortingCollection<T> implements Iterable<T> {
    private static final Logger log = LogManager.getLogger(SortingCollection.class);

    /**
     * Client must implement this class, which defines the way in which records are written to and
     * read from file.
     */
    public interface Codec<T> extends Cloneable {
        /**
         * Where to write encoded output
         */
        void setOutputStream(OutputStream os);

        void flushOutput();

        /**
         * Where to read encoded input from
         */
        void setInputStream(InputStream is);

        /**
         * Write object to output stream
         *
         * @param val what to write
         */
        void encode(T val);

        /**
         * Read the next record from the input stream and convert into a java object.
         *
         * @return null if no more records.  Should throw exception if EOF is encountered in the middle of
         * a record.
         */
        T decode();

        /**
         * Must return a cloned copy of the codec that can be used independently of
         * the original instance.  This is required so that multiple codecs can exist simultaneously
         * that each is reading a separate file.
         */
        Codec<T> clone();
    }

    /**
     * Directories where files of sorted records go.
     */
    private final Path[] tmpDirs;

    /**
     * Used to write records to file, and used as a prototype to create codecs for reading.
     */
    private final SortingCollection.Codec<T> codec;

    /**
     * For sorting, both when spilling records to file, and merge sorting.
     */
    private final Comparator<T> comparator;
    private final int maxRecordsInRam;
    private int numRecordsInRam = 0;
    private T[] ramRecords;
    private boolean iterationStarted = false;
    private boolean doneAdding = false;

    /**
     * Set to true when all temp files have been cleaned up
     */
    private boolean cleanedUp = false;

    /**
     * List of files in tmpDir containing sorted records
     */
    private final List<Path> files = new ArrayList<>();

    private boolean destructiveIteration = true;

    private final TempStreamFactory tempStreamFactory = new TempStreamFactory();

    private final boolean printRecordSizeSampling;

    /**
     * Prepare to accumulate records to be sorted
     *
     * @param componentType   Class of the record to be sorted.  Necessary because of Java generic lameness.
     * @param codec           For writing records to file and reading them back into RAM
     * @param comparator      Defines output sort order
     * @param maxRecordsInRam how many records to accumulate before spilling to disk
     * @param printRecordSizeSampling If true record size will be sampled and output at DEBUG log level
     * @param tmpDir          Where to write files of records that will not fit in RAM
     */
    private SortingCollection(final Class<T> componentType, final SortingCollection.Codec<T> codec,
                              final Comparator<T> comparator, final int maxRecordsInRam,
                              final boolean printRecordSizeSampling, final Path... tmpDir) {
        if (maxRecordsInRam <= 0) {
            throw new IllegalArgumentException("maxRecordsInRam must be > 0");
        }

        if (tmpDir == null || tmpDir.length == 0) {
            throw new IllegalArgumentException("At least one temp directory must be provided.");
        }

        this.tmpDirs = tmpDir;
        this.codec = codec;
        this.comparator = comparator;
        this.maxRecordsInRam = maxRecordsInRam;
        @SuppressWarnings("unchecked")
        T[] ramRecords = (T[]) Array.newInstance(componentType, maxRecordsInRam);
        this.ramRecords = ramRecords;
        this.printRecordSizeSampling = printRecordSizeSampling;
    }

    public void add(final T rec) {
        if (doneAdding) {
            throw new IllegalStateException("Cannot add after calling doneAdding()");
        }
        if (iterationStarted) {
            throw new IllegalStateException("Cannot add after calling iterator()");
        }
        if (numRecordsInRam == maxRecordsInRam) {

            long startMem = 0;
            if (printRecordSizeSampling) {
                // Garbage collect and get free memory
                Runtime.getRuntime().gc();
                startMem = Runtime.getRuntime().freeMemory();
            }

            spillToDisk();

            if (printRecordSizeSampling) {
                //Garbage collect again and get free memory
                Runtime.getRuntime().gc();
                long endMem = Runtime.getRuntime().freeMemory();

                long usedBytes = endMem - startMem;
                log.debug(String.format("%d records in ram required approximately %s memory or %s per record. ", maxRecordsInRam,
                        StringUtil.humanReadableByteCount(usedBytes),
                        StringUtil.humanReadableByteCount(usedBytes / maxRecordsInRam)));

            }
        }
        ramRecords[numRecordsInRam++] = rec;
    }

    /**
     * This method can be called after caller is done adding to collection, in order to possibly free
     * up memory.  If iterator() is called immediately after caller is done adding, this is not necessary,
     * because iterator() triggers the same freeing.
     */
    public void doneAdding() {
        if (this.cleanedUp) {
            throw new IllegalStateException("Cannot call doneAdding() after cleanup() was called.");
        }
        if (doneAdding) {
            return;
        }

        doneAdding = true;

        if (this.files.isEmpty()) {
            return;
        }

        if (this.numRecordsInRam > 0) {
            spillToDisk();
        }

        // Facilitate GC
        this.ramRecords = null;
    }

    /**
     * @return True if this collection is allowed to discard data during iteration in order to reduce memory
     * footprint, precluding a second iteration over the collection.
     */
    public boolean isDestructiveIteration() {
        return destructiveIteration;
    }

    /**
     * Tell this collection that it is allowed to discard data during iteration in order to reduce memory footprint,
     * precluding a second iteration.  This is true by default.
     */
    public void setDestructiveIteration(boolean destructiveIteration) {
        this.destructiveIteration = destructiveIteration;
    }

    /**
     * Sort the records in memory, write them to a file, and clear the buffer of records in memory.
     */
    public void spillToDisk() {
        try {
            Arrays.parallelSort(this.ramRecords, 0, this.numRecordsInRam, this.comparator);

            final Path f = newTempFile();
            try (OutputStream os
                         = tempStreamFactory.wrapTempOutputStream(Files.newOutputStream(f), Defaults.BUFFER_SIZE)) {
                this.codec.setOutputStream(os);
                for (int i = 0; i < this.numRecordsInRam; ++i) {
                    this.codec.encode(ramRecords[i]);
                    // Facilitate GC
                    this.ramRecords[i] = null;
                }
                this.codec.flushOutput();
            } catch (RuntimeIOException ex) {
                throw new RuntimeIOException("Problem writing temporary file " + f.toUri() +
                        ".  Try setting TMP_DIR to a file system with lots of space.", ex);
            }

            this.numRecordsInRam = 0;
            this.files.add(f);
        } catch (IOException e) {
            throw new RuntimeIOException(e);
        }
    }


    /**
     * Creates a new tmp file on one of the available temp filesystems, registers it for deletion
     * on JVM exit and then returns it.
     */
    private Path newTempFile() throws IOException {
        /* The minimum amount of space free on a temp filesystem to write a file there. */
        return IOUtil.newTempPath("sortingcollection.", ".tmp", this.tmpDirs, IOUtil.FIVE_GBS);
    }

    public void printTempFileStats() {
        log.info(String.format("Wrote %d sorted temp files\n", this.files.size()));

        if ( this.files.isEmpty() ) {
            return;
        }

        long maxSize = Long.MIN_VALUE;
        long minSize = Long.MAX_VALUE;
        long totalSize = 0;
        for ( final Path tempFile : this.files ) {
            try {
                long size = Files.size(tempFile);
                totalSize += size;
                if ( size > maxSize ) {
                    maxSize = size;
                }
                if ( size < minSize ) {
                    minSize = size;
                }
            } catch ( IOException e ) {
                log.warn(String.format("Error getting size of temp file %s\n", tempFile));
            }
        }

        log.info(String.format("Temp files total size: %d bytes\n", totalSize));
        log.info(String.format("Temp files max size: %d bytes\n", maxSize));
        log.info(String.format("Temp files min size: %d bytes\n", minSize));
        log.info(String.format("Temp files average size: %.2f bytes\n", (totalSize / (double)this.files.size())));
    }

    /**
     * Prepare to iterate through the records in order.  This method may be called more than once,
     * but add() may not be called after this method has been called.
     */
    @Override
    public CloseableIterator<T> iterator() {
        if (this.cleanedUp) {
            throw new IllegalStateException("Cannot call iterator() after cleanup() was called.");
        }
        doneAdding();

        this.iterationStarted = true;
        if (this.files.isEmpty()) {
            return new InMemoryIterator();
        } else {
            return new MergingIterator();
        }
    }

    /**
     * Delete any temporary files.  After this method is called, iterator() may not be called.
     */
    public void cleanup() {
        this.iterationStarted = true;
        this.cleanedUp = true;

        IOUtil.deletePaths(this.files);
    }

    /**
     * Syntactic sugar around the ctor, to save some typing of type parameters
     *
     * @param componentType   Class of the record to be sorted.  Necessary because of Java generic lameness.
     * @param codec           For writing records to file and reading them back into RAM
     * @param comparator      Defines output sort order
     * @param maxRecordsInRAM how many records to accumulate in memory before spilling to disk
     * @param tmpDir          Where to write files of records that will not fit in RAM
     * @deprecated since 2017-09. Use {@link #newInstance(Class, Codec, Comparator, int, Path...)} instead
     */
    @Deprecated
    public static <T> SortingCollection<T> newInstance(final Class<T> componentType,
                                                       final SortingCollection.Codec<T> codec,
                                                       final Comparator<T> comparator,
                                                       final int maxRecordsInRAM,
                                                       final File... tmpDir) {
        return new SortingCollection<>(componentType, codec, comparator, maxRecordsInRAM, false, Arrays.stream(tmpDir).map(File::toPath).toArray(Path[]::new));

    }

    /**
     * Syntactic sugar around the ctor, to save some typing of type parameters
     *
     * @param componentType   Class of the record to be sorted.  Necessary because of Java generic lameness.
     * @param codec           For writing records to file and reading them back into RAM
     * @param comparator      Defines output sort order
     * @param maxRecordsInRAM how many records to accumulate in memory before spilling to disk
     * @param tmpDirs         Where to write files of records that will not fit in RAM
     * @deprecated since 2017-09. Use {@link #newInstanceFromPaths(Class, Codec, Comparator, int, Collection)} instead
     */
    @Deprecated
    public static <T> SortingCollection<T> newInstance(final Class<T> componentType,
                                                       final SortingCollection.Codec<T> codec,
                                                       final Comparator<T> comparator,
                                                       final int maxRecordsInRAM,
                                                       final Collection<File> tmpDirs) {
        return new SortingCollection<>(componentType,
                codec,
                comparator,
                maxRecordsInRAM,
                false,
                tmpDirs.stream().map(File::toPath).toArray(Path[]::new));

    }

    /**
     * Syntactic sugar around the ctor, to save some typing of type parameters.  Writes files to java.io.tmpdir
     *
     * @param componentType    Class of the record to be sorted.  Necessary because of Java generic lameness.
     * @param codec            For writing records to file and reading them back into RAM
     * @param comparator       Defines output sort order
     * @param maxRecordsInRAM  how many records to accumulate in memory before spilling to disk
     * @param printRecordSizeSampling If true record size will be sampled and output at DEBUG log level
     */
    public static <T> SortingCollection<T> newInstance(final Class<T> componentType,
                                                       final SortingCollection.Codec<T> codec,
                                                       final Comparator<T> comparator,
                                                       final int maxRecordsInRAM,
                                                       final boolean printRecordSizeSampling) {
        final Path tmpDir = Paths.get(System.getProperty("java.io.tmpdir"));
        return new SortingCollection<>(componentType, codec, comparator, maxRecordsInRAM, printRecordSizeSampling, tmpDir);
    }

    /**
     * Syntactic sugar around the ctor, to save some typing of type parameters
     *
     * @param componentType    Class of the record to be sorted.  Necessary because of Java generic lameness.
     * @param codec            For writing records to file and reading them back into RAM
     * @param comparator       Defines output sort order
     * @param maxRecordsInRAM  how many records to accumulate in memory before spilling to disk
     * @param printRecordSizeSampling If true record size will be sampled and output at DEBUG log level
     * @param tmpDir           Where to write files of records that will not fit in RAM
     */
    public static <T> SortingCollection<T> newInstance(final Class<T> componentType,
                                                       final SortingCollection.Codec<T> codec,
                                                       final Comparator<T> comparator,
                                                       final int maxRecordsInRAM,
                                                       final boolean printRecordSizeSampling,
                                                       final Path... tmpDir) {
        return new SortingCollection<>(componentType, codec, comparator, maxRecordsInRAM, printRecordSizeSampling, tmpDir);
    }

    /**
     * Syntactic sugar around the ctor, to save some typing of type parameters.  Writes files to java.io.tmpdir
     *
     * @param componentType   Class of the record to be sorted.  Necessary because of Java generic lameness.
     * @param codec           For writing records to file and reading them back into RAM
     * @param comparator      Defines output sort order
     * @param maxRecordsInRAM how many records to accumulate in memory before spilling to disk
     */
    public static <T> SortingCollection<T> newInstance(final Class<T> componentType,
                                                       final SortingCollection.Codec<T> codec,
                                                       final Comparator<T> comparator,
                                                       final int maxRecordsInRAM) {
        final Path tmpDir = Paths.get(System.getProperty("java.io.tmpdir"));
        return new SortingCollection<>(componentType, codec, comparator, maxRecordsInRAM, false, tmpDir);
    }

    /**
     * Syntactic sugar around the ctor, to save some typing of type parameters
     *
     * @param componentType   Class of the record to be sorted.  Necessary because of Java generic lameness.
     * @param codec           For writing records to file and reading them back into RAM
     * @param comparator      Defines output sort order
     * @param maxRecordsInRAM how many records to accumulate in memory before spilling to disk
     * @param tmpDir          Where to write files of records that will not fit in RAM
     */
    public static <T> SortingCollection<T> newInstance(final Class<T> componentType,
                                                       final SortingCollection.Codec<T> codec,
                                                       final Comparator<T> comparator,
                                                       final int maxRecordsInRAM,
                                                       final Path... tmpDir) {
        return new SortingCollection<>(componentType, codec, comparator, maxRecordsInRAM, false, tmpDir);
    }

    /**
     * Syntactic sugar around the ctor, to save some typing of type parameters
     *
     * @param componentType   Class of the record to be sorted.  Necessary because of Java generic lameness.
     * @param codec           For writing records to file and reading them back into RAM
     * @param comparator      Defines output sort order
     * @param maxRecordsInRAM how many records to accumulate in memory before spilling to disk
     * @param tmpDirs         Where to write files of records that will not fit in RAM
     */
    public static <T> SortingCollection<T> newInstanceFromPaths(final Class<T> componentType,
                                                                final SortingCollection.Codec<T> codec,
                                                                final Comparator<T> comparator,
                                                                final int maxRecordsInRAM,
                                                                final Collection<Path> tmpDirs) {
        return new SortingCollection<>(componentType,
                codec,
                comparator,
                maxRecordsInRAM,
                false,
                tmpDirs.toArray(new Path[tmpDirs.size()]));
    }

    /**
     * For iteration when number of records added is less than the threshold for spilling to disk.
     */
    class InMemoryIterator implements CloseableIterator<T> {
        private int iterationIndex = 0;

        InMemoryIterator() {
            Arrays.parallelSort(SortingCollection.this.ramRecords,
                    0,
                    SortingCollection.this.numRecordsInRam,
                    SortingCollection.this.comparator);
        }

        @Override
        public void close() {
            // nothing to do
        }

        @Override
        public boolean hasNext() {
            return this.iterationIndex < SortingCollection.this.numRecordsInRam;
        }

        @Override
        public T next() {
            if (!hasNext()) {
                throw new NoSuchElementException();
            }
            T ret = SortingCollection.this.ramRecords[iterationIndex];
            if (destructiveIteration) SortingCollection.this.ramRecords[iterationIndex] = null;
            ++iterationIndex;
            return ret;
        }

        @Override
        public void remove() {
            throw new UnsupportedOperationException();
        }
    }

    /**
     * For iteration when spilling to disk has occurred.
     * Each file is has records in sort order within the file.
     * This iterator automatically closes when it iterates to the end, but if not iterating
     * to the end it is a good idea to call close().
     * <p>
     * Algorithm: MergingIterator maintains a PriorityQueue of PeekFileRecordIterators.
     * Each PeekFileRecordIterator iterates through a file in which the records are sorted.
     * The comparator for PeekFileRecordIterator used by the PriorityQueue peeks at the next record from
     * the file, so the first element in the PriorityQueue is the file that has the next record to be emitted.
     * In order to get the next record, the first PeekFileRecordIterator in the PriorityQueue is popped,
     * the record is obtained from that iterator, and then if that iterator is not empty, it is pushed back into
     * the PriorityQueue.  Because it now has a different record as its next element, it may go into another
     * location in the PriorityQueue
     */
    class MergingIterator implements CloseableIterator<T> {
        private final TreeSet<PeekFileRecordIterator> queue;

        MergingIterator() {
            this.queue = new TreeSet<>(new PeekFileRecordIteratorComparator());
            int n = 0;
            log.debug(String.format("Creating merging iterator from %d files", files.size()));
            int suggestedBufferSize = checkMemoryAndAdjustBuffer(files.size());
            for (final Path f : files) {
                final FileRecordIterator it = new FileRecordIterator(f, suggestedBufferSize);
                if (it.hasNext()) {
                    this.queue.add(new PeekFileRecordIterator(it, n++));
                } else {
                    it.close();
                }
            }
        }

        // Since we need to open and buffer all temp files in the sorting collection at once it is important
        // to have enough memory left to do this. This method checks to make sure that, given the number of files and
        // the size of the buffer, we can reasonably open all files. If we can't it will return a buffer size that
        // is appropriate given the number of temp files and the amount of memory left on the heap. If there isn't
        // enough memory for buffering it will return zero and all reading will be unbuffered.
        private int checkMemoryAndAdjustBuffer(int numFiles) {
            int bufferSize = Defaults.BUFFER_SIZE;
            // garbage collect so that our calculation is accurate.
            Runtime.getRuntime().gc();

            // There is ~20k in overhead per file.
            final long freeMemory = Runtime.getRuntime().freeMemory() - (numFiles * 20 * 1024);
            // use the floor value from the divide
            final int memoryPerFile = (int) (freeMemory / numFiles);

            if (memoryPerFile < 0) {
                log.warn("There is not enough memory per file for buffering. Reading will be unbuffered.");
                bufferSize = 0;
            } else if (bufferSize > memoryPerFile) {
                log.warn(String.format("Default io buffer size of %s is larger than available memory per file of %s.",
                        StringUtil.humanReadableByteCount(bufferSize),
                        StringUtil.humanReadableByteCount(memoryPerFile)));
                bufferSize = memoryPerFile;
            }
            return bufferSize;
        }

        @Override
        public boolean hasNext() {
            return !this.queue.isEmpty();
        }

        @Override
        public T next() {
            if (!hasNext()) {
                throw new NoSuchElementException();
            }

            final PeekFileRecordIterator fileIterator = queue.pollFirst();
            final T ret = fileIterator.next();
            if (fileIterator.hasNext()) {
                this.queue.add(fileIterator);
            } else {
                ((CloseableIterator<T>) fileIterator.getUnderlyingIterator()).close();
            }

            return ret;
        }

        @Override
        public void remove() {
            throw new UnsupportedOperationException();
        }

        @Override
        public void close() {
            while (!this.queue.isEmpty()) {
                final PeekFileRecordIterator it = this.queue.pollFirst();
                ((CloseableIterator<T>) it.getUnderlyingIterator()).close();
            }
        }
    }

    /**
     * Read a file of records in format defined by the codec
     */
    class FileRecordIterator implements CloseableIterator<T> {
        private final Path file;
        private final InputStream is;
        private final Codec<T> codec;
        private T currentRecord = null;

        FileRecordIterator(final Path file, final int bufferSize) {
            this.file = file;
            try {
                this.is = Files.newInputStream(file);
                this.codec = SortingCollection.this.codec.clone();
                this.codec.setInputStream(tempStreamFactory.wrapTempInputStream(this.is, bufferSize));
                advance();
            } catch (IOException e) {
                throw new RuntimeIOException(e);
            }
        }

        @Override
        public boolean hasNext() {
            return this.currentRecord != null;
        }

        @Override
        public T next() {
            if (!hasNext()) {
                throw new NoSuchElementException();
            }
            final T ret = this.currentRecord;
            advance();
            return ret;
        }

        @Override
        public void remove() {
            throw new UnsupportedOperationException();
        }

        private void advance() {
            this.currentRecord = this.codec.decode();
        }

        @Override
        public void close() {
            CloserUtil.close(this.is);
        }
    }


    /**
     * Just a typedef
     */
    class PeekFileRecordIterator extends PeekIterator<T> {
        final int n; // A serial number used for tie-breaking in the sort

        PeekFileRecordIterator(final Iterator<T> underlyingIterator, final int n) {
            super(underlyingIterator);
            this.n = n;
        }
    }

    class PeekFileRecordIteratorComparator implements Comparator<PeekFileRecordIterator> {
        @Override
        public int compare(final PeekFileRecordIterator lhs, final PeekFileRecordIterator rhs) {
            final int result = comparator.compare(lhs.peek(), rhs.peek());
            if (result == 0) return lhs.n - rhs.n;
            else return result;
        }
    }
}
