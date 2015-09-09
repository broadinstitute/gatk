package htsjdk.samtools.reference;

import com.google.common.io.Files;
import htsjdk.samtools.Defaults;
import htsjdk.samtools.SAMException;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FSDataInputStream;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.FileUtil;
import org.apache.hadoop.fs.Path;

import java.io.File;
import java.io.IOException;
import java.nio.ByteBuffer;

/**
 * <b>Hellbender: Extension of ReferenceSequenceFileFactory to allow reading references from Hadoop paths. See https://github.com/broadinstitute/hellbender/issues/831.</b>
 */
public class HadoopReferenceSequenceFileFactory {

    public static ReferenceSequenceFile getReferenceSequenceFile(Path fastaFile) throws IOException {
        Path fastaIndexFile = getFastaIndexFileName(fastaFile);
        Path fastaDictFile = getFastaDictFileName(fastaFile);

        File tempDir = Files.createTempDir();
        File localFakeFastaFile = new File(tempDir, fastaFile.getName());
        File localFastaIndexFile = new File(tempDir, fastaIndexFile.getName());
        File localFastaDictFile = new File(tempDir, fastaDictFile.getName());

        localFakeFastaFile.createNewFile(); // create empty file to satisfy IndexedFastaSequenceFile; we'll actually read from HDFS

        FileSystem fs = fastaFile.getFileSystem(new Configuration());
        fs.copyToLocalFile(fastaIndexFile, new Path(localFastaIndexFile.toString()));
        if (fs.exists(fastaDictFile)) {
            fs.copyToLocalFile(fastaDictFile, new Path(localFastaDictFile.toString()));
        }

        HadoopFastaSequenceIndex index = new HadoopFastaSequenceIndex(localFastaIndexFile);
        FastaSequenceIndex indexCopy = new FastaSequenceIndex(localFastaIndexFile);
        return new HadoopIndexedFastaSequenceFile(localFakeFastaFile, index, indexCopy, tempDir, fs, fastaFile);
    }

    private static Path getFastaIndexFileName(Path fastaFile) {
        return new Path(fastaFile.toUri().toString() + ".fai");
    }

    private static Path getFastaDictFileName(Path fastaFile) {
        return new Path(fastaFile.toUri().toString().replaceFirst(".fasta$", ".dict"));
    }

    static class HadoopIndexedFastaSequenceFile extends IndexedFastaSequenceFile {

        private HadoopFastaSequenceIndex index;
        private File tempDir;
        private FSDataInputStream in;

        HadoopIndexedFastaSequenceFile(File file, HadoopFastaSequenceIndex index, FastaSequenceIndex indexCopy, File tempDir, FileSystem fs, Path fastaFile) throws IOException {
            super(file, indexCopy);
            this.index = index;
            this.tempDir = tempDir;
            this.in = fs.open(fastaFile);
        }

        /**
         * <b>Hellbender: Copy of IndexedFastaSequenceFile#getSubsequenceAt, updated to use a Hadoop FSDataInputStream
         // rather than a local FileChannel. See comments inline for changes.</b>
         * Gets the subsequence of the contig in the range [start,stop]
         *
         * @param contig Contig whose subsequence to retrieve.
         * @param start  inclusive, 1-based start of region.
         * @param stop   inclusive, 1-based stop of region.
         * @return The partial reference sequence associated with this range.
         */
        @Override
        public ReferenceSequence getSubsequenceAt(String contig, long start, long stop) {
            if (start > stop + 1)
                throw new SAMException(String.format("Malformed query; start point %d lies after end point %d", start, stop));

            HadoopFastaSequenceIndexEntry indexEntry = index.getIndexEntry(contig);

            if (stop > indexEntry.getSize())
                throw new SAMException("Query asks for data past end of contig");

            int length = (int) (stop - start + 1);

            byte[] target = new byte[length];
            ByteBuffer targetBuffer = ByteBuffer.wrap(target);

            final int basesPerLine = indexEntry.getBasesPerLine();
            final int bytesPerLine = indexEntry.getBytesPerLine();
            final int terminatorLength = bytesPerLine - basesPerLine;

            long startOffset = ((start - 1) / basesPerLine) * bytesPerLine + (start - 1) % basesPerLine;

            // Allocate a 128K buffer for reading in sequence data.
            ByteBuffer channelBuffer = ByteBuffer.allocate(Defaults.NON_ZERO_BUFFER_SIZE);
            byte[] channelBytes = new byte[Defaults.NON_ZERO_BUFFER_SIZE];

            while (targetBuffer.position() < length) {
                // If the bufferOffset is currently within the eol characters in the string, push the bufferOffset forward to the next printable character.
                startOffset += Math.max((int) (startOffset % bytesPerLine - basesPerLine + 1), 0);

                try {
                    // Hellbender: read from `in` rather than `channel`
                    //startOffset += channel.read(channelBuffer, indexEntry.getLocation() + startOffset);
                    int bytesRead = in.read(indexEntry.getLocation() + startOffset, channelBytes, 0, channelBytes.length);
                    startOffset += bytesRead;
                    channelBuffer = ByteBuffer.wrap(channelBytes, 0, bytesRead);
                    // End Hellbender changes
                } catch (IOException ex) {
                    throw new SAMException("Unable to load " + contig + "(" + start + ", " + stop + ") from " + file);
                }

                // Reset the buffer for outbound transfers.
                channelBuffer.flip();

                // Calculate the size of the next run of bases based on the contents we've already retrieved.
                final int positionInContig = (int) start - 1 + targetBuffer.position();
                final int nextBaseSpan = Math.min(basesPerLine - positionInContig % basesPerLine, length - targetBuffer.position());
                // Cap the bytes to transfer by limiting the nextBaseSpan to the size of the channel buffer.
                int bytesToTransfer = Math.min(nextBaseSpan, channelBuffer.capacity());

                channelBuffer.limit(channelBuffer.position() + bytesToTransfer);

                while (channelBuffer.hasRemaining()) {
                    targetBuffer.put(channelBuffer);

                    bytesToTransfer = Math.min(basesPerLine, length - targetBuffer.position());
                    channelBuffer.limit(Math.min(channelBuffer.position() + bytesToTransfer + terminatorLength, channelBuffer.capacity()));
                    channelBuffer.position(Math.min(channelBuffer.position() + terminatorLength, channelBuffer.capacity()));
                }

                // Reset the buffer for inbound transfers.
                channelBuffer.flip();
            }

            return new ReferenceSequence(contig, indexEntry.getSequenceIndex(), target);
        }

        @Override
        public void close() throws IOException {
            try {
                super.close();
            } finally {
                in.close();
                FileUtil.fullyDelete(tempDir);
            }
        }
    }
}
