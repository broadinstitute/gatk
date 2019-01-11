package org.broadinstitute.hellbender.tools.diagnostics;

import htsjdk.samtools.util.BlockCompressedStreamConstants;
import htsjdk.samtools.util.IOUtil;
import org.apache.commons.lang3.StringUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import picard.cmdline.programgroups.OtherProgramGroup;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Path;

/**
 * A diagnostic tool that prints information about the compressed blocks in a BGZF format file,
 * such as a .vcf.gz file.
 * <p>
 * The output looks like this:
 * </p>
 * <pre>
 *     Block at file offset 0
 *         - compressed size: 12932
 *         - uncompressed size: 65280
 *
 *     Block at file offset 12932
 *         - compressed size: 9978
 *         - uncompressed size: 65280
 *     ...
 *     etc.
 * </pre>
 * <p>
 * The output can be redirected to a file using the -O option.
 * </p>
 */
@ExperimentalFeature
@CommandLineProgramProperties(
        summary = "Print information about the compressed blocks in a BGZF format file",
        oneLineSummary = "Print information about the compressed blocks in a BGZF format file",
        programGroup = OtherProgramGroup.class
)
public class PrintBGZFBlockInformation extends CommandLineProgram {

    @Argument(fullName = "bgzf-file", doc = "The BGZF-format file for which to print block information", optional = false)
    private String bgzfPathString;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "File to which to write block information (if not specified, prints to standard output)", optional = true)
    private String output;

    private Path bgzfPath;

    private long streamOffset = 0l;

    private PrintStream outStream;

    @Override
    protected void onStartup() {
        super.onStartup();

        bgzfPath = IOUtils.getPath(bgzfPathString);

        if ( ! Files.exists(bgzfPath) ) {
            throw new UserException.CouldNotReadInputFile("File " + bgzfPathString + " does not exist");
        }

        if ( ! IOUtil.hasBlockCompressedExtension(bgzfPathString) ) {
            throw new UserException.CouldNotReadInputFile("File " + bgzfPathString + " does not end in a recognized BGZF file extension (" +
                    StringUtils.join(IOUtil.BLOCK_COMPRESSED_EXTENSIONS, ",") + ")");
        }

        try {
            // Check that the file is in BGZF format. This catches the "regular GZIP" case as well:
            if ( ! IOUtil.isBlockCompressed(bgzfPath) ) {
                throw new UserException.CouldNotReadInputFile(bgzfPath, "File is not a valid BGZF file. Could possibly be a regular GZIP file?");
            }
        }
        catch ( IOException e ) {
            throw new UserException.CouldNotReadInputFile(bgzfPath, "Unable to determine whether file is a valid BGZF file", e);
        }

        if ( output != null ) {
            try {
                outStream = new PrintStream(output);
            } catch (FileNotFoundException e) {
                throw new UserException.CouldNotCreateOutputFile(output, "Unable to open output file", e);
            }
        } else {
            outStream = System.out;
        }
    }

    @Override
    protected Object doWork() {
        BGZFBlockMetadata previousBlockInfo = null;
        boolean sawNonFinalTerminatorBlock = false;

        try ( InputStream bgzfInputStream = Files.newInputStream(bgzfPath) ) {
            outStream.printf("BGZF block information for file: %s\n\n", bgzfPath.getFileName());

            BGZFBlockMetadata blockInfo;

            while ( (blockInfo = processNextBlock(bgzfInputStream, bgzfPathString)) != null ) {

                // If we saw a 0-byte terminator block that was not the final block in the file,
                // emit a warning
                if ( previousBlockInfo != null && previousBlockInfo.uncompressedSize == 0 ) {
                    sawNonFinalTerminatorBlock = true;

                    outStream.println("*************************************************************************");
                    outStream.println("WARNING: Premature BGZF 0-byte terminator block found before final block!");
                    outStream.println("*************************************************************************");
                    outStream.println();
                }

                outStream.printf("Block at file offset %d\n", blockInfo.blockOffset);
                outStream.printf("\t- compressed size: %d\n", blockInfo.compressedSize);
                outStream.printf("\t- uncompressed size: %d\n", blockInfo.uncompressedSize);
                outStream.println();

                previousBlockInfo = blockInfo;
            }
        } catch ( IOException e ) {
            throw new UserException.CouldNotReadInputFile("Error while parsing BGZF file.", e);
        }

        // Check whether the last block in the file was a 0-byte BGZF terminator block
        if ( previousBlockInfo == null || previousBlockInfo.uncompressedSize != 0 ) {
            outStream.println("********************************************************");
            outStream.println("WARNING: Final BGZF 0-byte terminator block was MISSING!");
            outStream.println("********************************************************");
            outStream.println();
        } else {
            outStream.println("*****************************************");
            outStream.println("Final BGZF 0-byte terminator block FOUND!");
            outStream.println("*****************************************");
            outStream.println();
        }

        // Emit a warning at the end if we encountered any terminator blocks before the final block:
        if ( sawNonFinalTerminatorBlock ) {
            outStream.println("*************************************************************************");
            outStream.println("WARNING: Premature BGZF 0-byte terminator block found before final block!");
            outStream.println("*************************************************************************");
            outStream.println();
        }

        return 0;
    }

    @Override
    protected void onShutdown() {
        if ( outStream != null && outStream != System.out ) {
            outStream.close();
        }
    }

    // Code adapted from HTSJDK's BlockCompressedInputStream class
    private BGZFBlockMetadata processNextBlock(InputStream stream, String streamSource) throws IOException {
        final byte[] buffer = new byte[BlockCompressedStreamConstants.MAX_COMPRESSED_BLOCK_SIZE];
        long blockAddress = streamOffset;

        final int headerByteCount = readBytes(stream, buffer, 0, BlockCompressedStreamConstants.BLOCK_HEADER_LENGTH);

        // Return null when we hit EOF
        if ( headerByteCount <= 0 ) {
            return null;
        }
        if (headerByteCount != BlockCompressedStreamConstants.BLOCK_HEADER_LENGTH) {
            throw new IOException("Incorrect header size for file: " + streamSource);
        }
        streamOffset += headerByteCount;

        final int blockLength = unpackInt16(buffer, BlockCompressedStreamConstants.BLOCK_LENGTH_OFFSET) + 1;

        if (blockLength < BlockCompressedStreamConstants.BLOCK_HEADER_LENGTH || blockLength > buffer.length) {
            throw new IOException("Unexpected compressed block length: " + blockLength + " for " + streamSource);
        }

        final int remaining = blockLength - BlockCompressedStreamConstants.BLOCK_HEADER_LENGTH;
        final int dataByteCount = readBytes(stream, buffer, BlockCompressedStreamConstants.BLOCK_HEADER_LENGTH,
                remaining);

        if (dataByteCount != remaining) {
            throw new IOException("Premature end of file: " + streamSource);
        }
        streamOffset += dataByteCount;

        final int uncompressedLength = unpackInt32(buffer, blockLength - 4);

        if (uncompressedLength < 0) {
            throw new IOException(streamSource + " has invalid uncompressed length: " + uncompressedLength);
        }

        return new BGZFBlockMetadata(blockAddress, blockLength, uncompressedLength);
    }

    private static int unpackInt16(final byte[] buffer, final int offset) {
        return ((buffer[offset] & 0xFF) |
                ((buffer[offset+1] & 0xFF) << 8));
    }

    private static int unpackInt32(final byte[] buffer, final int offset) {
        return ((buffer[offset] & 0xFF) |
                ((buffer[offset+1] & 0xFF) << 8) |
                ((buffer[offset+2] & 0xFF) << 16) |
                ((buffer[offset+3] & 0xFF) << 24));
    }

    private static int readBytes(final InputStream stream, final byte[] buffer, final int offset, final int length) throws IOException {
        int bytesRead = 0;
        while (bytesRead < length) {
            final int count = stream.read(buffer, offset + bytesRead, length - bytesRead);

            // Return EOF if we get EOF from read() and we haven't read any bytes
            if ( count < 0 && bytesRead == 0 ) {
                return count;
            // Otherwise if we hit EOF and we have read something, return the bytes read
            } else if (count <= 0) {
                break;
            }

            bytesRead += count;
        }
        return bytesRead;
    }

    private static final class BGZFBlockMetadata {
        private final long blockOffset;
        private final int compressedSize;
        private final int uncompressedSize;

        public BGZFBlockMetadata(final long blockOffset, final int compressedSize, final int uncompressedSize) {
            this.blockOffset = blockOffset;
            this.compressedSize = compressedSize;
            this.uncompressedSize = uncompressedSize;
        }
    }
}
