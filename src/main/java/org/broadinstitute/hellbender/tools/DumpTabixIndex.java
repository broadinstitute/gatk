package org.broadinstitute.hellbender.tools;

import htsjdk.samtools.util.FileExtensions;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.exceptions.UserException;
import picard.cmdline.programgroups.DiagnosticsAndQCProgramGroup;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.GZIPInputStream;

@CommandLineProgramProperties(
        summary = "Prints a text file describing the contents of the tabix index input file.",
        oneLineSummary = "Dumps a tabix index file.",
        usageExample = "gatk DumpTabixIndex -I tabixIndex.tbi -O output.txt",
        programGroup = DiagnosticsAndQCProgramGroup.class
)
@DocumentedFeature
public class DumpTabixIndex extends CommandLineProgram {
    @Argument( doc = "Tabix index file.",
            fullName = "tabix-index",
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME)
    GATKPath tabixIndexFile;

    @Argument( doc = "Output file.",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME)
    GATKPath outputFile;

    @Override
    protected Object doWork() {
        if ( !tabixIndexFile.hasExtension(FileExtensions.TABIX_INDEX) ) {
            throw new UserException("Expected a " + FileExtensions.TABIX_INDEX + " file as input.");
        }

        final PrintStream writer = new PrintStream(outputFile.getOutputStream());
        try ( final PushbackInputStream is =
                      new PushbackInputStream(new GZIPInputStream(tabixIndexFile.getInputStream())) ) {
            dumpTabixIndex(is, writer);
        } catch ( final IOException ioe ) {
            throw new UserException("Trouble reading index.", ioe);
        }
        if ( writer.checkError() ) {
            throw new UserException("Trouble writing output.");
        }

        return null;
    }

    public static void dumpTabixIndex( final PushbackInputStream is,
                                       final PrintStream writer ) throws IOException {
        if ( readByte(is) != 'T' || readByte(is) != 'B' ||
                readByte(is) != 'I' || readByte(is) != 1 ) {
            throw new UserException("Incorrect magic number for tabix index");
        }
        final int nTigs = readInt(is);
        final int format = readInt(is);
        final int seqCol = readInt(is);
        final int begCol = readInt(is);
        final int endCol = readInt(is);
        final char meta = (char)readInt(is);
        final int skip = readInt(is);
        final int namesLen = readInt(is);
        writer.println("#tigs\tformat\tseqCol\tbegCol\tendCol\tmetaChr\tskip");
        writer.println(nTigs + "\t" + format + "\t" + seqCol + "\t" + begCol + "\t" +
                endCol + "\t" + meta + "\t" + skip + "\n");
        final List<String> tigNames = readContigNames(is, nTigs, namesLen);
        for ( final String tigName : tigNames ) {
            writer.println(tigName + " binned index:");
            int nBins = readInt(is);
            while ( nBins-- > 0 ) {
                final int binNo = readInt(is);
                int nChunks = readInt(is);
                if ( binNo > 37448 ) {
                    if ( nChunks != 2 ) {
                        throw new UserException("pseudobin has " + nChunks + " chunks");
                    }
                    final long chunkStart = readLong(is);
                    final long chunkEnd = readLong(is);
                    final long nMapped = readLong(is);
                    final long nUnmapped = readLong(is);
                    writer.println(tigName + " summary: mapped=" + nMapped +
                            "\tplaced=" + nUnmapped +
                            "\tstart=" + Long.toHexString(chunkStart >>> 16) +
                            ":" + Long.toHexString(chunkStart & 0xffff) +
                            "\tend=" + Long.toHexString(chunkEnd >>> 16) +
                            ":" + Long.toHexString(chunkStart & 0xffff));
                    continue;
                }
                if ( binNo == 0 ) {
                    writer.print(binNo + "\t" + tigName + ":0M-512M\t");
                } else if ( binNo <= 8 ) {
                    final int binStart = (binNo - 1)*64;
                    writer.print(binNo + "\t" + tigName + ":" + binStart + "M-" + (binStart+64) + "M\t");
                } else if ( binNo <= 72 ) {
                    final int binStart = (binNo - 9)*8;
                    writer.print(binNo + "\t" + tigName + ":" + binStart + "M-" + (binStart+8) + "M\t");
                } else if ( binNo <= 584 ) {
                    final int binStart = binNo - 73;
                    writer.print(binNo + "\t" + tigName + ":" + binStart + "M-" + (binStart + 1) + "M\t");
                } else if ( binNo <= 4680 ) {
                    final int binStart = (binNo - 585)*128;
                    writer.print(binNo + "\t" + tigName + ":" + binStart + "K-" + (binStart + 128) + "K\t");
                } else {
                    final int binStart = (binNo - 4681)*16;
                    writer.print(binNo + "\t" + tigName + ":" + binStart + "K-" + (binStart + 16) + "K\t");
                }
                while ( nChunks-- > 0 ) {
                    final long chunkStart = readLong(is);
                    final long chunkEnd = readLong(is);
                    writer.print("\t" + Long.toHexString(chunkStart >>> 16) +
                            ":" + Long.toHexString(chunkStart & 0xffff) +
                            "->" + Long.toHexString(chunkEnd >>> 16) +
                            ":" + Long.toHexString(chunkEnd & 0xffff));
                }
                writer.println();
            }
            int nIntervals = readInt(is);
            int nK = 0;
            writer.println();
            writer.println(tigName + " linear index:");
            while ( nIntervals-- > 0 ) {
                final long chunkOffset = readLong(is);
                writer.println(nK + "K\t" + Long.toHexString(chunkOffset >>> 16) +
                        ":" + Long.toHexString(chunkOffset & 0xffff));
                nK += 16;
            }
        }
        int nextByte;
        if ( (nextByte = is.read()) != -1 ) {
            is.unread(nextByte);
            final long nUnmapped = readLong(is);
            writer.println(nUnmapped + " unplaced reads.");
        }
        if ( (nextByte = is.read()) != -1 ) {
            throw new UserException("Unexpected data follows index.");
        }
    }

    public static List<String> readContigNames( final InputStream is,
                                                final int nNames,
                                                final int namesLen ) throws IOException {
        final List<String> names = new ArrayList<>(nNames);
        final StringBuilder sb = new StringBuilder();
        int namesProcessed = 0;
        int nextByte;
        int nBytesRead = 0;
        while ( namesProcessed++ < nNames ) {
            while ( (nextByte = readByte(is)) != 0 ) {
                sb.append((char)nextByte);
                nBytesRead += 1;
            }
            nBytesRead += 1;
            names.add(sb.toString());
            sb.setLength(0);
        }
        if ( nBytesRead != namesLen ) {
            throw new UserException("Contig names didn't have the correct length.");
        }
        return names;
    }

    public static int readByte( final InputStream is ) throws IOException {
        final int value = is.read();
        if ( value == -1 ) {
            throw new IOException("Tried to read past EOF");
        }
        return value;
    }

    // tabix integers are little-endian
    public static int readShort( final InputStream is ) throws IOException {
        return readByte(is) | (readByte(is) << 8);
    }

    public static int readInt( final InputStream is ) throws IOException {
        return readShort(is) | (readShort(is) << 16);
    }

    public static long readLong( final InputStream is ) throws IOException {
        return (readInt(is) & 0xffffffffL) | ((long)readInt(is) << 32);
    }
}
