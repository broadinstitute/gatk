package org.broadinstitute.hellbender.tools;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import picard.cmdline.programgroups.DiagnosticsAndQCProgramGroup;

import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.GZIPInputStream;

@CommandLineProgramProperties(
        summary = "Prints a tabix index file as text.",
        oneLineSummary = "Dumps a tabix index file.",
        programGroup = DiagnosticsAndQCProgramGroup.class
)
@ExperimentalFeature
public class DumpTabixIndex extends CommandLineProgram {
    @Argument( doc = "Tabix index file.",
            fullName = "tabix-index",
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME)
    String tabixIndexFile;

    @Override
    protected Object doWork() {
        try ( final InputStream is = new GZIPInputStream(new FileInputStream(tabixIndexFile)) ) {
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
            System.out.println("#tigs\tformat\tseqCol\tbegCol\tendCol\tmetaChr\tskip");
            System.out.println(nTigs + "\t" + format + "\t" + seqCol + "\t" + begCol + "\t" +
                    endCol + "\t" + meta + "\t" + skip + "\n");
            final List<String> tigNames = readContigNames(is, nTigs, namesLen);
            for ( final String tigName : tigNames ) {
                System.out.println(tigName + " binned index:");
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
                        System.out.println(tigName + " summary: mapped=" + nMapped +
                                "\tplaced=" + nUnmapped +
                                "\tstart=" + Long.toHexString(chunkStart >>> 16) +
                                ":" + Long.toHexString(chunkStart & 0xffff) +
                                "\tend=" + Long.toHexString(chunkEnd >>> 16) +
                                ":" + Long.toHexString(chunkStart & 0xffff));
                        continue;
                    }
                    if ( binNo == 0 ) {
                        System.out.print(binNo + "\t" + tigName + ":0M-512M\t");
                    } else if ( binNo <= 8 ) {
                        final int binStart = (binNo - 1)*64;
                        System.out.print(binNo + "\t" + tigName + ":" + binStart + "M-" + (binStart+64) + "M\t");
                    } else if ( binNo <= 72 ) {
                        final int binStart = (binNo - 9)*8;
                        System.out.print(binNo + "\t" + tigName + ":" + binStart + "M-" + (binStart+8) + "M\t");
                    } else if ( binNo <= 584 ) {
                        final int binStart = binNo - 73;
                        System.out.print(binNo + "\t" + tigName + ":" + binStart + "M-" + (binStart + 1) + "M\t");
                    } else if ( binNo <= 4680 ) {
                        final int binStart = (binNo - 585)*128;
                        System.out.print(binNo + "\t" + tigName + ":" + binStart + "K-" + (binStart + 128) + "K\t");
                    } else {
                        final int binStart = (binNo - 4681)*16;
                        System.out.print(binNo + "\t" + tigName + ":" + binStart + "K-" + (binStart + 16) + "K\t");
                    }
                    while ( nChunks-- > 0 ) {
                        final long chunkStart = readLong(is);
                        final long chunkEnd = readLong(is);
                        System.out.print("\t" + Long.toHexString(chunkStart >>> 16) +
                                ":" + Long.toHexString(chunkStart & 0xffff) +
                                "->" + Long.toHexString(chunkEnd >>> 16) +
                                ":" + Long.toHexString(chunkEnd & 0xffff));
                    }
                    System.out.println();
                }
                int nIntervals = readInt(is);
                int nK = 0;
                System.out.println();
                System.out.println(tigName + " linear index:");
                while ( nIntervals-- > 0 ) {
                    final long chunkOffset = readLong(is);
                    System.out.println(nK + "K\t" + Long.toHexString(chunkOffset >>> 16) +
                            ":" + Long.toHexString(chunkOffset & 0xffff));
                    nK += 16;
                }
            }
            if ( is.available() > 0 ) {
                final long nUnmapped = readLong(is);
                System.out.println(nUnmapped + " unplaced reads.");
            }
            if ( is.available() != 0 ) {
                throw new UserException("Unexpected data follows index.");
            }
        } catch ( final IOException ioe ) {
            throw new UserException("Can't read " + tabixIndexFile + " as gzip input stream", ioe);
        }
        return null;
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
