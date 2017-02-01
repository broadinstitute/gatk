package org.broadinstitute.hellbender.utils.fermi;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.NativeUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

public class FermiLiteAssembler implements Assembler {
    private static volatile boolean nativeLibLoaded = false;

    public FermiLiteAssembler() {
        loadNativeLibrary();
    }

    public Assembly createAssembly( final Collection<GATKRead> reads ) {
        final ByteBuffer assemblyData = createAssemblyData(makeReadData(reads));
        if ( assemblyData == null ) throw new IllegalStateException("Unable to create assembly.");
        try {
            return interpretAssemblyData(assemblyData);
        } finally {
            destroyAssemblyData(assemblyData);
        }
    }

    private static void loadNativeLibrary() {
        if ( !nativeLibLoaded ) {
            synchronized (FermiLiteAssembler.class) {
                if ( !nativeLibLoaded ) {
                    final String libNameOverride = System.getProperty("LIBFML_PATH");
                    final String libName;
                    if ( libNameOverride != null ) libName = libNameOverride;
                    else if ( NativeUtils.runningOnPPCArchitecture() ) libName = null;
                    else if (NativeUtils.runningOnMac()) libName = "/libfml.Darwin.dylib";
                    else if (NativeUtils.runningOnLinux()) libName = "/libfml.Linux.so";
                    else libName = null;
                    if ( libName == null ) {
                        throw new UserException.HardwareFeatureException("We have a JNI binding for fermi-lite only for Linux and Mac.");
                    }
                    if ( !NativeUtils.loadLibraryFromClasspath(libName) ) {
                        throw new UserException.HardwareFeatureException("Misconfiguration: Unable to load fermi-lite native library "+libName);
                    }
                }
            }
        }
    }

    // Writes the number of reads, and then seqs and quals for each into a ByteBuffer.
    private static ByteBuffer makeReadData( final Collection<GATKRead> reads ) {
        // capacity calculation:  for each GATKRead we need two bytes (one for the base, one for the qual)
        //  for the length of the GATKRead (plus 1 for the null terminators), plus 4 bytes to give the array length.
        int capacity = reads.stream().mapToInt(GATKRead -> 2*(GATKRead.getLength()+1)).sum() + 4;
        final ByteBuffer readData = ByteBuffer.allocateDirect(capacity);
        readData.order(ByteOrder.nativeOrder());
        readData.putInt(reads.size()); // array length
        for ( final GATKRead read : reads ) {
            readData.put(read.getBases());
            readData.put((byte)0);
            readData.put(read.getBaseQualities());
            readData.put((byte)0);
        }
        readData.flip();
        return readData;
    }

    // expects a direct ByteBuffer containing:
    //  the number of contigs (4-byte int)
    //  the offset to the beginning of a byte pool containing sequence and per-base coverage bytes (4 byte int)
    //  N.B.: the sequence and per-base coverage bytes are NOT null terminated.
    //  for each contig, an fml_utg_t structure minus the seq and cov pointers:
    //    the length of the sequence (and per-base support) data (4-byte int)
    //    the number of supporting reads (4-byte int)
    //    the number of connections (4-byte int)
    //    a variable number (given by # of connections, above) of fml_ovlp_t's (8 bytes each)
    //  a byte pool containing the seq and cov data
    private static Assembly interpretAssemblyData( final ByteBuffer assemblyData ) {
        assemblyData.order(ByteOrder.nativeOrder()).position(0).limit(assemblyData.capacity());

        // make the contigs
        final int nContigs = assemblyData.getInt();
        int seqOffset = assemblyData.getInt();
        final List<Assembly.Contig> contigs = new ArrayList<>(nContigs);
        for ( int idx = 0; idx != nContigs; ++idx ) {
            final int seqLen = assemblyData.getInt();
            final int nSupportingReads = assemblyData.getInt();
            final int nConnections = assemblyData.getInt();
            final int mark = assemblyData.position()+8*nConnections; // sizeof(fml_ovlp_t) is 8
            assemblyData.position(seqOffset);
            final byte[] seq = new byte[seqLen];
            assemblyData.get(seq);
            final byte[] coverage = new byte[seqLen];
            assemblyData.get(coverage);
            contigs.add(new Assembly.Contig(seq,nSupportingReads,coverage));
            assemblyData.position(mark);
            seqOffset += 2*seqLen;
        }
        // connect the contigs
        assemblyData.position(8); // skip past nContigs and seqOffset
        for ( int idx = 0; idx != nContigs; ++idx ) {
            final Assembly.Contig contig = contigs.get(idx);
            assemblyData.getInt(); // skip seqLen
            assemblyData.getInt(); // skip # of supporting reads
            int nConnections = assemblyData.getInt();
            final List<Assembly.Connection> connections = new ArrayList<>(nConnections);
            while ( nConnections-- > 0 ) {
                int overlapLen = assemblyData.getInt();
                final boolean isRC = overlapLen < 0;
                overlapLen = overlapLen << 1 >> 1;
                int contigId = assemblyData.getInt();
                final boolean isTargetRC = contigId < 0;
                contigId = contigId << 1 >> 1;
                connections.add(new Assembly.Connection(contigs.get(contigId), overlapLen, isRC, isTargetRC));
            }
            contig.setConnections(connections);
        }

        return new Assembly(contigs);
    }

    // these should be called in succession by the same thread
    private static native ByteBuffer createAssemblyData( final ByteBuffer readData );
    private static native void destroyAssemblyData( final ByteBuffer assemblyData );
}
