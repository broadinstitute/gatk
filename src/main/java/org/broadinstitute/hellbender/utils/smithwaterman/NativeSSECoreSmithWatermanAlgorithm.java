package org.broadinstitute.hellbender.utils.smithwaterman;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.TextCigarCodec;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import ssw.Aligner;
import ssw.Alignment;

/**
 * Adapter for the Native vectorized SSE implementation of SmithWaterman
 */
public final class NativeSSECoreSmithWatermanAlgorithm implements CoreSmithWatermanAlgorithm{

    private static final Logger logger = LogManager.getLogger(NativeSSECoreSmithWatermanAlgorithm.class);

    private static boolean loadedOK; //Note: Can't make it final because compiler complains in the initializer block

    static {
        try {
            System.loadLibrary("sswjni");
            loadedOK = true;
        } catch (final UnsatisfiedLinkError e) {
            loadedOK = false;
            logger.warn(String.format("Cannot find libsswjni.so. Has the library been built and LD_LIBRARY_PATH or -Djava.library.path set appropriately?\n%s", e));
            logger.warn("java.library.path=" + System.getProperty("java.library.path"));
            logger.warn("LD_LIBRARY_PATH=" + System.getenv("LD_LIBRARY_PATH"));
        }
    }

    private Alignment result;
    private SmithWatermanParameters params;
    private int[][] matrix;

    @Override
    public boolean initialize(final SmithWatermanParameters params) {
        this.params = params;
        this.matrix = params.toSubstitutionMatrix();
        return true;//FIXME
    }

    @Override
    public void align(byte[] reference, byte[] alternate) {
        //Note the order swap of ref and alt for this implementation
        //Note that the sign of penalties is different in this implementation
        this.result = Aligner.align(alternate, reference, this.matrix, -1 * params.w_open, -1 * params.w_extend, false);
    }

    @Override
    public Cigar getCigar() {
        return TextCigarCodec.decode(result.cigar);
    }

    @Override
    public int getAlignmentStart2wrt1() {
        return result.ref_begin1;
    }
}
