package org.broadinstitute.hellbender.utils.smithwaterman;

import com.intel.gkl.smithwaterman.IntelSmithWaterman;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWAlignmentResult;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWAlignerNativeBinding;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWAlignerNativeArguments;
import org.broadinstitute.hellbender.exceptions.UserException;
import htsjdk.samtools.*;

/**
 * create a GKL smithwaterman
 *
 * @param implementation    which implementation to use (AVX or OMP)
 * @param args              arguments to the native GKL implementation
 */
public class SWAlignerAVX extends SWAligner {

    protected SWAlignerNativeArguments args;

    protected TextCigarCodec cigarCodec;


    public SWAlignerAVX(SWAlignerNativeArguments args) throws UserException.HardwareFeatureException {
        this.args = args;
    }

    /**
     * Aligns the alternate sequence to the reference sequence
     *
     * @param ref reference sequence
     * @param alt alternate sequence
     */

    public SWAligner.SWAlignerResult align(final byte[] ref, final byte[] alt)
    {
        final boolean isSupported;
        final SWAlignmentResult alignResult;
        IntelSmithWaterman smithwaterman = new IntelSmithWaterman();

        isSupported =  smithwaterman.load(null);
        if (!isSupported) {
            throw new UserException.HardwareFeatureException("Machine does not support AVX SmithWaterman");
        }

        smithwaterman.initialize(args);
        alignResult = smithwaterman.align(ref,alt);

        /*** code converts the cigar string to a cigar **/

        SWAlignerResult result = new SWAlignerResult(cigarCodec.decode(alignResult.cigar), alignResult.alignment_offset);


        return result;
    }
}
