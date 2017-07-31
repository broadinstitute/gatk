package org.broadinstitute.hellbender.utils.smithwaterman;


import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import com.intel.gkl.smithwaterman.IntelSmithWaterman;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWAlignerNativeArguments;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWAlignmentResult;
import org.broadinstitute.hellbender.exceptions.UserException;
import java.util.Map;
import java.util.function.Function;

/**
 * Pairwise discrete smith-waterman alignment
 *
 * ************************************************************************
 * ****                    IMPORTANT NOTE:                             ****
 * ****  This class assumes that all bytes come from UPPERCASED chars! ****
 * ************************************************************************
 */

public class SWAligner {

    protected static final Logger logger = LogManager.getLogger(SWAligner.class);

    protected SWAlignerNativeArguments args;


    public static final class SWAlignerResult {
        public Cigar cigar;
        public int alignment_offset;
        SWAlignerResult(final Cigar cigar, final int alignment_offset) {
            this.cigar = cigar;
            this.alignment_offset = alignment_offset;
        }

    }


    public enum Implementation {
        FASTEST_AVAILABLE,
        AVX_SMITHWATERMAN,
        ORIGINAL_SMITHWATERMAN,
    }

    public static final SWAlignerNativeArguments ORIGINAL_DEFAULT = new SWAlignerNativeArguments(SWAlignerNativeArguments.OverhangStrategy.SOFTCLIP,3,-1,-4,-3);

    public SWAligner() {
        args = ORIGINAL_DEFAULT;
    }

    public SWAligner(SWAlignerNativeArguments args) {
        this.args = args;
    }

    public SWAligner(SWAlignerNativeArguments args, Implementation implementation) {
        this.args = args;
    }

    public SWAligner getAligner(SWAlignerNativeArguments args, Implementation implementation) {

        /* Selection of SmithWaterman implemetation */
        switch (implementation) {

            case FASTEST_AVAILABLE:{
                  /*  AVX2 implementation of SmithWaterman */
                try {
                    final SWAlignerAVX swAligner = new SWAlignerAVX(args);
                    logger.info("AVX optmized smithwaterman implementation");
                    return swAligner;
                }
                catch  (UserException.HardwareFeatureException e ){
                    logger.info("AVX smithwaterman not supported");
                }
                return new SWAlignerJava(args);
            }
            case AVX_SMITHWATERMAN: {
                 /*  AVX2 implementation of SmithWaterman */
                try {
                    final SWAlignerAVX swAligner = new SWAlignerAVX(args);
                    logger.info("AVX optmized smithwaterman implementation");
                    return swAligner;
                }
                catch  (UserException.HardwareFeatureException e ){
                    logger.info("AVX smithwaterman not supported");
                }

            }
        /* Oroginal slower Java  implementation of SmithWaterman */
            case ORIGINAL_SMITHWATERMAN: {

                final SWAlignerJava swAligner = new SWAlignerJava(args);
                logger.info("Using the slower smithwaterman implementation");
                return swAligner;
            }

            default:
                throw new UserException.HardwareFeatureException("Unknown SmithWaterman implementation.");
        }

    }
}

