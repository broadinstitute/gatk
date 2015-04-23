package org.broadinstitute.hellbender.tools.walkers.vqsr;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.exceptions.GATKException;

public final class VQSRUtils {
    private VQSRUtils() {/*no instance*/}

    public static boolean checkVariationClass( final VariantContext evalVC, final VariantRecalibratorArgumentCollection.Mode mode ) {
        switch( mode ) {
            case SNP:
                return evalVC.isSNP() || evalVC.isMNP();
            case INDEL:
                return evalVC.isStructuralIndel() || evalVC.isIndel() || evalVC.isMixed() || evalVC.isSymbolic();
            case BOTH:
                return true;
            default:
                throw new GATKException( "Encountered unknown recal mode: " + mode );
        }
    }
}
