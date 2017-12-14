package org.broadinstitute.hellbender.utils.smithwaterman;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.testng.SkipException;
import java.util.*;


public class SmithWatermanIntelAlignerUnitTest extends SmithWatermanAlignerAbstractUnitTest {

    /*
    *Test the Intel Aligner Native Interface
    */


    @Override
    protected SmithWatermanIntelAligner getAligner() {

        boolean loaded = true;
        SmithWatermanIntelAligner aligner = null;
        try {
            aligner = new SmithWatermanIntelAligner();
        }
        catch (UserException.HardwareFeatureException e ) {
            loaded = false;
        }
        if(!loaded) {
            //Skip test if correct version of AVX is not supported
            throw new SkipException("AVX SmithWaterman is not supported on this system or the library is not available");
        }
        return aligner;
    }

}
