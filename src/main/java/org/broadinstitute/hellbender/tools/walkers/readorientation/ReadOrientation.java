package org.broadinstitute.hellbender.tools.walkers.readorientation;

/**
 * Created by tsato on 3/28/18.
 */

public enum ReadOrientation {
    F1R2,
    F2R1;

    public static ReadOrientation getOtherOrientation(final ReadOrientation orientation){
        return orientation == F1R2 ? F2R1 : F1R2;
    }

    public static int SIZE = ReadOrientation.values().length;

}

