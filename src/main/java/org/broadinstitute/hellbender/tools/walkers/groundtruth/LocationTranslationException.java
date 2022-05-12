package org.broadinstitute.hellbender.tools.walkers.groundtruth;

import org.broadinstitute.hellbender.exceptions.GATKException;

public class LocationTranslationException extends GATKException {

    static final private long        serialVersionUID = 0;

    LocationTranslationException(String msg) {

        super(msg);
    }
}
