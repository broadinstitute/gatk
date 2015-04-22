package org.broadinstitute.hellbender.tools.picard.illumina.parser.readers;

import org.broadinstitute.hellbender.exceptions.GATKException;

public final class IlluminaReaderException extends GATKException {
    public IlluminaReaderException(String message) {
        super(message);
    }

    public IlluminaReaderException(String message, Throwable throwable) {
        super(message, throwable);
    }
}
