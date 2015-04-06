package org.broadinstitute.hellbender.tools.picard.illumina.parser;

import org.broadinstitute.hellbender.exceptions.GATKException;

public class IlluminaParserException extends GATKException {
    public IlluminaParserException(String message) {
        super(message);
    }

    public IlluminaParserException(String message, Throwable throwable) {
        super(message, throwable);
    }
}
