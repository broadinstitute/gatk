package org.broadinstitute.hellbender.tools.picard.illumina.parser;

/**
 * List of data types of interest when parsing Illumina data.  Because different Illumina versions
 * splatter data of these types across different files, by specifying only the data types of interest,
 * the number of files read can be reduced.
 *
 * @author jburke@broadinstitute.org
 */
public enum IlluminaDataType {
    Position, BaseCalls, QualityScores, PF, Barcodes
}
