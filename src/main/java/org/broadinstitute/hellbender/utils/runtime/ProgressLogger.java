package org.broadinstitute.hellbender.utils.runtime;

import htsjdk.samtools.util.AbstractProgressLogger;
import org.apache.logging.log4j.Logger;

import java.text.DecimalFormat;
import java.text.NumberFormat;

/**
 * Facilitate consistent logging output when progressing through a stream of SAM records.
 *
 * Implements the SAMTools public htsjdk.samtools.util.AbstractProgressLogger interface to enable integration with SAMWriters.
 *
 */
public class ProgressLogger extends AbstractProgressLogger {

    private final Logger logger;

    /**
     * Construct a progress logger.
     * @param logger the Logger object to write output to
     * @param n the frequency with which to output (i.e. every N records)
     * @param verb the verb to log, e.g. "Processed, Read, Written".
     * @param noun the noun to use when logging, e.g. "Records, Variants, Loci"
     */
    public ProgressLogger(final Logger logger, final int n, final String verb, final String noun) {
        super(noun, verb, n);
        this.logger = logger;
    }

    /**
     * Construct a progress logger.
     * @param logger the Logger object to write outputs to
     * @param n the frequency with which to output (i.e. every N records)
     * @param verb the verb to log, e.g. "Processed, Read, Written".
     */
    public ProgressLogger(final Logger logger, final int n, final String verb) {
        this(logger, n, verb, "records");
    }

    /**
     * Construct a progress logger with the desired log and frequency and the verb "Processed".
     * @param logger the Logger object to write outputs to
     * @param n the frequency with which to output (i.e. every N records)
     */
    public ProgressLogger(final Logger logger, final int n) { this(logger, n, "Processed"); }

    /**
     * Construct a progress logger with the desired log, the verb "Processed" and a period of 1m records.
     * @param logger the Logger object to write outputs to
     */
    public ProgressLogger(final Logger logger) { this(logger, 1000000); }

    /**
     * Log a message to whatever logger is being used
     *
     * @param message a message to be logged by the logger (recommended output level is INFO or the equivalent)
     */
    @Override
    protected void log(String... message) {
        StringBuilder stringBuilder = new StringBuilder();
        for (String s : message) {
            stringBuilder.append(s);
        }
        logger.info(stringBuilder.toString());
    }
}
