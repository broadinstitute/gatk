package org.broadinstitute.hellbender.tools.picard.interval;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.liftover.LiftOver;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.Log;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.PicardCommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.IntervalProgramGroup;

import java.io.File;
import java.util.*;

import static htsjdk.samtools.SamReaderFactory.makeDefault;
import static htsjdk.samtools.liftover.LiftOver.DEFAULT_LIFTOVER_MINMATCH;
import static htsjdk.samtools.util.IOUtil.assertFileIsReadable;
import static htsjdk.samtools.util.IOUtil.assertFileIsWritable;
import static htsjdk.samtools.util.IntervalList.fromFile;
import static htsjdk.samtools.util.Log.getInstance;
import static org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions.INPUT_SHORT_NAME;
import static org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions.OUTPUT_SHORT_NAME;
import static org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions.SEQUENCE_DICTIONARY_SHORT_NAME;

/**
 * @author alecw@broadinstitute.org
 */
@CommandLineProgramProperties(
        usage = "Lifts over an interval list from one reference build to another. Based on UCSC liftOver." +
                " Uses a UCSC chain file to guide the liftOver.",
        usageShort = "Lifts over an interval list from one reference build to another",
        programGroup = IntervalProgramGroup.class
)
public final class LiftOverIntervalList extends PicardCommandLineProgram {

    private static final Log LOG = getInstance(LiftOverIntervalList.class);

    @Argument(doc = "Interval list to be lifted over.", shortName = INPUT_SHORT_NAME)
    public File INPUT;

    @Argument(doc = "Where to write lifted-over interval list.", shortName = OUTPUT_SHORT_NAME)
    public File OUTPUT;

    @Argument(doc = "Sequence dictionary to write into the output interval list.",
            shortName = SEQUENCE_DICTIONARY_SHORT_NAME)
    public File SEQUENCE_DICTIONARY;

    @Argument(doc = "Chain file that guides LiftOver.")
    public File CHAIN;

    @Argument(doc = "Minimum percentage of bases in each input interval that must map to output interval.")
    public double MIN_LIFTOVER_PCT = DEFAULT_LIFTOVER_MINMATCH;

    /**
     * Do the work after command line has been parsed. RuntimeException may be
     * thrown by this method, and are reported appropriately.
     */
    @Override
    protected Object doWork() {
        assertFileIsReadable(INPUT);
        assertFileIsReadable(SEQUENCE_DICTIONARY);
        assertFileIsReadable(CHAIN);
        assertFileIsWritable(OUTPUT);

        final LiftOver liftOver = new LiftOver(CHAIN);
        liftOver.setLiftOverMinMatch(MIN_LIFTOVER_PCT);

        final IntervalList fromIntervals = fromFile(INPUT);
        final SAMFileHeader toHeader = makeDefault().getFileHeader(SEQUENCE_DICTIONARY);
        liftOver.validateToSequences(toHeader.getSequenceDictionary());
        final IntervalList toIntervals = new IntervalList(toHeader);
        boolean anyFailed = false;
        for (final Interval fromInterval : fromIntervals) {
            final Interval toInterval = liftOver.liftOver(fromInterval);
            if (toInterval != null) {
                toIntervals.add(toInterval);
            } else {
                anyFailed = true;
                LOG.warn("Liftover failed for ", fromInterval, "(len ", fromInterval.length(), ")");
                final List<LiftOver.PartialLiftover> partials = liftOver.diagnosticLiftover(fromInterval);
                for (final LiftOver.PartialLiftover partial : partials) {
                    LOG.info(partial);
                }
            }
        }

        toIntervals.sorted();
        toIntervals.write(OUTPUT);
        return anyFailed ? 1 : 0;
    }
}
