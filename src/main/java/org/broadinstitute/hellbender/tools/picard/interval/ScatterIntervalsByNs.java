package org.broadinstitute.hellbender.tools.picard.interval;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.StringUtil;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.PicardCommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.IntervalProgramGroup;

import java.io.File;
import java.util.*;

/**
 * A CLP for breaking up a reference into intervals of Ns and ACGTs bases.
 * Used for creating a broken-up interval list for calling WGS.
 *
 * @author Yossi Farjoun
 */
@CommandLineProgramProperties(
        usage = "Writes an interval list based on splitting the reference by Ns.",
        usageShort = "Writes an interval list based on splitting the reference by Ns",
        programGroup = IntervalProgramGroup.class
)
public class ScatterIntervalsByNs extends PicardCommandLineProgram {

    @Argument(shortName = StandardArgumentDefinitions.REFERENCE_SHORT_NAME, doc = "Reference sequence to use.")
    public File REFERENCE;

    @Argument(shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "Output file for interval list.")
    public File OUTPUT;

    @Argument(shortName = "OT", doc = "Type of intervals to output.", optional = true)
    public OutputType OUTPUT_TYPE = OutputType.BOTH;

    @Argument(shortName = "N", doc = "Maximal number of contiguous N bases to tolerate, thereby continuing the current ACGT interval.", optional = true)
    public int MAX_TO_MERGE = 1;

    //not using an enum since Interval.name is a String, and am using that to define the type of the Interval
    static final String
            ACGTmer = "ACGTmer",
            Nmer    = "Nmer";

    //an enum to determine which types of intervals get outputted
    private enum OutputType {
        N(Nmer),
        ACGT(ACGTmer),
        BOTH(Nmer, ACGTmer);

        private final Set acceptedTypes;

        public Boolean accepts(final String string) {return acceptedTypes.contains(string);}

        OutputType(final String... strings) {
            acceptedTypes = new HashSet<String>();
            Collections.addAll(acceptedTypes, strings);
        }
    }

    private static final Log log = Log.getInstance(ScatterIntervalsByNs.class);
    final ProgressLogger locusProgress = new ProgressLogger(log, (int) 1e7, "examined", "loci");
    final ProgressLogger intervalProgress = new ProgressLogger(log, (int) 10, "found", "intervals");

    @Override
    protected Object doWork() {
        IOUtil.assertFileIsReadable(REFERENCE);
        IOUtil.assertFileIsWritable(OUTPUT);

        final ReferenceSequenceFile refFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(REFERENCE, true);

        // get the intervals
        final IntervalList intervals = segregateReference(refFile, MAX_TO_MERGE);

        log.info(String.format("Found %d intervals in %d loci during %s seconds", intervalProgress.getCount(), locusProgress.getCount(), locusProgress.getElapsedSeconds()));

        /**********************************
         * Now output regions for calling *
         **********************************/

        final IntervalList outputIntervals = new IntervalList(intervals.getHeader().clone());
        log.info(String.format("Collecting requested type of intervals (%s)", OUTPUT_TYPE));

        for (final Interval i : intervals.getIntervals()) {
            if (OUTPUT_TYPE.accepts(i.getName())) {
                outputIntervals.add(i);
            }
        }

        log.info("Writing Intervals.");
        outputIntervals.write(OUTPUT);

        log.info(String.format("Execution ending. Total time %d seconds", locusProgress.getElapsedSeconds()));

        return null;
    }

    /**
     * ****************************************************************
     * Generate an interval list that alternates between Ns and ACGTs *
     * ****************************************************************
     */
    public static IntervalList segregateReference(final ReferenceSequenceFile refFile, final int maxNmerToMerge) {
        final List<Interval> preliminaryIntervals = new LinkedList<>();
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(refFile.getSequenceDictionary());
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        final IntervalList finalIntervals = new IntervalList(header);

        //iterate over all the sequences in the dictionary
        for (final SAMSequenceRecord rec : refFile.getSequenceDictionary().getSequences()) {
            final ReferenceSequence ref = refFile.getSequence(rec.getSequenceName());
            final byte[] bytes = ref.getBases();
            StringUtil.toUpperCase(bytes);

            boolean nBlockIsOpen = (bytes[0] == 'N');
            int start = 0;

            for (int i = 0; i < bytes.length; ++i) {
                final boolean currentBaseIsN = (bytes[i] == 'N');

                //create intervals when switching, i.e "nBlockIsOpen" disagrees with "currentBaseIsN"
                if (nBlockIsOpen != currentBaseIsN) {
                    preliminaryIntervals.add(new Interval(rec.getSequenceName(), start + 1, i, false, nBlockIsOpen ? Nmer : ACGTmer));
                    start = i;
                    nBlockIsOpen = !nBlockIsOpen;
                }
            }
            // Catch the last block of chromosome
            preliminaryIntervals.add(new Interval(rec.getSequenceName(), start + 1, bytes.length, false, nBlockIsOpen ? Nmer : ACGTmer));
        }

        // now that we have the whole list, we need to remove the short Nmers.
        // process the list, replacing trios with short Nmers in the middle with longer intervals:
        while (!preliminaryIntervals.isEmpty()) {

            //if top trio match the bill, replace them with a merged interval,
            // and push it back the top of the list (we expect alternating Nmers and ACGTmers, but
            // not assuming it in the logic)

            //(I want this to be fast and the strings are all copies of the static prototypes Nmer and ACGTmer )
            //noinspection StringEquality
            if (preliminaryIntervals.size() >= 3 && // three or more intervals
                    preliminaryIntervals.get(0).getName() == ACGTmer &&   //an N-mer
                    preliminaryIntervals.get(1).getName() == Nmer &&      //between two
                    preliminaryIntervals.get(2).getName() == ACGTmer &&   //ACGT-mers
                    preliminaryIntervals.get(0).abuts(preliminaryIntervals.get(1)) && // all abutting
                    preliminaryIntervals.get(1).abuts(preliminaryIntervals.get(2)) && // each other (there are many contigs...)
                    preliminaryIntervals.get(1).length() <= maxNmerToMerge) //and the N-mer is of length N or less
            {
                // create the new ACGTmer interval
                final Interval temp = new Interval(
                        preliminaryIntervals.get(0).getContig(),
                        preliminaryIntervals.get(0).getStart(),
                        preliminaryIntervals.get(2).getEnd(), false, ACGTmer);

                //remove the first 3 elements of the list
                for (int i = 0; i < 3; ++i) {
                    preliminaryIntervals.remove(0);
                }
                //and replace them with the newly created one
                preliminaryIntervals.add(0, temp);
            } else { //if cannot merge top three intervals, transfer the top intervals to finalIntervals
                finalIntervals.add(preliminaryIntervals.remove(0));
            }
        }
        return finalIntervals;
    }
}
