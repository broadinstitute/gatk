package org.broadinstitute.hellbender.tools.walkers.misc;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CoverageAnalysisProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.iterators.IntervalLocusIterator;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Created by tsato on 1/24/18.
 */
@CommandLineProgramProperties(
        summary = "",
        oneLineSummary = "",
        programGroup = CoverageAnalysisProgramGroup.class
)
public class CreateHomopolymerIntervals extends GATKTool {
    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "")
    private File outputFile = null;

    private PrintStream outputStream = null;

    private static final int MIN_HOMOPOLYMER_LENGTH = 5;

    private static final int LOOK_AHEAD_LENGTH = 100;

    @Override
    public void onTraversalStart() {
        try {
            outputStream = outputFile != null ? new PrintStream(outputFile) : System.out;
        }
        catch ( final FileNotFoundException e ) {
            throw new UserException.CouldNotReadInputFile(outputFile, e);
        }
    }

    @Override
    public void traverse() {
        final ReferenceDataSource reference = new ReferenceFileSource(referenceArguments.getReferencePath());

        final SAMSequenceDictionary sequenceDictionary = getBestAvailableSequenceDictionary();
        final List<SimpleInterval> intervals = new ArrayList<>();

        for (SAMSequenceRecord record : sequenceDictionary.getSequences()){
            final String contig = record.getSequenceName();
            final int length = record.getSequenceLength();
            intervals.add(new SimpleInterval(contig, 1, length));
        }

        IntervalLocusIterator iterator = new IntervalLocusIterator(intervals.iterator());

        while (iterator.hasNext()){
            final SimpleInterval locus = iterator.next();
            final ReferenceContext referenceContext = new ReferenceContext(reference, locus);
            if (referenceContext.getBase() == 'N'){
                referenceContext.setWindow(0, LOOK_AHEAD_LENGTH);
                final byte[] upcomingBases = referenceContext.getBases();

                if (upcomingBases.length == LOOK_AHEAD_LENGTH + 1 && upcomingBases[LOOK_AHEAD_LENGTH] == 'N'){
                    // we're in an ocean of N's - skip ahead this many bases
                    for (int i = 0; i < LOOK_AHEAD_LENGTH; i++){
                        iterator.next();
                    }
                }

                continue;
            } else {
                final byte baseAtLocus = referenceContext.getBase();

                final byte[] precedingBases = getPrecedingBases(referenceContext, MIN_HOMOPOLYMER_LENGTH);
                final byte[] trailingBases = getTrailingBases(referenceContext, MIN_HOMOPOLYMER_LENGTH);

                Utils.validate(precedingBases.length == MIN_HOMOPOLYMER_LENGTH,
                        String.format("preceding bases must have the length %d but got %d",
                                MIN_HOMOPOLYMER_LENGTH, precedingBases.length));

                Utils.validate(trailingBases.length == MIN_HOMOPOLYMER_LENGTH,
                        String.format("trailing bases must have the length %d but got %d",
                                MIN_HOMOPOLYMER_LENGTH, trailingBases.length));

                if (neighboringBasesAreAHomopolymer(precedingBases, baseAtLocus)){
                    outputStream.println(String.format("%s\t%s\t%s", locus.getContig(), locus.getStart(), locus.getEnd()));
                } else if (neighboringBasesAreAHomopolymer(trailingBases, baseAtLocus)){
                    outputStream.println(String.format("%s\t%s\t%s", locus.getContig(), locus.getStart(), locus.getEnd()));
                }
            }
        }
    }

    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    public void closeTool() {
        if ( outputStream != null ) {
            outputStream.close();
        }
    }


    public static byte[] getTrailingBases(final ReferenceContext referenceContext, final int numberOfBases){
        referenceContext.setWindow(0, MIN_HOMOPOLYMER_LENGTH);
        return Arrays.copyOfRange(referenceContext.getBases(), 1, MIN_HOMOPOLYMER_LENGTH + 1);

    }

    public static byte[] getPrecedingBases(final ReferenceContext referenceContext, final int numberOfBases){
        referenceContext.setWindow(numberOfBases, 0);
        return Arrays.copyOf(referenceContext.getBases(), MIN_HOMOPOLYMER_LENGTH);

    }

    private boolean neighboringBasesAreAHomopolymer(final byte[] neighboringBases, final byte baseAtLocus){
        Utils.nonNull(neighboringBases, "the array of bases cannot be null");

        if (baseAtLocus == neighboringBases[0]){
            // If the current base matches any of the preceding bases, then it doesn't immediately follow a homopolymer
            return false;
        }

        final byte nucleotide = neighboringBases[0];
        for (byte precedingBase : neighboringBases){
            if (precedingBase != nucleotide || precedingBase == 'N'){
                return false;
            }
        }
        return true;
    }

}
