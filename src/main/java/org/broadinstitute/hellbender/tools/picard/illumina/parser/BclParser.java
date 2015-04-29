package org.broadinstitute.hellbender.tools.picard.illumina.parser;

import htsjdk.samtools.util.CloseableIterator;
import org.broadinstitute.hellbender.tools.picard.illumina.parser.readers.BclQualityEvaluationStrategy;
import org.broadinstitute.hellbender.tools.picard.illumina.parser.readers.BclReader;

import java.io.File;
import java.util.*;

import static htsjdk.samtools.util.CollectionUtil.makeSet;
import static java.lang.Integer.MIN_VALUE;
import static java.lang.Math.max;
import static java.util.Collections.unmodifiableSet;
import static org.broadinstitute.hellbender.tools.picard.illumina.parser.IlluminaDataType.BaseCalls;
import static org.broadinstitute.hellbender.tools.picard.illumina.parser.IlluminaDataType.QualityScores;

/**
 * BclParser parses a number of BclFiles equal to the total of all the values in outputLengths and returns a BclData object
 * segmented based on these lengths.  The only client of this class should be IlluminaDataProvider and an test classes.  See BclReader for
 * more information on BclFiles.  BclParser provides support for reading BaseCalls and QualityScores.
 */
class BclParser extends PerTileCycleParser<BclData> {
    private static final int EAMSS_M2_GE_THRESHOLD = 30;
    private static final int EAMSS_S1_LT_THRESHOLD = 15; //was 15
    public static final byte MASKING_QUALITY = (byte) 0x02;

    private static final Set<IlluminaDataType> SUPPORTED_TYPES = unmodifiableSet(makeSet(BaseCalls, QualityScores));

    protected final BclQualityEvaluationStrategy bclQualityEvaluationStrategy;
    private final boolean applyEamssFilter;

    public BclParser(final File directory, final int lane, final CycleIlluminaFileMap tilesToCycleFiles, final OutputMapping outputMapping, final BclQualityEvaluationStrategy bclQualityEvaluationStrategy) {
        this(directory, lane, tilesToCycleFiles, outputMapping, true, bclQualityEvaluationStrategy);
        this.initialize();
    }

    public BclParser(final File directory, final int lane, final CycleIlluminaFileMap tilesToCycleFiles, final OutputMapping outputMapping, final boolean applyEamssFilter, final BclQualityEvaluationStrategy bclQualityEvaluationStrategy) {
        super(directory, lane, tilesToCycleFiles, outputMapping);
        this.bclQualityEvaluationStrategy = bclQualityEvaluationStrategy;
        this.applyEamssFilter = applyEamssFilter;
        this.initialize();
    }

    /**
     * Create a Bcl parser for an individual cycle and wrap it with the CycleFilesParser interface which populates
     * the correct cycle in BclData.
     *
     * @param files The files to parse.
     * @return A CycleFilesParser that populates a BclData object with data for a single cycle
     */
    @Override
    protected CycleFilesParser<BclData> makeCycleFileParser(final List<File> files) {
        return new BclDataCycleFileParser(files);
    }

    @Override
    public void initialize() {
        seekToTile(currentTile);
    }

    @Override
    public BclData next() {
        final BclData bclData = super.next();

        final byte[][] bases = bclData.getBases();
        final byte[][] qualities = bclData.getQualities();

        //first run EAMSS
        if (this.applyEamssFilter) {
            for (int i = 0; i < bases.length; i++) {
                runEamssForReadInPlace(bases[i], qualities[i]);
            }
        }

        return bclData;
    }

    /**
     * EAMSS is an Illumina Developed Algorithm for detecting reads whose quality has deteriorated towards
     * their end and revising the quality to the masking quality (2) if this is the case.  This algorithm
     * works as follows (with one exception):
     * <p>
     * Start at the end (high indices, at the right below) of the read and calculate an EAMSS tally at each
     * location as follow:
     * if(quality[i] < 15) tally += 1
     * if(quality[i] >= 15 and < 30) tally = tally
     * if(quality[i] >= 30) tally -= 2
     * <p>
     * <p>
     * For each location, keep track of this tally (e.g.)
     * Read Starts at <- this end
     * Cycle:       1  2  3  4  5  6  7  8  9
     * Bases:       A  C  T  G  G  G  T  C  A
     * Qualities:   32 32 16 15 8  10 32 2  2
     * Cycle Score: -2 -2 0  0  1  1  -2 1  1           //The EAMSS Score determined for this cycle alone
     * EAMSS TALLY: 0  0  2  2  2  1  0  2  1
     * X - Earliest instance of Max-Score
     * <p>
     * You must keep track of the maximum EAMSS tally (in this case 2) and the earliest(lowest) cycle at which
     * it occurs.  If and only if, the max EAMSS tally >= 1 then from there until the end(highest cycle) of the
     * read reassign these qualities as 2 (the masking quality).  The output qualities would therefore be
     * transformed from:
     * <p>
     * Original Qualities: 32 32 16 15 8  10 32 2  2    to
     * Final    Qualities: 32 32 2  2  2  2  2  2  2
     * X - Earliest instance of max-tally/end of masking
     * <p>
     * IMPORTANT:
     * The one exception is: If the max EAMSS Tally is preceded by a long string of G basecalls (10 or more, with a single basecall exception per10 bases)
     * then the masking continues to the beginning of that string of G's. E.g.:
     * <p>
     * Cycle:       1  2  3  4  5  6  7  8   9  10 11 12 13 14 15 16 17 18
     * Bases:       C  T  A  C  A  G  A  G   G  G  G  G  G  G  G  C  A  T
     * Qualities:   30 22 26 27 28 30 7  34  20 19 38 15 32 32 10 4  2  5
     * Cycle Score: -2  0  0  0  0 -2 1  -2  0  0  -2 0  -2 -2  1 1  1  1
     * EAMSS TALLY: -2 -5 -5 -5 -5 -5 -3 -4 -2 -2  -2 0   0  2  4 3  2  1
     * X- Earliest instance of Max-Tally
     * <p>
     * Resulting Transformation:
     * Bases:                C  T  A  C  A  G  A   G   G  G  G  G  G  G  G  C  A  T
     * Original Qualities:   30 22 26 27 28 30 7  34  20 19 38 15 32 32 10  4  2  5
     * Final    Qualities:   30 22 26 27 28  2 2   2   2  2  2  2  2  2  2  2  2  2
     * X- Earliest instance of Max-Tally
     * X - Start of EAMSS masking due to G-Run
     * <p>
     * To further clarify the exception rule here are a few examples:
     * A C G A C G G G G G G G G G G G G G G G G G G G G A C T
     * X - Earliest instance of Max-Tally
     * X - Start of EAMSS masking (with a two base call jump because we have 20 bases in the run already)
     * <p>
     * T T G G A G G G G G G G G G G G G G G G G G G A G A C T
     * X - Earliest instance of Max-Tally
     * X - We can skip this A as well as the earlier A because we have 20 or more bases in the run already
     * X - Start of EAMSS masking (with a two base call jump because we have 20 bases in the run)
     * <p>
     * T T G G G A A G G G G G G G G G G G G G G G G G G T T A T
     * X - Earliest instance of Max-Tally
     * X X - WE can skip these bases because the first A counts as the first skip and as far as the length of the string of G's is
     * concerned, these are both counted like G's
     * X - This A is the 20th base in the string of G's and therefore can be skipped
     * X - Note that the A's previous to the G's are only included because there are G's further on that are within the number
     * of allowable exceptions away (i.e. 2 in this instance), if there were NO G's after the A's you CANNOT count the A's
     * as part of the G strings (even if no exceptions have previously occured) In other words, the end of the string of G's
     * MUST end in a G not an "exception"
     * <p>
     * However, if the max-tally occurs to the right of the run of Gs then this is still part of the string of G's but does count towards
     * the number of exceptions allowable
     * (e.g.)
     * T T G G G G G G G G G G A C G
     * X - Earliest instance of Max-tally
     * The first index CAN be considered as an exception, the above would be masked to
     * the following point:
     * T T G G G G G G G G G G A C G
     * X - End of EAMSS masking due to G-Run
     * <p>
     * To sum up the final points, a string of G's CAN START with an exception but CANNOT END in an exception.
     *
     * @param bases     Bases for a single read in the cluster ( not the entire cluster )
     * @param qualities Qualities for a single read in the cluster ( not the entire cluster )
     */
    protected static void runEamssForReadInPlace(final byte[] bases, final byte[] qualities) {
        int eamssTally = 0;
        int maxTally = MIN_VALUE;
        int indexOfMax = -1;

        for (int i = bases.length - 1; i >= 0; i--) {
            final int quality = (0xff & qualities[i]);

            if (quality >= EAMSS_M2_GE_THRESHOLD) {
                eamssTally -= 2;
            } else if (quality < EAMSS_S1_LT_THRESHOLD) {
                eamssTally += 1;
            }

            if (eamssTally >= maxTally) {
                indexOfMax = i;
                maxTally = eamssTally;
            }
        }

        if (maxTally >= 1) {
            int numGs = 0;
            int exceptions = 0;

            for (int i = indexOfMax; i >= 0; i--) {
                if (bases[i] == 'G') {
                    ++numGs;
                } else {
                    final Integer skip = skipBy(i, numGs, exceptions, bases);
                    if (skip != null) {
                        exceptions += skip;
                        numGs += skip;
                        i -= (skip - 1);
                    } else {
                        break;
                    }
                }
            }

            if (numGs >= 10) {
                indexOfMax = (indexOfMax + 1) - numGs;
            }

            for (int i = indexOfMax; i < qualities.length; i++) {
                qualities[i] = MASKING_QUALITY;
            }
        }
    }

    /**
     * Determine whether or not the base at index is part of a skippable section in a run of G's, if so
     * return the number of bases that the section is composed of.
     *
     * @param index          Current index, which should be the index of a non-'G' base
     * @param numGs          The number of bases in the current string of G's for this read
     * @param prevExceptions The number of exceptions previously detected in this string by this method
     * @param bases          The bases of this read
     * @return If we have not reached our exception limit (1/every 10bases) and a G is within exceptionLimit(numGs/10)
     * indices before the current index then return index - (index of next g), else return null  Null indicates this is
     * NOT a skippable region, if we run into index 0 without finding a g then NULL is also returned
     */
    private static Integer skipBy(final int index, final int numGs, final int prevExceptions, final byte[] bases) {
        Integer skip = null;
        for (int backup = 1; backup <= index; backup++) {
            final int exceptionLimit = max((numGs + backup) / 10, 1);
            if (prevExceptions + backup > exceptionLimit) {
                break;
            }
            if (bases[index - backup] == 'G') {
                skip = backup;
                break;
            }
        }

        return skip;
    }

    private class BclDataCycleFileParser implements CycleFilesParser<BclData> {
        final CloseableIterator<BclData> reader;

        public BclDataCycleFileParser(final List<File> files) {
            reader = new BclReader(files, outputMapping.getOutputReadLengths(),
                    bclQualityEvaluationStrategy, false);
        }

        @Override
        public void close() {
            reader.close();
        }

        @Override
        public BclData next() {
            if (!hasNext()) {
                throw new NoSuchElementException();
            }
            return reader.next();
        }

        @Override
        public boolean hasNext() {
            try {
                return reader.hasNext();
            } catch (final NullPointerException npe) {
                return false;
            }
        }
    }
}
