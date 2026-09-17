package org.broadinstitute.hellbender.utils.smithwaterman;

import org.apache.commons.lang3.ArrayUtils;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWOverhangStrategy;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWParameters;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.RandomSequenceTestUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.testng.Assert;
import org.testng.SkipException;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Random;

/**
 * Runs the shared aligner tests against the native Smith-Waterman aligner, and checks that it produces exactly the
 * alignments the Java aligner does on many random sequence pairs.
 */
public class SmithWatermanIntelAlignerUnitTest extends SmithWatermanAlignerAbstractUnitTest {

    private static final int PAIRS_PER_SCENARIO = 100;

    /** Every parameter set GATK aligns with. */
    private static final SWParameters[] PARAMETER_SETS = {
            SmithWatermanAlignmentConstants.ORIGINAL_DEFAULT,
            SmithWatermanAlignmentConstants.STANDARD_NGS,
            SmithWatermanAlignmentConstants.NEW_SW_PARAMETERS,
            SmithWatermanAlignmentConstants.ALIGNMENT_TO_BEST_HAPLOTYPE_SW_PARAMETERS,
    };

    /**
     * Every overhang strategy GATK aligns with. IGNORE is left out: nothing in GATK aligns with it, and the native
     * aligner soft-clips a trailing overhang under IGNORE where the Java aligner reports it as matched bases.
     */
    private static final SWOverhangStrategy[] OVERHANG_STRATEGIES = {
            SWOverhangStrategy.SOFTCLIP,
            SWOverhangStrategy.INDEL,
            SWOverhangStrategy.LEADING_INDEL,
    };

    /** Kinds of random sequence pair, each exercising a different part of the alignment and overhang handling. */
    private enum Scenario {
        /** A reference of 30 to 400 bases and an alternate that is a slice of it with a few substitutions and indels. */
        EDITED_SLICE {
            @Override
            byte[][] pair(final Random rng) {
                final byte[] reference = RandomSequenceTestUtils.randomBases(rng, 30 + rng.nextInt(371));
                return new byte[][]{reference, editedSliceOf(rng, reference)};
            }
        },
        /** An edited slice with up to 30 random bases hanging off each end. */
        EDITED_SLICE_WITH_OVERHANGS {
            @Override
            byte[][] pair(final Random rng) {
                final byte[][] pair = EDITED_SLICE.pair(rng);
                final byte[] alternate = ArrayUtils.addAll(ArrayUtils.addAll(RandomSequenceTestUtils.randomBases(rng, rng.nextInt(31)), pair[1]),
                                                           RandomSequenceTestUtils.randomBases(rng, rng.nextInt(31)));
                return new byte[][]{pair[0], alternate};
            }
        },
        /** An alternate of 30 to 400 bases and a reference that is an edited slice of it. */
        ALTERNATE_LONGER_THAN_REFERENCE {
            @Override
            byte[][] pair(final Random rng) {
                final byte[] alternate = RandomSequenceTestUtils.randomBases(rng, 30 + rng.nextInt(371));
                return new byte[][]{editedSliceOf(rng, alternate), alternate};
            }
        },
        /** Tandem repeats of a one to three base unit, where many alignments score the same. */
        REPETITIVE {
            @Override
            byte[][] pair(final Random rng) {
                final byte[] unit = RandomSequenceTestUtils.randomBases(rng, 1 + rng.nextInt(3));
                final byte[] reference = new byte[30 + rng.nextInt(171)];
                for (int i = 0; i < reference.length; i++) {
                    reference[i] = unit[i % unit.length];
                }
                return new byte[][]{reference, editedSliceOf(rng, reference)};
            }
        },
        /** Two independent random sequences. */
        UNRELATED {
            @Override
            byte[][] pair(final Random rng) {
                return new byte[][]{RandomSequenceTestUtils.randomBases(rng, 30 + rng.nextInt(371)), RandomSequenceTestUtils.randomBases(rng, 10 + rng.nextInt(291))};
            }
        };

        abstract byte[][] pair(final Random rng);

        /** A slice of at least ten bases of {@code sequence} with up to three edits, each indel at most eight bases. */
        static byte[] editedSliceOf(final Random rng, final byte[] sequence) {
            final int length = 10 + rng.nextInt(Math.min(sequence.length, 300) - 9);
            final int start = rng.nextInt(sequence.length - length + 1);
            return RandomSequenceTestUtils.withRandomEdits(rng, Arrays.copyOfRange(sequence, start, start + length), rng.nextInt(4), 8);
        }
    }

    @Override
    protected SmithWatermanIntelAligner getAligner() {
        try {
            return new SmithWatermanIntelAligner();
        } catch (final UserException.HardwareFeatureException e) {
            throw new SkipException("AVX SmithWaterman is not supported on this system or the library is not available");
        }
    }

    @DataProvider
    public Object[][] scenarios() {
        return new Object[][]{
                {Scenario.EDITED_SLICE},
                {Scenario.EDITED_SLICE_WITH_OVERHANGS},
                {Scenario.ALTERNATE_LONGER_THAN_REFERENCE},
                {Scenario.REPETITIVE},
                {Scenario.UNRELATED},
        };
    }

    @Test(dataProvider = "scenarios")
    public void testNativeAlignerMatchesJavaAligner(final Scenario scenario) {
        Utils.resetRandomGenerator();
        final Random rng = Utils.getRandomGenerator();
        final SmithWatermanAligner javaAligner = SmithWatermanJavaAligner.getInstance();
        try (final SmithWatermanAligner nativeAligner = getAligner()) {
            for (int i = 0; i < PAIRS_PER_SCENARIO; i++) {
                final byte[][] pair = scenario.pair(rng);
                for (final SWParameters parameters : PARAMETER_SETS) {
                    for (final SWOverhangStrategy strategy : OVERHANG_STRATEGIES) {
                        final SmithWatermanAlignment expected = javaAligner.align(pair[0], pair[1], parameters, strategy);
                        final SmithWatermanAlignment actual = nativeAligner.align(pair[0], pair[1], parameters, strategy);
                        final String description = String.format("%s: reference %s, alternate %s, match %d mismatch %d gap open %d gap extend %d, %s",
                                scenario, new String(pair[0]), new String(pair[1]), parameters.getMatchValue(), parameters.getMismatchPenalty(),
                                parameters.getGapOpenPenalty(), parameters.getGapExtendPenalty(), strategy);
                        Assert.assertEquals(actual.getAlignmentOffset(), expected.getAlignmentOffset(), "alignment offset for " + description);
                        Assert.assertEquals(actual.getCigar().toString(), expected.getCigar().toString(), "cigar for " + description);
                    }
                }
            }
        }
    }
}
