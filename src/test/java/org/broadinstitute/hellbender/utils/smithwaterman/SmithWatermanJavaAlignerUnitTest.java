package org.broadinstitute.hellbender.utils.smithwaterman;

import org.broadinstitute.gatk.nativebindings.smithwaterman.SWOverhangStrategy;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWParameters;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.List;
import java.util.Random;

public final class SmithWatermanJavaAlignerUnitTest extends SmithWatermanAlignerAbstractUnitTest {

    @Override
    protected SmithWatermanJavaAligner getAligner() {
        return SmithWatermanJavaAligner.getInstance();
    }

    @Test
    public void testShortcutsAgreeWithFullSmithWaterman() {
        // reads from anywhere in short, often repetitive haplotypes, including overhanging their ends, with a few
        // substitutions, insertions and deletions, so that exact and one-mismatch matches, ties between them, and nearly
        // exact matches with overhangs all occur.  Some parameters allow the one-mismatch shortcut and some don't.
        final List<SWParameters> parametersList = Arrays.asList(SmithWatermanAlignmentConstants.ALIGNMENT_TO_BEST_HAPLOTYPE_SW_PARAMETERS,
                SmithWatermanAlignmentConstants.STANDARD_NGS, SmithWatermanAlignmentConstants.NEW_SW_PARAMETERS,
                SmithWatermanAlignmentConstants.ORIGINAL_DEFAULT);
        final Random random = new Random(13);
        int numOneMismatchShortcuts = 0;
        for (int trial = 0; trial < 20000; trial++) {
            final String alphabet = trial % 3 == 0 ? "AC" : "ACGT";
            final byte[] haplotype = randomBases(random, 8 + random.nextInt(60), alphabet);
            final byte[] read = mutate(random, readFrom(random, haplotype, alphabet), alphabet);
            for (final SWParameters parameters : parametersList) {
                for (final SWOverhangStrategy strategy : Arrays.asList(SWOverhangStrategy.SOFTCLIP, SWOverhangStrategy.IGNORE)) {
                    final SmithWatermanAlignment expected = SmithWatermanJavaAligner.alignWithoutShortcuts(haplotype, read, parameters, strategy);
                    final SmithWatermanAlignment actual = getAligner().align(haplotype, read, parameters, strategy);
                    final String message = new String(haplotype) + " " + new String(read) + " " + parameters + " " + strategy;
                    Assert.assertEquals(actual.getCigar(), expected.getCigar(), message);
                    Assert.assertEquals(actual.getAlignmentOffset(), expected.getAlignmentOffset(), message);
                }
                if (SmithWatermanJavaAligner.lastIndexOfOneMismatchAlignment(haplotype, read, parameters) != -1) {
                    numOneMismatchShortcuts++;
                }
            }
        }
        // make sure the shortcut was exercised
        Assert.assertTrue(numOneMismatchShortcuts > 1000, "only " + numOneMismatchShortcuts + " one-mismatch shortcuts");
    }

    private static byte[] randomBases(final Random random, final int length, final String alphabet) {
        final byte[] bases = new byte[length];
        for (int n = 0; n < length; n++) {
            bases[n] = (byte) alphabet.charAt(random.nextInt(alphabet.length()));
        }
        return bases;
    }

    // a stretch of the haplotype, sometimes overhanging one of its ends by a few random bases, and sometimes longer
    private static byte[] readFrom(final Random random, final byte[] haplotype, final String alphabet) {
        final int length = 3 + random.nextInt(haplotype.length + 3);
        final int start = random.nextInt(haplotype.length + 7) - 3;
        final byte[] read = new byte[length];
        for (int n = 0; n < length; n++) {
            final int position = start + n;
            read[n] = position >= 0 && position < haplotype.length ? haplotype[position] :
                    (byte) alphabet.charAt(random.nextInt(alphabet.length()));
        }
        return read;
    }

    // up to three edits, mostly substitutions
    private static byte[] mutate(final Random random, final byte[] read, final String alphabet) {
        byte[] result = read;
        final int numEdits = random.nextInt(4);
        for (int edit = 0; edit < numEdits && result.length > 1; edit++) {
            final int position = random.nextInt(result.length);
            final byte base = (byte) alphabet.charAt(random.nextInt(alphabet.length()));
            final int type = random.nextInt(5);
            if (type < 3) {
                result[position] = base;
            } else if (type == 3) {
                final byte[] inserted = new byte[result.length + 1];
                System.arraycopy(result, 0, inserted, 0, position);
                inserted[position] = base;
                System.arraycopy(result, position, inserted, position + 1, result.length - position);
                result = inserted;
            } else {
                final byte[] deleted = new byte[result.length - 1];
                System.arraycopy(result, 0, deleted, 0, position);
                System.arraycopy(result, position + 1, deleted, position, result.length - position - 1);
                result = deleted;
            }
        }
        return result;
    }
}
