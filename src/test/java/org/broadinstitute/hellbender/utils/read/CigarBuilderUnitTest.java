package org.broadinstitute.hellbender.utils.read;

import htsjdk.samtools.TextCigarCodec;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

public class CigarBuilderUnitTest {

    @DataProvider(name = "simple_concatenation")
    public Iterator<Object[]> cigarAlgebra() {
        final List<List<String>> leadingClips = Arrays.asList(
                Collections.emptyList(),
                Collections.singletonList("10H"),
                Collections.singletonList("10S"),
                Arrays.asList("10H","10S")
        );

        final List<List<String>> middleOperators = Arrays.asList(
                Collections.singletonList("10M"),
                Arrays.asList("10M","10I", "10M"),
                Arrays.asList("10M","10D", "10M")
        );

        final List<List<String>> trailingClips = Arrays.asList(
                Collections.emptyList(),
                Collections.singletonList("10H"),
                Collections.singletonList("10S"),
                Arrays.asList("10S","10H")
        );

        final List<Object[]> result = new ArrayList<>();
        for (final List<String> leadingString : leadingClips) {
            for (final List<String> middleString : middleOperators) {
                for (final List<String> trailingString : trailingClips) {
                    final List<String> allElements = new ArrayList<>();
                    allElements.addAll(leadingString);
                    allElements.addAll(middleString);
                    allElements.addAll(trailingString);
                    result.add(new Object[] {allElements});
                }
            }
        }

        return result.iterator();
    }

    @Test(dataProvider = "simple_concatenation")
    public void testSimpleConcatenation(final List<String> cigarElementStrings) {
        final CigarBuilder builder = new CigarBuilder();
        for (final String elementString : cigarElementStrings) {
            builder.add(TextCigarCodec.decode(elementString).getFirstCigarElement());
        }

        final String expected = String.join("", cigarElementStrings);
        Assert.assertEquals(builder.make().toString(), expected);
    }

    @DataProvider(name = "initial_and_final_deletions")
    public Object[][] initialAndFinalDeletions() {
        return new Object[][] {
                {Arrays.asList("10M", "10D"), "10M"},
                {Arrays.asList("10D", "10M"), "10M"},
                {Arrays.asList("10H", "10D", "10M"), "10H10M"},
                {Arrays.asList("10S", "10D", "10M"), "10S10M"},
                {Arrays.asList("10S", "10D", "10M", "10S"), "10S10M10S"},
                {Arrays.asList("10M", "10D", "10S"), "10M10S"},
                {Arrays.asList("10M", "10D", "10H"), "10M10H"},
                {Arrays.asList("10S", "10M", "10D", "10H"), "10S10M10H"}
        };
    }

    @Test(dataProvider = "initial_and_final_deletions")
    public void testInitialAndFinalDeletions(final List<String> cigarElementStrings, final String expected) {
        final CigarBuilder builder = new CigarBuilder();
        for (final String elementString : cigarElementStrings) {
            builder.add(TextCigarCodec.decode(elementString).getFirstCigarElement());
        }

        Assert.assertEquals(builder.make().toString(), expected);
    }

    @DataProvider(name = "retain_deletions")
    public Object[][] retainDeletions() {
        return new Object[][] {
                // terminal deletions should be kept
                {Arrays.asList("10M", "10D"), "10M10D"},
                {Arrays.asList("10D", "10M"), "10D10M"},
                {Arrays.asList("10H", "10D", "10M"), "10H10D10M"},
                {Arrays.asList("10S", "10D", "10M"), "10S10D10M"},
                {Arrays.asList("10S", "10D", "10M", "10S"), "10S10D10M10S"},
                {Arrays.asList("10M", "10D", "10S"), "10M10D10S"},
                {Arrays.asList("10M", "10D", "10H"), "10M10D10H"},
                {Arrays.asList("10S", "10M", "10D", "10H"), "10S10M10D10H"},

                // merging consecutive elements should still work
                {Arrays.asList("10M", "10D", "10D"), "10M20D"},
                {Arrays.asList("10M", "10M", "10D", "10D"), "20M20D"}
        };
    }

    @Test(dataProvider = "retain_deletions")
    public void testRetainDeletions(final List<String> cigarElementStrings, final String expected) {
        final CigarBuilder builder = new CigarBuilder(false);
        for (final String elementString : cigarElementStrings) {
            builder.add(TextCigarCodec.decode(elementString).getFirstCigarElement());
        }

        Assert.assertEquals(builder.make().toString(), expected);
    }

    @DataProvider(name = "merge_consecutive")
    public Object[][] mergeConsecutive() {
        return new Object[][] {
                {Arrays.asList("10H", "10H", "10M"), "20H10M"},
                {Arrays.asList("10S", "10M", "10M"), "10S20M"},
                {Arrays.asList("10S", "10M", "10S", "10S"), "10S10M20S"},
                {Arrays.asList("10S", "10M", "10I", "10I", "10I", "10S", "10H"), "10S10M30I10S10H"},
                {Arrays.asList("10S", "10S", "10M", "10M", "10I", "10I", "10S", "10H"), "20S20M20I10S10H"},
        };
    }

    @Test(dataProvider = "merge_consecutive")
    public void testMergeConsecutive(final List<String> cigarElementStrings, final String expected) {
        final CigarBuilder builder = new CigarBuilder();
        for (final String elementString : cigarElementStrings) {
            builder.add(TextCigarCodec.decode(elementString).getFirstCigarElement());
        }

        Assert.assertEquals(builder.make().toString(), expected);
    }

    @DataProvider(name = "tricky")
    public Object[][] tricky() {
        return new Object[][] {
                {Arrays.asList("10H", "10H", "10D", "10D", "10M"), "20H10M"},
        };
    }

    @Test(dataProvider = "tricky")
    public void testTrickyCases(final List<String> cigarElementStrings, final String expected) {
        final CigarBuilder builder = new CigarBuilder();
        for (final String elementString : cigarElementStrings) {
            builder.add(TextCigarCodec.decode(elementString).getFirstCigarElement());
        }

        Assert.assertEquals(builder.make().toString(), expected);
    }

    @DataProvider(name = "indel_sandwich")
    public Object[][] indelSandwich() {
        return new Object[][] {
                {Arrays.asList("10M", "10I", "10D", "10M"), "10M10D10I10M"},
                {Arrays.asList("10M", "10D", "10I", "10M"), "10M10D10I10M"},
                {Arrays.asList("10M", "10I", "10D", "10I", "10M"), "10M10D20I10M"},
                {Arrays.asList("10M", "10I", "10D", "10I", "10D", "10I", "10M"), "10M20D30I10M"},
                {Arrays.asList("10M", "10I", "10D", "10I", "10M", "10D", "10I", "10M"), "10M10D20I10M10D10I10M"},

                //does the indel sandwich logic interfere with removing leading/trailing deletions
                {Arrays.asList("10D", "10I", "10M"), "10I10M"},
                {Arrays.asList("10M", "10I", "10D"), "10M10I"},
                {Arrays.asList("10M", "10D", "10I"), "10M10I"},
                {Arrays.asList("10M", "10D", "10I", "10S"), "10M10I10S"},
                {Arrays.asList("10S", "10D", "10I", "10M"), "10S10I10M"},
                {Arrays.asList("10S", "10I", "10D", "10I", "10M"), "10S20I10M"},
        };
    }

    @Test(dataProvider = "indel_sandwich")
    public void testIndelSandwich(final List<String> cigarElementStrings, final String expected) {
        final CigarBuilder builder = new CigarBuilder();
        for (final String elementString : cigarElementStrings) {
            builder.add(TextCigarCodec.decode(elementString).getFirstCigarElement());
        }

        Assert.assertEquals(builder.make().toString(), expected);
    }

    @DataProvider(name = "invalid")
    public Object[][] invalid() {
        return new Object[][] {
                // completely soft-clipped
                {Arrays.asList("10S")},
                {Arrays.asList("10S", "10S")},

                // also completely clipped
                {Arrays.asList("10S", "10D")},
                {Arrays.asList("10S", "10D", "10S")},
                {Arrays.asList("10S", "10D", "10D", "10S")},

                // wrong order of hard and soft clips
                {Arrays.asList("10S", "10H", "10M")},
                {Arrays.asList("10M", "10H", "10S")},

                // clipping in middle of read
                {Arrays.asList("10M", "10H", "10M")},
                {Arrays.asList("10M", "10S", "10M")},
        };
    }

    @Test(dataProvider = "invalid", expectedExceptions = IllegalStateException.class)
    public void testInvalid(final List<String> cigarElementStrings) {
        final CigarBuilder builder = new CigarBuilder();
        for (final String elementString : cigarElementStrings) {
            builder.add(TextCigarCodec.decode(elementString).getFirstCigarElement());
        }

        builder.make();
    }

    @DataProvider(name = "removed_deletions")
    public Object[][] removedDeletions() {
        return new Object[][] {
                {Arrays.asList("10M"), 0, 0},
                {Arrays.asList("10S", "10M"), 0, 0},
                {Arrays.asList("10M","10S"), 0, 0},
                {Arrays.asList("10M", "10I", "10D", "10M"), 0, 0},
                {Arrays.asList("10M", "10D", "10I", "10M"), 0, 0},

                {Arrays.asList("10D", "10I", "10M"), 10, 0},
                {Arrays.asList("10D", "10D", "10I", "10M"), 20, 0},
                {Arrays.asList("10D", "10D", "10I", "10D", "10M"), 30, 0},
                {Arrays.asList("10S", "10D", "10D", "10I", "10D", "10M"), 30, 0},

                {Arrays.asList("10M", "10I", "10D"), 0, 10},
                {Arrays.asList("10M", "10D", "10I"), 0, 10},
                {Arrays.asList("10M", "10D", "10I", "10D"), 0, 20},
                {Arrays.asList("10M", "10D", "10I", "10D", "10S", "10H"), 0, 20},

                {Arrays.asList("10H", "10S", "10D", "10M", "10D", "10I", "10D", "10S", "10H"), 10, 20},
        };
    }

    @Test(dataProvider = "removed_deletions")
    public void testRemovedDeletions(final List<String> cigarElementStrings, final int removedLeading, final int removedTrailing) {
        final CigarBuilder builder = new CigarBuilder();
        for (final String elementString : cigarElementStrings) {
            builder.add(TextCigarCodec.decode(elementString).getFirstCigarElement());
        }

        builder.make();
        Assert.assertEquals(builder.getLeadingDeletionBasesRemoved(), removedLeading);
        Assert.assertEquals(builder.getTrailingDeletionBasesRemoved(), removedTrailing);
    }

    @DataProvider(name = "removed_deletions_two_makes")
    public Object[][] removedDeletionsTwoMakes() {
        return new Object[][] {
                {Arrays.asList("10M"), Arrays.asList("10M"), 0, 0},
                {Arrays.asList("10M", "10I"), Arrays.asList("10D", "10M"), 0, 0},
                {Arrays.asList("10M", "10D"), Arrays.asList("10I", "10M"), 0, 0},

                {Arrays.asList("10D", "10I"), Arrays.asList("10M"), 10, 0},
                {Arrays.asList("10D", "10D", "10I"), Arrays.asList("10D", "10M"), 30, 0},

                {Arrays.asList("10H", "10S", "10D", "10M"), Arrays.asList("10D", "10I", "10D", "10S", "10H"), 10, 20},
        };
    }

    @Test(dataProvider = "removed_deletions_two_makes")
    public void testRemovedDeletionsWithTwoMakes(final List<String> cigarElementStrings1, final List<String> cigarElementStrings2, final int removedLeading, final int removedTrailing) {
        final CigarBuilder builder = new CigarBuilder();
        for (final String elementString : cigarElementStrings1) {
            builder.add(TextCigarCodec.decode(elementString).getFirstCigarElement());
        }

        builder.make();

        for (final String elementString : cigarElementStrings2) {
            builder.add(TextCigarCodec.decode(elementString).getFirstCigarElement());
        }

        builder.make();


        Assert.assertEquals(builder.getLeadingDeletionBasesRemoved(), removedLeading);
        Assert.assertEquals(builder.getTrailingDeletionBasesRemoved(), removedTrailing);
    }
}