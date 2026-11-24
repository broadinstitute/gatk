package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.positional;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAlignmentConstants;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.List;

public class PositionalEventUnitTest extends GATKBaseTest {

    @DataProvider(name = "normalize")
    public Object[][] normalizeData() {
        return new Object[][] {
                // reference, position, ref, alt, floor, expected events
                {"ACGTA", 0, "ACGTA", "ACCTA", 0, List.of(new PositionalEvent(2, "G", "C"))},
                {"ACGTA", 0, "ACGTA", "ACGTA", 0, List.of()},
                // two adjacent mismatches are cheaper than an insertion and a deletion
                {"ACGTA", 0, "ACGT", "AGCT", 0, List.of(new PositionalEvent(1, "C", "G"), new PositionalEvent(2, "G", "C"))},
                // a deletion in a homopolymer is left-aligned ...
                {"GAAAAT", 3, "A", "", 0, List.of(new PositionalEvent(1, "A", ""))},
                // ... but not past the floor
                {"GAAAAT", 3, "A", "", 2, List.of(new PositionalEvent(2, "A", ""))},
                // an insertion in a dinucleotide repeat is left-aligned, rotating its bases
                {"TCACAG", 5, "", "CA", 0, List.of(new PositionalEvent(1, "", "CA"))},
                // the same insertion described with flanking reference bases
                {"TCACAG", 4, "AG", "ACAG", 0, List.of(new PositionalEvent(1, "", "CA"))},
                // two deletions around a matching base
                {"AAACGTTTGGG", 2, "ACGTTTG", "AGTTG", 0, List.of(new PositionalEvent(3, "C", ""), new PositionalEvent(5, "T", ""))},
        };
    }

    @Test(dataProvider = "normalize")
    public void testNormalize(final String reference, final int position, final String ref, final String alt, final int floor,
                              final List<PositionalEvent> expected) {
        Assert.assertEquals(PositionalEvent.normalize(reference, position, ref, alt, floor, SmithWatermanAlignmentConstants.NEW_SW_PARAMETERS), expected);
    }

    @Test
    public void testLeftAligned() {
        // a deletion that a floor stopped short of the start of a homopolymer
        Assert.assertEquals(new PositionalEvent(3, "A", "").leftAligned("GAAAAT"), new PositionalEvent(1, "A", ""));
        // an insertion at the end of a dinucleotide repeat, rotating its bases
        Assert.assertEquals(new PositionalEvent(5, "", "CA").leftAligned("TCACAG"), new PositionalEvent(1, "", "CA"));
        // already leftmost, and a substitution
        Assert.assertEquals(new PositionalEvent(1, "A", "").leftAligned("GAAAAT"), new PositionalEvent(1, "A", ""));
        Assert.assertEquals(new PositionalEvent(2, "A", "C").leftAligned("GAAAAT"), new PositionalEvent(2, "A", "C"));
    }

    @Test
    public void testApply() {
        final List<PositionalEvent> events = List.of(new PositionalEvent(12, "", "TT"), new PositionalEvent(11, "C", "G"));
        Assert.assertEquals(PositionalEvent.apply("ACGT", 10, events), "AGTTGT");
        Assert.assertEquals(PositionalEvent.apply("ACGT", 10, List.of(new PositionalEvent(13, "T", ""))), "ACG");
        Assert.assertNull(PositionalEvent.apply("ACGT", 10, List.of(new PositionalEvent(11, "CG", ""), new PositionalEvent(12, "G", "A"))));
        Assert.assertNull(PositionalEvent.apply("ACGT", 10, List.of(new PositionalEvent(13, "TA", ""))));
    }
}
