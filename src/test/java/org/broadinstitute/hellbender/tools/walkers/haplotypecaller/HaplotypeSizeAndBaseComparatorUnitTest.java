package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public final class HaplotypeSizeAndBaseComparatorUnitTest extends BaseTest {
    @Test
    public void testComparison() {
        // desired ordering is by size first, subordered by lexacographic relationship between bases
        final List<String> rawStrings = Arrays.asList("A", "C", "AC", "CC", "CT", "AAT", "ACT", "GAT", "ACGT");
        final List<String> lexStrings = new ArrayList<>(rawStrings);

        for ( final List<String> seqs : Utils.makePermutations(lexStrings, lexStrings.size(), false) ) {
            final List<Haplotype> haps = new ArrayList<>(seqs.size());
            for ( final String seq : seqs ) {
                haps.add(new Haplotype(seq.getBytes(), false));
            }

            Collections.sort(haps, Haplotype.SIZE_AND_BASE_ORDER);
            for ( int i = 0; i < lexStrings.size(); i++ )
                Assert.assertEquals(haps.get(i).getBaseString(), lexStrings.get(i), "Failed sort " + haps + " expected " + lexStrings);
        }
    }
}