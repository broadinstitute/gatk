package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import org.apache.commons.lang3.tuple.Pair;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.List;

public class InferOriginalReadEngineUnitTest {

    @Test
    public void testVarianceToIndelQual(){
        Assert.assertEquals(InferOriginalReadEngine.varianceToIndelQuality(30), (byte)15);
    }

    @Test
    public void testComputeVarianceAroundMostCommon(){
    }

    @Test
    public void testCigarIndex(){
        final Cigar cigar = new Cigar();
        final List<CigarElement> elements = Arrays.asList(
                new CigarElement(95, CigarOperator.M),
                new CigarElement(2, CigarOperator.I),
                new CigarElement(54, CigarOperator.M));
        for (CigarElement element : elements){
            cigar.add(element);
        }

        final int[] expectedIndexStarts = new int[]{0, 95, 97};
        final Pair<List<Integer>, List<Integer>> result = InferOriginalReadEngine.getInsertionStartIndicesAndLengths(cigar);

        Assert.assertEquals(result.getLeft().size(), 1);
        Assert.assertEquals(result.getRight().size(), 1);
        Assert.assertEquals((int) result.getLeft().get(0), 95);
        Assert.assertEquals((int) result.getRight().get(0), 2);
    }

}