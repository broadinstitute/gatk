package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collection;
import java.util.List;

public final class GraphUtilsUnitTest extends BaseTest {

    @Test
    public void compPrefixLen1(){
        final List<byte[]> bytes = Arrays.asList("ABC".getBytes(), "CDE".getBytes());
        final int pref = GraphUtils.commonMaximumPrefixLength(bytes);
        Assert.assertEquals(pref, 0);
    }

    @Test
    public void compPrefixLen2(){
        final List<byte[]> bytes = Arrays.asList("ABC".getBytes(), "ABD".getBytes());
        final int pref = GraphUtils.commonMaximumPrefixLength(bytes);
        Assert.assertEquals(pref, 2);
    }

    @Test
    public void compSuffixLen1(){
        final List<byte[]> bytes = Arrays.asList(reverse("ABC".getBytes()), reverse("CDE".getBytes()));
        final int pref = GraphUtils.commonMaximumSuffixLength(bytes, 3);
        Assert.assertEquals(pref, 0);
    }

    @Test
    public void compSuffixLen2(){
        final List<byte[]> bytes = Arrays.asList(reverse("ABC".getBytes()), reverse("ABD".getBytes()));
        final int pref = GraphUtils.commonMaximumSuffixLength(bytes, 3);
        Assert.assertEquals(pref, 2);
    }

    @Test
    public void compSuffixLen3(){
        final List<byte[]> bytes = Arrays.asList(reverse("ABC".getBytes()), reverse("ABD".getBytes()));
        final int pref = GraphUtils.commonMaximumSuffixLength(bytes, 1);
        Assert.assertEquals(pref, 1);
    }

    @Test
    public void kMers(){
        final SeqVertex v1 = new SeqVertex("fred");
        final SeqVertex v2 = new SeqVertex("frodo");
        final Collection<SeqVertex> vertices = Arrays.asList(v1, v2);
        final List<byte[]> kmers = GraphUtils.getKmers(vertices);
        Assert.assertEquals(kmers.size(), 2);
        Assert.assertTrue(Arrays.equals(kmers.get(0), "fred".getBytes()));
        Assert.assertTrue(Arrays.equals(kmers.get(1), "frodo".getBytes()));
    }

    private static byte[] reverse(final byte[] arr){
        final byte[] result = new byte[arr.length];
        for (int i = 0; i < arr.length; i++) {
            result[arr.length - 1 - i] = arr[i];
        }
        return result;
    }
}
