package org.broadinstitute.hellbender.tools.spark.sv;

import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.junit.Assert;
import org.testng.annotations.Test;

/**
 * Created by shuang on 9/16/16.
 */
public final class LocalAssemblyContigUnitTest extends BaseTest {

    @Test
    public void testEquals(){

        final LocalAssemblyContig contig = new LocalAssemblyContig(1, "contig-1", "CTTGCGGTGGGAGAGGGGAACCTGCGCCCGCGATTCTCTCTCTCTCCCCCTCTCTCCCTCTCTCTTTCTCGCCCCCCCCCTCCCCCTCTCTCTCCCTCCCTCCCTCTCTCTTTCTCTCTCCCCCTCCCTCCCTCTCTCTCCCTCCCTCCCT");
        final LocalAssemblyContig identicalContig = new LocalAssemblyContig(1, "contig-1", "CTTGCGGTGGGAGAGGGGAACCTGCGCCCGCGATTCTCTCTCTCTCCCCCTCTCTCCCTCTCTCTTTCTCGCCCCCCCCCTCCCCCTCTCTCTCCCTCCCTCCCTCTCTCTTTCTCTCTCCCCCTCCCTCCCTCTCTCTCCCTCCCTCCCT", null);

        Assert.assertEquals(contig, identicalContig);
    }

    @Test
    public void testHashCode(){
        final LocalAssemblyContig contig = new LocalAssemblyContig(1, "contig-1", "CTTGCGGTGGGAGAGGGGAACCTGCGCCCGCGATTCTCTCTCTCTCCCCCTCTCTCCCTCTCTCTTTCTCGCCCCCCCCCTCCCCCTCTCTCTCCCTCCCTCCCTCTCTCTTTCTCTCTCCCCCTCCCTCCCTCTCTCTCCCTCCCTCCCT");
        final LocalAssemblyContig identicalContig = new LocalAssemblyContig(1, "contig-1", "CTTGCGGTGGGAGAGGGGAACCTGCGCCCGCGATTCTCTCTCTCTCCCCCTCTCTCCCTCTCTCTTTCTCGCCCCCCCCCTCCCCCTCTCTCTCCCTCCCTCCCTCTCTCTTTCTCTCTCCCCCTCCCTCCCTCTCTCTCCCTCCCTCCCT", null);

        Assert.assertEquals(contig.hashCode(), identicalContig.hashCode());
    }
}
