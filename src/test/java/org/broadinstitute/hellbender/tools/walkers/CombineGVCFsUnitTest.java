package org.broadinstitute.hellbender.tools.walkers;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.testng.Assert;
import org.testng.annotations.Test;
import org.mockito.Mockito;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static org.mockito.Mockito.when;

/**
 * Created by emeryj on 8/25/17.
 */
public class CombineGVCFsUnitTest extends BaseTest{
    public static final Allele REF = Allele.create("A", true);
    public static final Allele ALT = Allele.create("C", false);

//    // This test is disabled until mocking classes with mockito is allowed in GATK4.
//    @Test(enabled = false )
//    public void testUpdatePostiionalState() {
//        ReferenceContext referenceContext1 = Mockito.mock(ReferenceContext.class);
//        when(referenceContext1.getWindow()).thenReturn(new SimpleInterval("20",1000,1100));
//        when(referenceContext1.getBases()).thenReturn(new byte[]{'A'});
//
//        ReferenceContext referenceContext2 = Mockito.mock(ReferenceContext.class);
//        when(referenceContext2.getWindow()).thenReturn(new SimpleInterval("20",1000,1050));
//        when(referenceContext2.getBases()).thenReturn(new byte[]{'C'});
//
//        ReferenceContext referenceContext3 = Mockito.mock(ReferenceContext.class);
//        when(referenceContext3.getWindow()).thenReturn(new SimpleInterval("20",1050,1100));
//        when(referenceContext3.getBases()).thenReturn(new byte[]{'G'});
//
//        ReferenceContext referenceContext4 = Mockito.mock(ReferenceContext.class);
//        when(referenceContext4.getWindow()).thenReturn(new SimpleInterval("20",1050,1150));
//        when(referenceContext4.getBases()).thenReturn(new byte[]{'T'});
//
//        List<VariantContext> VCs = new ArrayList<>();
//        VariantContext A = new VariantContextBuilder().loc("20",1001,1099).alleles(Arrays.asList(REF, GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE)).attribute(VCFConstants.END_KEY, 1099).make();
//
//        VCs.add(A);
//        CombineGVCFs.QueuedContextState contextState = new CombineGVCFs.QueuedContextState(VCs, referenceContext1, new SimpleInterval("20",1001,1099));
//        Assert.assertEquals(contextState.refBasesLoc.getEnd(),1100);
//        Assert.assertEquals(contextState.refBases[0],'A');
//
//        VCs.add(A);
//        contextState = CombineGVCFs.updatePositionalState(contextState, Arrays.asList(A,A) ,referenceContext1);
//        Assert.assertEquals(contextState.refBasesLoc.getEnd(),1100);
//        Assert.assertEquals(contextState.refBases[0],'A');
//        Assert.assertEquals(contextState.variantContexts.size(), 2);
//
//        VCs.clear();
//        VCs.add(A);
//        contextState = CombineGVCFs.updatePositionalState(contextState, Arrays.asList(A) ,referenceContext4);
//        Assert.assertEquals(contextState.refBasesLoc.getEnd(),1150);
//        Assert.assertEquals(contextState.refBases[0],'T');
//        Assert.assertEquals(contextState.variantContexts.size(), 1);
//
//        VCs.add(A);
//        contextState = CombineGVCFs.updatePositionalState(contextState, Arrays.asList(A, A) ,referenceContext3);
//        Assert.assertEquals(contextState.refBasesLoc.getEnd(),1150);
//        Assert.assertEquals(contextState.refBases[0],'T');
//        Assert.assertEquals(contextState.variantContexts.size(), 2);
//    }

}