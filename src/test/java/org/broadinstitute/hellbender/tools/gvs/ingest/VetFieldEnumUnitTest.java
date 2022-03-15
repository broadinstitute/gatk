package org.broadinstitute.hellbender.tools.gvs.ingest;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.testng.annotations.Test;

import java.util.Arrays;
import org.testng.Assert;

import static htsjdk.variant.vcf.VCFConstants.DEPTH_KEY;
import static org.broadinstitute.hellbender.utils.variant.GATKVCFConstants.AS_RAW_MAP_QUAL_RANK_SUM_KEY;
import static org.broadinstitute.hellbender.utils.variant.GATKVCFConstants.AS_RAW_QUAL_APPROX_KEY;
import static org.broadinstitute.hellbender.utils.variant.GATKVCFConstants.AS_RAW_READ_POS_RANK_SUM_KEY;
import static org.broadinstitute.hellbender.utils.variant.GATKVCFConstants.AS_RAW_RMS_MAPPING_QUALITY_KEY;
import static org.broadinstitute.hellbender.utils.variant.GATKVCFConstants.AS_SB_TABLE_KEY;
import static org.broadinstitute.hellbender.utils.variant.GATKVCFConstants.AS_VARIANT_DEPTH_KEY;
import static org.broadinstitute.hellbender.utils.variant.GATKVCFConstants.MAP_QUAL_RANK_SUM_KEY;
import static org.broadinstitute.hellbender.utils.variant.GATKVCFConstants.RAW_MAPPING_QUALITY_WITH_DEPTH_KEY;
import static org.broadinstitute.hellbender.utils.variant.GATKVCFConstants.RAW_QUAL_APPROX_KEY;
import static org.broadinstitute.hellbender.utils.variant.GATKVCFConstants.READ_POS_RANK_SUM_KEY;
import static org.broadinstitute.hellbender.utils.variant.GATKVCFConstants.STRAND_BIAS_BY_SAMPLE_KEY;
import static org.broadinstitute.hellbender.utils.variant.GATKVCFConstants.VARIANT_DEPTH_KEY;

public class VetFieldEnumUnitTest {
    private static final String SAMPLE_1 = "NA1";

    @Test
    public void testAlleleSpecificOneAlt() {
        // Original Record:
        //
        // chr1	10329	.	C	A,<NON_REF>	43	.	AS_QUALapprox=|74|0;AS_RAW_BaseQRankSum=|-0.2,1|1.8,1;AS_RAW_MQ=29707.00|39366.00|2405.00;AS_RAW_MQRankSum=|-0.2,1|-2.5,1;AS_RAW_ReadPosRankSum=|2.4,1|1.5,1;AS_SB_TABLE=1,21|3,39|3,11;AS_VarDP=22|42|0;BaseQRankSum=0.452;DP=95;MQRankSum=-0.951;QUALapprox=74;RAW_GT_COUNT=0,1,0;RAW_MQandDP=83487,95;ReadPosRankSum=2.470;VarDP=64	GT:AD:DP:GP:GQ:PG:PL:SB	0/1:22,42,0:64:43,0,37.01,706,420,469.01:36:0,31,34.01,30,61,33.01:74,0,34,707,390,467:1,21,6,50
        //
        VariantContextBuilder builderA =
                new VariantContextBuilder("a","1",10329,10329,
                        Arrays.asList(Allele.REF_C,Allele.ALT_A,Allele.NON_REF_ALLELE));


        Genotype g = new GenotypeBuilder(SAMPLE_1)
                .alleles(Arrays.asList(Allele.REF_C, Allele.ALT_A))
                .PL(new int[]{74,0,34,707,390,467})
                .DP(64)
                .GQ(36)
                .AD(new int[]{22,42,0})
                .attribute(STRAND_BIAS_BY_SAMPLE_KEY, "1,21,6,50")
                .make();

        builderA.attribute(AS_RAW_RMS_MAPPING_QUALITY_KEY,"29707.00|39366.00|2405.00")
                .attribute(AS_RAW_MAP_QUAL_RANK_SUM_KEY,"|-0.2,1|-2.5,1")
                .attribute(RAW_QUAL_APPROX_KEY,"74")
                .attribute(AS_RAW_QUAL_APPROX_KEY,"|74|0")
                .attribute(AS_RAW_READ_POS_RANK_SUM_KEY,"|2.4,1|1.5,1")
                .attribute(AS_SB_TABLE_KEY,"1,21|3,39|3,11")
                .attribute(AS_VARIANT_DEPTH_KEY,"22|42|0")
                .genotypes(Arrays.asList(g));

        VariantContext vc = builderA.make();

        // generally expect
        //     leading | to be removed
        //     NON_REF allele specific vales
        Assert.assertEquals(VetFieldEnum.ref.getColumnValue(vc, false), "C");
        Assert.assertEquals(VetFieldEnum.alt.getColumnValue(vc, false), "A");
        Assert.assertEquals(VetFieldEnum.AS_RAW_MQ.getColumnValue(vc, false), "29707|39366");
        Assert.assertEquals(VetFieldEnum.AS_RAW_MQRankSum.getColumnValue(vc, false), "-0.2,1");
        Assert.assertEquals(VetFieldEnum.QUALapprox.getColumnValue(vc, false), "74");
        Assert.assertEquals(VetFieldEnum.AS_QUALapprox.getColumnValue(vc, false), "74");
        Assert.assertEquals(VetFieldEnum.AS_RAW_ReadPosRankSum.getColumnValue(vc, false), "2.4,1");
        Assert.assertEquals(VetFieldEnum.AS_SB_TABLE.getColumnValue(vc, false), "1,21|3,39");
        Assert.assertEquals(VetFieldEnum.AS_VarDP.getColumnValue(vc, false), "22|42");

        Assert.assertEquals(VetFieldEnum.call_GQ.getColumnValue(vc, false), "36");
        Assert.assertEquals(VetFieldEnum.call_AD.getColumnValue(vc, false), "22,42");
        Assert.assertEquals(VetFieldEnum.call_GT.getColumnValue(vc, false), "0/1");

        Assert.assertEquals(VetFieldEnum.call_PL.getColumnValue(vc, false), "74,0,34,707,390,467");

    }

    @Test
    public void testAlleleSpecificTwoAlts() {
        // just focusing on the difference for AS annotations for non 0/x genotypes

        VariantContextBuilder builderA =
                new VariantContextBuilder("a","1",10329,10329,
                        Arrays.asList(Allele.REF_C,Allele.ALT_A,Allele.ALT_T, Allele.NON_REF_ALLELE));


        Genotype g = new GenotypeBuilder(SAMPLE_1)
                .alleles(Arrays.asList(Allele.ALT_A, Allele.ALT_T))
                .make();

        builderA.attribute(AS_RAW_MAP_QUAL_RANK_SUM_KEY,"|-0.2,1|-2.5,1")
                .attribute(AS_RAW_READ_POS_RANK_SUM_KEY,"|2.4,1|1.5,1")
                .genotypes(Arrays.asList(g));

        VariantContext vc = builderA.make();

        Assert.assertEquals(VetFieldEnum.AS_RAW_MQRankSum.getColumnValue(vc, false), "");
        Assert.assertEquals(VetFieldEnum.AS_RAW_ReadPosRankSum.getColumnValue(vc, false), "");
    }

    @Test
    public void testNonAlleleSpecificSingleAlt() {
        // Original Record:
        //
        // chr1	12796	.	T	A,<NON_REF>	125.77	.	AS_QUALapprox=|154|0;AS_VarDP=2|6|0;BaseQRankSum=-0.572;DP=8;MQRankSum=-0.572;MQ_DP=8;QUALapprox=154;RAW_GT_COUNT=0,1,0;RAW_MQ=5024.00;RAW_MQandDP=5024,8;ReadPosRankSum=1.067;VarDP=8	GT:AD:DP:GQ:PL:SB	0/1:2,6,0:8:35:154,0,35,161,53,214:2,0,6,0
        //
        VariantContextBuilder builderA =
                new VariantContextBuilder("a","1",12796,12796,
                        Arrays.asList(Allele.REF_T,Allele.ALT_A,Allele.NON_REF_ALLELE));

        Genotype g = new GenotypeBuilder(SAMPLE_1)
                .alleles(Arrays.asList(Allele.REF_T, Allele.ALT_A))
                .PL(new int[]{154,0,35,161,53,214})
                .DP(8)
                .GQ(35)
                .AD(new int[]{2,6,0})
                .attribute(STRAND_BIAS_BY_SAMPLE_KEY, "2,0,6,0")
                .make();

        builderA.attribute(AS_RAW_QUAL_APPROX_KEY,"|154|0")
                .attribute(AS_VARIANT_DEPTH_KEY,"2|6|0")
                .attribute(DEPTH_KEY,"8")
                .attribute(MAP_QUAL_RANK_SUM_KEY,"-0.572")
                .attribute(RAW_QUAL_APPROX_KEY,"154")
                .attribute(RAW_MAPPING_QUALITY_WITH_DEPTH_KEY,"5024,8")
                .attribute(READ_POS_RANK_SUM_KEY,"1.067")
                .attribute(VARIANT_DEPTH_KEY,"8")
                .genotypes(Arrays.asList(g))
        ;


        VariantContext vc = builderA.make();

        Assert.assertEquals(VetFieldEnum.ref.getColumnValue(vc, false), "T");
        Assert.assertEquals(VetFieldEnum.alt.getColumnValue(vc, false), "A");
        Assert.assertEquals(VetFieldEnum.AS_RAW_MQ.getColumnValue(vc, false), "0|5024");
        Assert.assertEquals(VetFieldEnum.AS_RAW_MQRankSum.getColumnValue(vc, false), "-0.572,1");
        Assert.assertEquals(VetFieldEnum.QUALapprox.getColumnValue(vc, false), "154");
        Assert.assertEquals(VetFieldEnum.AS_QUALapprox.getColumnValue(vc, false), "154");
        Assert.assertEquals(VetFieldEnum.AS_RAW_ReadPosRankSum.getColumnValue(vc, false), "1.067,1");
        Assert.assertEquals(VetFieldEnum.AS_SB_TABLE.getColumnValue(vc, false), "2,0|6,0");
        Assert.assertEquals(VetFieldEnum.AS_VarDP.getColumnValue(vc, false), "2|6");

        Assert.assertEquals(VetFieldEnum.call_GQ.getColumnValue(vc, false), "35");
        Assert.assertEquals(VetFieldEnum.call_AD.getColumnValue(vc, false), "2,6");
        Assert.assertEquals(VetFieldEnum.call_GT.getColumnValue(vc, false), "0/1");

         Assert.assertEquals(VetFieldEnum.call_PL.getColumnValue(vc, false), "154,0,35,161,53,214");
    }

    @Test
    public void testNonAlleleSpecificTwoAlts() {
        Allele ALT_CTTT = Allele.create("CTTT".getBytes(), false);
        Allele ALT_CTT = Allele.create("CTT".getBytes(), false);

        //
        // Original
        // chr1	134065	.	C	CTTT,CTT,<NON_REF>	46.38	.	AS_QUALapprox=|31|29|0;AS_VarDP=1|6|3|0;DP=18;MQ=12.49;MQRankSum=0.696;QUALapprox=52;RAW_GT_COUNT=0,0,1;RAW_MQandDP=2808,18;ReadPosRankSum=-0.696;VarDP=10	GT:AD:AF:DP:F1R2:F2R1:GP:GQ:ICNT:MB:PL:PRI:SB:SPL	1/2:1,6,3,0:0.600,0.300:10:1,3,1,0:0,3,2,0:4.6376e+01,2.6376e+01,2.1241e+01,7.6870e+01,6.4666e-02,2.3065e+01,1.5561e+02,7.1204e+01,1.0968e+02,1.4905e+02:18:1,5:1,0,5,4:52,29,21,80,0,23,127,39,78,117:0.00,3.00,6.00,3.00,6.00,6.00,34.77,37.77,37.77,37.77:1,0,4,5:178,0,255
        VariantContextBuilder builderA =
                new VariantContextBuilder("a","1",134065,134065,
                        Arrays.asList(Allele.REF_C,ALT_CTTT, ALT_CTT, Allele.NON_REF_ALLELE));


        Genotype g = new GenotypeBuilder(SAMPLE_1)
                .alleles(Arrays.asList(ALT_CTTT, ALT_CTT))
                .PL(new int[]{52,29,21,80,0,23,127,39,78,117})
                .DP(10)
                .GQ(18)
                .AD(new int[]{1,6,3,0})
                .attribute(STRAND_BIAS_BY_SAMPLE_KEY, "1,0,4,5")
                .make();

        builderA.attribute(AS_RAW_QUAL_APPROX_KEY,"|31|29|0")
                .attribute(AS_VARIANT_DEPTH_KEY,"1|6|3|0")
                .attribute(DEPTH_KEY,"18")
                .attribute(MAP_QUAL_RANK_SUM_KEY,"0.696")
                .attribute(RAW_QUAL_APPROX_KEY,"52")
                .attribute(RAW_MAPPING_QUALITY_WITH_DEPTH_KEY,"2808,18")
                .attribute(READ_POS_RANK_SUM_KEY,"-0.696")
                .attribute(VARIANT_DEPTH_KEY,"10")
                .genotypes(Arrays.asList(g));


        VariantContext vc = builderA.make();

        Assert.assertEquals(VetFieldEnum.ref.getColumnValue(vc, false), "C");
        Assert.assertEquals(VetFieldEnum.alt.getColumnValue(vc, false), "CTTT,CTT");
        Assert.assertEquals(VetFieldEnum.AS_RAW_MQ.getColumnValue(vc, false), "0|2808|2808");
        Assert.assertEquals(VetFieldEnum.AS_RAW_MQRankSum.getColumnValue(vc, false), "");
        Assert.assertEquals(VetFieldEnum.QUALapprox.getColumnValue(vc, false), "52");
        Assert.assertEquals(VetFieldEnum.AS_QUALapprox.getColumnValue(vc, false), "31|29");
        Assert.assertEquals(VetFieldEnum.AS_RAW_ReadPosRankSum.getColumnValue(vc, false), "");
        Assert.assertEquals(VetFieldEnum.AS_SB_TABLE.getColumnValue(vc, false), "1,0|2,2|2,2");
        Assert.assertEquals(VetFieldEnum.AS_VarDP.getColumnValue(vc, false), "1|6|3");

        Assert.assertEquals(VetFieldEnum.call_GQ.getColumnValue(vc, false), "18");
        Assert.assertEquals(VetFieldEnum.call_AD.getColumnValue(vc, false), "1,6,3");
        Assert.assertEquals(VetFieldEnum.call_GT.getColumnValue(vc, false), "1/2");

        Assert.assertEquals(VetFieldEnum.call_PL.getColumnValue(vc, false), "52,29,21,80,0,23,127,39,78,117");
    }

    @Test
    public void testForceNonASLoading(){
        VariantContextBuilder builderA =
                new VariantContextBuilder("a","1",1729859,1729859,
                        Arrays.asList(Allele.REF_C,Allele.ALT_G,Allele.NON_REF_ALLELE));


        Genotype g = new GenotypeBuilder(SAMPLE_1)
                .alleles(Arrays.asList(Allele.REF_C, Allele.ALT_G))
                .AD(new int[]{22,42,0})
                .attribute(STRAND_BIAS_BY_SAMPLE_KEY, "6,7,8,9")
                .make();

        builderA.attribute(AS_RAW_RMS_MAPPING_QUALITY_KEY,"29707.00|39366.00|2405.00")
                .attribute(AS_RAW_MAP_QUAL_RANK_SUM_KEY,"|-0.2,1|-2.5,1")
                .attribute(AS_RAW_QUAL_APPROX_KEY,"|74|0")
                .attribute(AS_RAW_READ_POS_RANK_SUM_KEY,"|2.4,1|1.5,1")
                .attribute(AS_SB_TABLE_KEY,"1,1|2,2|3,3")
                .attribute(MAP_QUAL_RANK_SUM_KEY,"-3.1")
                .attribute(RAW_QUAL_APPROX_KEY,"52")
                .attribute(RAW_MAPPING_QUALITY_WITH_DEPTH_KEY,"1000,18")
                .attribute(READ_POS_RANK_SUM_KEY,"-4.1")

                .genotypes(Arrays.asList(g));

        VariantContext vc = builderA.make();

        Assert.assertEquals(VetFieldEnum.AS_RAW_MQ.getColumnValue(vc, true), "0|1000");
        Assert.assertEquals(VetFieldEnum.AS_RAW_MQRankSum.getColumnValue(vc, true), "-3.1,1");
        Assert.assertEquals(VetFieldEnum.AS_RAW_ReadPosRankSum.getColumnValue(vc, true), "-4.1,1");
        Assert.assertEquals(VetFieldEnum.AS_SB_TABLE.getColumnValue(vc, true), "6,7|8,9");
    }

}
