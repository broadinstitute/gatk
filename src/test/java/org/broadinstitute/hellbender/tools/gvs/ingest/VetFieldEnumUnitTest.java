package org.broadinstitute.hellbender.tools.gvs.ingest;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.testng.annotations.Test;

import java.util.Arrays;
import org.testng.Assert;

import static org.broadinstitute.hellbender.utils.variant.GATKVCFConstants.STRAND_BIAS_BY_SAMPLE_KEY;

public class VetFieldEnumUnitTest {
    private static final String SAMPLE_1 = "NA1";

    @Test
    public void testAlleleSpecific() {
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
                .attribute(STRAND_BIAS_BY_SAMPLE_KEY, new int[]{1,21,6,50})
                .make();

        // TODO: use constants!!!
        builderA.attribute("AS_RAW_MQ","29707.00|39366.00|2405.00")
                .attribute("AS_RAW_MQRankSum","|-0.2,1|-2.5,1")
                .attribute("QUALapprox","74")
                .attribute("AS_QUALapprox","|74|0")
                .attribute("AS_RAW_ReadPosRankSum","|2.4,1|1.5,1")
                .attribute("AS_SB_TABLE","1,21|3,39|3,11")
                .attribute("AS_VarDP","22|42|0")
                .genotypes(Arrays.asList(g));

        VariantContext vc = builderA.make();

        // generally expect
        //     leading | to be removed
        //     NON_REF allele specific vales
        Assert.assertEquals(VetFieldEnum.ref.getColumnValue(vc), "C");
        Assert.assertEquals(VetFieldEnum.alt.getColumnValue(vc), "A");
        Assert.assertEquals(VetFieldEnum.AS_RAW_MQ.getColumnValue(vc), "29707|39366");
        Assert.assertEquals(VetFieldEnum.AS_RAW_MQRankSum.getColumnValue(vc), "-0.2,1");
        Assert.assertEquals(VetFieldEnum.QUALapprox.getColumnValue(vc), "74");
        Assert.assertEquals(VetFieldEnum.AS_QUALapprox.getColumnValue(vc), "74");
        Assert.assertEquals(VetFieldEnum.AS_RAW_ReadPosRankSum.getColumnValue(vc), "2.4,1");
        Assert.assertEquals(VetFieldEnum.AS_SB_TABLE.getColumnValue(vc), "1,21|3,39");
        Assert.assertEquals(VetFieldEnum.AS_VarDP.getColumnValue(vc), "22|42");

        Assert.assertEquals(VetFieldEnum.call_GQ.getColumnValue(vc), "36");
        Assert.assertEquals(VetFieldEnum.call_AD.getColumnValue(vc), "22,42");
        Assert.assertEquals(VetFieldEnum.call_GT.getColumnValue(vc), "0/1");

        Assert.assertEquals(VetFieldEnum.call_PL.getColumnValue(vc), "74,0,34,707,390,467");

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
                .attribute(STRAND_BIAS_BY_SAMPLE_KEY, "2,0,6,0") // KCIBUL: why doesn't this work with: new int[]{2,0,6,0}  What is right?
                .make();

        builderA.attribute("AS_QUALapprox","|154|0")
                .attribute("AS_VarDP","2|6|0")
                .attribute("BaseQRankSum","-0.572")
                .attribute("DP","8")
                .attribute("MQRankSum","-0.572")
                .attribute("MQ_DP","8")
                .attribute("QUALapprox","154")
                .attribute("RAW_GT_COUNT","0,1,0")
                .attribute("RAW_MQ","5024.00")
                .attribute("RAW_MQandDP","5024,8")
                .attribute("ReadPosRankSum","1.067")
                .attribute("VarDP","8")
                .genotypes(Arrays.asList(g))
        ;


        VariantContext vc = builderA.make();

        Assert.assertEquals(VetFieldEnum.ref.getColumnValue(vc), "T");
        Assert.assertEquals(VetFieldEnum.alt.getColumnValue(vc), "A");
        Assert.assertEquals(VetFieldEnum.AS_RAW_MQ.getColumnValue(vc), "2512|2512");
        Assert.assertEquals(VetFieldEnum.AS_RAW_MQRankSum.getColumnValue(vc), "-0.572,1");
        Assert.assertEquals(VetFieldEnum.QUALapprox.getColumnValue(vc), "154");
        Assert.assertEquals(VetFieldEnum.AS_QUALapprox.getColumnValue(vc), "154");
        Assert.assertEquals(VetFieldEnum.AS_RAW_ReadPosRankSum.getColumnValue(vc), "1.067,1");
        Assert.assertEquals(VetFieldEnum.AS_SB_TABLE.getColumnValue(vc), "2,0|6,0");
        Assert.assertEquals(VetFieldEnum.AS_VarDP.getColumnValue(vc), "2|6");

        Assert.assertEquals(VetFieldEnum.call_GQ.getColumnValue(vc), "35");
        Assert.assertEquals(VetFieldEnum.call_AD.getColumnValue(vc), "2,6");
        Assert.assertEquals(VetFieldEnum.call_GT.getColumnValue(vc), "0/1");

         Assert.assertEquals(VetFieldEnum.call_PL.getColumnValue(vc), "154,0,35,161,53,214");
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
                .attribute(STRAND_BIAS_BY_SAMPLE_KEY, "1,0,4,5") // KCIBUL: why doesn't this work with: new int[]{2,0,6,0}  What is right?
                .make();

        builderA.attribute("AS_QUALapprox","|31|29|0")
                .attribute("AS_VarDP","1|6|3|0")
//                .attribute("BaseQRankSum","-0.572")
                .attribute("DP","18")
                .attribute("MQRankSum","0.696")
                .attribute("QUALapprox","52")
                .attribute("RAW_GT_COUNT","0,0,1")
//                .attribute("RAW_MQ",".00")
                .attribute("RAW_MQandDP","2808,18")
                .attribute("ReadPosRankSum","-0.696")
                .attribute("VarDP","10")
                .genotypes(Arrays.asList(g));


        VariantContext vc = builderA.make();

        Assert.assertEquals(VetFieldEnum.ref.getColumnValue(vc), "C");
        Assert.assertEquals(VetFieldEnum.alt.getColumnValue(vc), "CTTT,CTT");
        Assert.assertEquals(VetFieldEnum.AS_RAW_MQ.getColumnValue(vc), "936|936|936");
        Assert.assertEquals(VetFieldEnum.AS_RAW_MQRankSum.getColumnValue(vc), "0.696,1|0.696,1");
        Assert.assertEquals(VetFieldEnum.QUALapprox.getColumnValue(vc), "52");
        Assert.assertEquals(VetFieldEnum.AS_QUALapprox.getColumnValue(vc), "31|29");
        Assert.assertEquals(VetFieldEnum.AS_RAW_ReadPosRankSum.getColumnValue(vc), "-0.696,1|-0.696,1");
        Assert.assertEquals(VetFieldEnum.AS_SB_TABLE.getColumnValue(vc), "1,0|2,2|2,2");
        Assert.assertEquals(VetFieldEnum.AS_VarDP.getColumnValue(vc), "1|6|3");

        Assert.assertEquals(VetFieldEnum.call_GQ.getColumnValue(vc), "18");
        Assert.assertEquals(VetFieldEnum.call_AD.getColumnValue(vc), "1,6,3");
        Assert.assertEquals(VetFieldEnum.call_GT.getColumnValue(vc), "1/2");

         Assert.assertEquals(VetFieldEnum.call_PL.getColumnValue(vc), "52,29,21,80,0,23,127,39,78,117");
    }

}
