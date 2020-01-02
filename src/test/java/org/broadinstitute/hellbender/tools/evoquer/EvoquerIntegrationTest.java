package org.broadinstitute.hellbender.tools.evoquer;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerIntegrationTest;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class EvoquerIntegrationTest extends CommandLineProgramTest {


    @Test(groups = {"cloud"})
    public void test3ExomesWithGnarlyGenotyper() throws IOException {
        final String projectID = "broad-dsde-dev";
        final String datasetMapString = "chr20   joint_genotyping_chr20_integration_test  pet_without_gq60_ir_c_sam_st vet";
        final String interval = "chr20:1-417395";
        final File outputVCF = createTempFile("output", ".vcf");

        final File combineGVCFOutput = new File("src/test/resources/large/integration_test_3_sample.vcf");

        final File datasetMapFile = createTempFile("testChr203Exomes", ".dataset_map");
        try ( final PrintWriter writer = new PrintWriter(datasetMapFile) ) {
            writer.println(datasetMapString);
        }

        final String[] args = {
                "--project-id", projectID,
                "--dataset-map", datasetMapFile.getAbsolutePath(),
                "-R", hg38Reference,
                "-L", interval,
                "-O", outputVCF.getAbsolutePath(),
                "--run-query-only", "false",
                "--disable-gnarly-genotyper", "false"
        };

        runCommandLine(args);

        final List<VariantContext> variants = VariantContextTestUtils.streamVcf(outputVCF)
                .sorted(IntervalUtils.LEXICOGRAPHICAL_ORDER_COMPARATOR)
                .collect(Collectors.toList());
        Assert.assertTrue(variants.size() > 0, "Output is empty.");

        final List<VariantContext> expectedVariants = VariantContextTestUtils.streamVcf(combineGVCFOutput)
                .sorted(IntervalUtils.LEXICOGRAPHICAL_ORDER_COMPARATOR)
                .collect(Collectors.toList());
        Assert.assertEquals(variants.size(), expectedVariants.size(), "Number of variants doesn't match expected.");

        HaplotypeCallerIntegrationTest.calculateConcordance(outputVCF, combineGVCFOutput, false, "NA12878");

        for (int i=0; i < variants.size(); i++) {
            VariantContext v = variants.get(i);
            VariantContext expectedV = expectedVariants.get(i);
            Assert.assertEquals(v.getStart(), expectedV.getStart(), "Positions don't match.");
            Assert.assertEquals(v.getReference(), expectedV.getReference(), "Position: " + v.getStart() + " REF allele doesn't match.");
            Assert.assertEquals(v.getAlleles(), expectedV.getAlleles(), "Position: " + v.getStart() + " Alt alleles don't match.");
            Assert.assertEquals(v.getAttribute("AC"), expectedV.getAttribute("AC"), "Position: " + v.getStart() + " does't have matching AC.");
            Assert.assertEquals(v.getAttribute("AS_MQ"), expectedV.getAttribute("AS_MQ"), "Position: " + v.getStart() + " does't have matching AS_MQ.");
            Assert.assertEquals(v.getAttribute("AS_MQRankSum"), expectedV.getAttribute("AS_MQRankSum"), "Position: " + v.getStart() + " does't have matching AS_MQRankSum.");
            //AS_QD doesn't currently match!
            //Assert.assertEquals(v.getAttribute("AS_QD"), expectedV.getAttribute("AS_QD"), "Position: " + v.getStart() + " does't have matching AS_QD.");
            Assert.assertEquals(v.getAttribute("AS_ReadPosRankSum"), expectedV.getAttribute("AS_ReadPosRankSum"), "Position: " + v.getStart() + " does't have matching AS_ReadPosRankSum.");
            //DP doesn't currently match!
            //Assert.assertEquals(v.getAttributeAsInt("DP", 0), expectedV.getAttributeAsInt("DP", 0), "Position: " + v.getStart() + " does't have matching DP.");
            Assert.assertEquals(v.getAttribute("AF"), expectedV.getAttribute("AF"), "Position: " + v.getStart() + " does't have matching AF.");
            Assert.assertEquals(v.getAttribute("AN"), expectedV.getAttribute("AN"), "Position: " + v.getStart() + " does't have matching AN.");
            Assert.assertEquals(v.getAttribute("AS_FS"), expectedV.getAttribute("AS_FS"), "Position: " + v.getStart() + " does't have matching AS_FS.");
            Assert.assertEquals(v.getAttribute("AS_SOR"), expectedV.getAttribute("AS_SOR"), "Position: " + v.getStart() + " does't have matching AS_SOR.");
            //Where is AS_VarDP in expected output?
            //Assert.assertEquals(v.getAttribute("AS_VarDP"), expectedV.getAttribute("AS_VarDP"), "Position: " + v.getStart() + " does't have matching AS_VAR_DP.");
            for (String s : v.getSampleNames()) {
                Genotype g = v.getGenotype(s);
                Genotype expectedG = expectedV.getGenotype(s);
                Assert.assertEquals(g.getGQ(), expectedG.getGQ(), "Position: " + v.getStart() + " Sample: " + s + " does't have matching GQ.");
                if (!g.isHomRef()) {
                    Assert.assertEquals(g.getAnyAttribute("GT"), expectedG.getAnyAttribute("GT"), "Position: " + v.getStart() + " Sample: " + s + " does't have matching GT.");
                    Assert.assertEquals(g.getAD(), expectedG.getAD(), "Position: " + v.getStart() + " Sample: " + s + " does't have matching AD.");
                    Assert.assertEquals(g.getDP(), expectedG.getDP(), "Position: " + v.getStart() + " Sample: " + s + " does't have matching DP.");
                    Assert.assertEquals(g.getPL(), expectedG.getPL(), "Position: " + v.getStart() + " Sample: " + s + " does't have matching PL.");
                    Assert.assertEquals(g.getGQ(), expectedG.getGQ(), "Position: " + v.getStart() + " Sample: " + s + " does't have matching GQ.");
                    if (g.hasAnyAttribute("PGT")) {
                        Assert.assertEquals(g.getAnyAttribute("PGT"), expectedG.getAnyAttribute("PGT"), "Position: " + v.getStart() + " Sample: " + s + " doesn't have matching PGT.");
                        Assert.assertEquals(g.getAnyAttribute("PID"), expectedG.getAnyAttribute("PID"), "Position: " + v.getStart() + " Sample: " + s + " doesn't have matching PID.");
                    }
                }
            }
        }
    }

    @Test(groups = {"cloud"})
    public void multiAllelicSite() throws IOException {
        final String projectID = "broad-dsde-dev";
        final String datasetMapString = "chr20   joint_genotyping_chr20_integration_test  pet_without_gq60_ir_c_sam_st vet";
        final String interval = "chr20:3264573";
        final File outputVCF = createTempFile("output", ".vcf");

        final File datasetMapFile = createTempFile("testChr203Exomes", ".dataset_map");
        try ( final PrintWriter writer = new PrintWriter(datasetMapFile) ) {
            writer.println(datasetMapString);
        }

        final String[] args = {
                "--project-id", projectID,
                "--dataset-map", datasetMapFile.getAbsolutePath(),
                "-R", hg38Reference,
                "-L", interval,
                "-O", outputVCF.getAbsolutePath(),
                "--run-query-only", "false",
                "--disable-gnarly-genotyper", "false"
        };

        runCommandLine(args);

        final List<VariantContext> variants = VariantContextTestUtils.streamVcf(outputVCF)
                .collect(Collectors.toList());
        Assert.assertEquals(variants.size(), 1);

        VariantContext v = variants.get(0);

        Assert.assertEquals(v.getStart(), 3264573, "Positions don't match.");
        Assert.assertEquals(v.getReference().toString(), "CTT*", "Position: " + v.getStart() + " REF allele doesn't match.");
        Assert.assertEquals(v.getAlleles().toString(), "[CTT*, C, CTTT]", "Position: " + v.getStart() + " Alt alleles don't match.");
        Assert.assertEquals(v.getAttribute("AC").toString(), "[1, 1]", "Position: " + v.getStart() + " does't have matching AC.");
        Assert.assertEquals(v.getAttribute("AS_MQ").toString(), "[60.00, 56.94]", "Position: " + v.getStart() + " does't have matching AS_MQ.");
        Assert.assertEquals(v.getAttribute("AS_MQRankSum").toString(), "[0.000, -1.500]", "Position: " + v.getStart() + " does't have matching AS_MQRankSum.");
        Assert.assertEquals(v.getAttribute("AS_QD").toString(), "[0.36, 1.38]", "Position: " + v.getStart() + " does't have matching AS_QD.");
        Assert.assertEquals(v.getAttribute("AS_ReadPosRankSum").toString(), "[1.300, 1.100]", "Position: " + v.getStart() + " does't have matching AS_ReadPosRankSum.");
        //DP doesn't currently match!
        //Assert.assertEquals(v.getAttributeAsInt("DP", 0), 110, "Position: " + v.getStart() + " does't have matching DP.");
        Assert.assertEquals(v.getAttribute("AF").toString(), "[0.167, 0.167]", "Position: " + v.getStart() + " does't have matching AF.");
        Assert.assertEquals(v.getAttribute("AN").toString(), "6", "Position: " + v.getStart() + " does't have matching AN.");
        Assert.assertEquals(v.getAttribute("AS_FS").toString(), "[5.185, 0.000]", "Position: " + v.getStart() + " does't have matching AS_FS.");
        Assert.assertEquals(v.getAttribute("AS_SOR").toString(), "[0.038, 0.671]", "Position: " + v.getStart() + " does't have matching AS_SOR.");
        //Where is AS_VarDP in expected output?
        //Assert.assertEquals(v.getAttribute("AS_VarDP"), ?, "Position: " + v.getStart() + " does't have matching AS_VAR_DP.");

        Genotype g78 = v.getGenotype("NA12878");
        Assert.assertEquals(g78.getAnyAttribute("GT").toString(), "[CTT*, CTT*]");
        Assert.assertEquals(g78.getGQ(), 60);

        Genotype g91 = v.getGenotype("NA12891");
        Assert.assertEquals(g91.getAnyAttribute("GT").toString(), "[CTT*, C]");
        Assert.assertEquals(g91.getGQ(), 20);
        Assert.assertEquals(g91.getPL(), new int[]{ 20, 0, 567, 71, 573, 644 });
        Assert.assertEquals(g91.getDP(), 19);
        Assert.assertEquals(g91.getAD(), new int[]{ 17, 2, 0 });

        Genotype g92 = v.getGenotype("NA12892");
        Assert.assertEquals(g92.getAnyAttribute("GT").toString(), "[CTT*, CTTT]");
        Assert.assertEquals(g92.getGQ(), 87);
        Assert.assertEquals(g92.getPL(), new int[]{ 87, 217, 1105, 0, 856, 832 });
        Assert.assertEquals(g92.getDP(), 51);
        Assert.assertEquals(g92.getAD(), new int[]{ 36, 0, 10 });
    }

    @Test(groups = {"cloud"})
    public void spanningDeletion() throws IOException {
        final String projectID = "broad-dsde-dev";
        final String datasetMapString = "chr20   joint_genotyping_chr20_integration_test  pet_without_gq60_ir_c_sam_st vet";
        final String interval = "chr20:1611703";
        final File outputVCF = createTempFile("output", ".vcf");

        final File datasetMapFile = createTempFile("testChr203Exomes", ".dataset_map");
        try ( final PrintWriter writer = new PrintWriter(datasetMapFile) ) {
            writer.println(datasetMapString);
        }

        final String[] args = {
                "--project-id", projectID,
                "--dataset-map", datasetMapFile.getAbsolutePath(),
                "-R", hg38Reference,
                "-L", interval,
                "-O", outputVCF.getAbsolutePath(),
                "--run-query-only", "false",
                "--disable-gnarly-genotyper", "false"
        };

        runCommandLine(args);

        final List<VariantContext> variants = VariantContextTestUtils.streamVcf(outputVCF)
                .collect(Collectors.toList());
        Assert.assertEquals(variants.size(), 1);

        //Spanning deletions don't work yet
        //Assert.assertEquals(variants.get(0).getAlleles().toString(), "[G*, A, *]");
    }

    @Test
    public void testAlleleSubsetting() {
        try {

            final String avrolist = toolsTestDir + "evoquer/subsetting-data.avro";
            final String interval = "chr20:62065822";
            final File outputVCF = createTempFile("output", ".vcf");

            final File avrolistFile = createTempFile("avrolist", "");
            try (final PrintWriter writer = new PrintWriter(avrolistFile)) {
                writer.println(avrolist);
            }

            final String[] args = {
                    "-R", hg38Reference,
                    "--precomputed-query-results", avrolistFile.getAbsolutePath(),
                    "--sample-list-file", toolsTestDir + "evoquer/subsetting-samples",
                    "-L", interval,
                    "-O", outputVCF.getAbsolutePath(),
                    "--run-query-only", "false",
                    "--subset-alleles", "6",
                    "--disable-gnarly-genotyper", "false"
            };

            runCommandLine(args);

            final List<VariantContext> variants = VariantContextTestUtils.streamVcf(outputVCF)
                    .collect(Collectors.toList());
            Assert.assertEquals(variants.size(), 1);
            VariantContext variant = variants.get(0);

            Allele ref = Allele.create("GGCCGCCGCC", true);
            Assert.assertEquals(ref, variant.getReference());
            List<Allele> expectedAlleles = Arrays.asList(
                    ref,
                    Allele.create("AGCCGCCGCC"),
                    Allele.create("G"),
                    Allele.create("GGCC"),
                    Allele.create("GGCCGCC"),
                    Allele.create("GGCCGCCGCCGCCGCC"),
                    Allele.create("GGCCGCCGCCGCCGCCGCC")
            );

            Assert.assertEquals(variant.getAlleles().stream().collect(Collectors.toSet()), expectedAlleles.stream().collect(Collectors.toSet()));

            Map<Allele, Integer> expectedToActualIndexMap = expectedAlleles.stream()
                    .collect(Collectors.toMap(Function.identity(), allele -> variant.getAlleleIndex(allele), (a, b)->a, () -> new HashMap<>()));


            // one allele subsetted
            verifyGenotype(variant.getGenotypes("M3399").get(0), "GGCCGCCGCC/.",
                    Arrays.asList(34,0,68,0,68,68,0,68,68,68,0,68,68,68,68,0,68,68,68,68,68,0,68,68,68,68,68,68),
                    Arrays.asList(3,0,0,0,0,0,0));
            // both alleles subsetted
            verifyGenotype(variant.getGenotypes("102185").get(0), "./.",
                    Arrays.asList(85,10,0,10,0,0,10,0,0,0,10,0,0,0,0,10,0,0,0,0,0,10,0,0,0,0,0,0),
                    Arrays.asList(0,0,0,0,0,0,0));
            // no alleles subsetted
            verifyGenotype(variant.getGenotypes("8917408").get(0), "GGCCGCCGCC/AGCCGCCGCC",
                    Arrays.asList(141,174,423,174,423,423,174,423,423,423,0,249,249,249,225,174,423,423,423,249,423,174,423,423,423,249,423,423),
                    Arrays.asList(11,0,0,0,8,0,0));

        } catch (Exception ex) {
            Assert.fail("caught exception", ex);
        }
    }


    private void verifyGenotype(Genotype sample, String expectedGenotype, List<Integer> expectedPLs, List<Integer> expectedADs) {
        Assert.assertEquals(sample.getGenotypeString(), expectedGenotype);
        List<Integer> PLs = Arrays.stream(sample.getPL()).sorted().boxed().collect(Collectors.toList());
        expectedPLs.sort(Integer::compareTo);
        Assert.assertEquals(PLs, expectedPLs);
        List<Integer> ADs = Arrays.stream(sample.getAD()).sorted().boxed().collect(Collectors.toList());
        expectedADs.sort(Integer::compareTo);
        Assert.assertEquals(ADs, expectedADs);
    }

}
