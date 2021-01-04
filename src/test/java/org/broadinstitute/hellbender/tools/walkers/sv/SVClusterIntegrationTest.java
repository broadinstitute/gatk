package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.StructuralVariantType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.cluster.SVClusterEngineArgumentsCollection;
import org.spark_project.guava.collect.Lists;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.List;

public class SVClusterIntegrationTest extends CommandLineProgramTest {

    private static final String SEQUENCE_DICT = GATKBaseTest.b38_reference_20_21.replace(".fasta", ".dict");

    @Test
    public void testDefragmentation() {
        final File output = createTempFile("defragmented", ".vcf");

        final String inputVcfPath = getToolTestDataDir() + "1kgp_test.cnvs.vcf.gz";
        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addOutput(output)
                .add(StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME, SEQUENCE_DICT)
                .addVCF(inputVcfPath)
                .add(SVCluster.ALGORITHM_LONG_NAME, SVCluster.CLUSTER_ALGORITHM.DEFRAGMENT_CNV)
                .add(SVCluster.DEFRAG_PADDING_FRACTION_LONG_NAME, 0.25)
                .add(SVClusterEngineArgumentsCollection.DEPTH_SAMPLE_OVERLAP_FRACTION_NAME, 0.5);

        runCommandLine(args, SVCluster.class.getSimpleName());

        final VCFHeader inputHeader =  VariantContextTestUtils.getVCFHeader(inputVcfPath);
        final Pair<VCFHeader, List<VariantContext>> vcf = VariantContextTestUtils.readEntireVCFIntoMemory(output.getAbsolutePath());
        final VCFHeader header = vcf.getKey();
        final List<VariantContext> records = vcf.getValue();

        Assert.assertEquals(header.getSampleNamesInOrder(), inputHeader.getSampleNamesInOrder());
        Assert.assertEquals(header.getSequenceDictionary(), inputHeader.getSequenceDictionary());

        Assert.assertEquals(records.size(), 408);

        // Check for one record
        boolean foundExpectedDefragmentedRecord = false;
        for (final VariantContext variant : records) {
            Assert.assertTrue(variant.hasAttribute(GATKSVVCFConstants.CLUSTER_MEMBER_IDS_KEY));
            Assert.assertFalse(variant.getAttributeAsStringList(GATKSVVCFConstants.CLUSTER_MEMBER_IDS_KEY, null).isEmpty());
            Assert.assertTrue(variant.hasAttribute(GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE));
            Assert.assertTrue(variant.getAttributeAsStringList(GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE, null).contains(GATKSVVCFConstants.DEPTH_ALGORITHM));
            if (variant.getStart() == 25928335 && variant.getEnd() == 29403000) {
                foundExpectedDefragmentedRecord = true;
                final List<Allele> alts = variant.getAlternateAlleles();
                Assert.assertEquals(alts.size(), 1);
                Assert.assertEquals(alts.get(0), Allele.SV_SIMPLE_DUP);
                Assert.assertEquals(variant.getStructuralVariantType(), StructuralVariantType.DUP);
                for (final Genotype g : variant.getGenotypes()) {
                    if (g.getSampleName().equals("HG00129")) {
                        Assert.assertTrue(g.isHet());
                    } else {
                        Assert.assertTrue(g.isHomRef());
                    }
                }
                break;
            }
        }
        Assert.assertTrue(foundExpectedDefragmentedRecord);
    }

    @DataProvider(name = "testMergeData")
    public Object[][] testMergeData() {
        return new Object[][]{
                {true},
                {false}
        };
    }

    @Test(dataProvider= "testMergeData")
    public void testMergeHelper(final boolean omitMembers) {
        final File output = createTempFile("merged", ".vcf");
        ArgumentsBuilder args = new ArgumentsBuilder()
                .addOutput(output)
                .addFlag(StandardArgumentDefinitions.DISABLE_SEQUENCE_DICT_VALIDATION_NAME)
                .add(StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME, SEQUENCE_DICT)
                .addVCF(getToolTestDataDir() + "1kgp_test.cnvs.vcf.gz")
                .addVCF(getToolTestDataDir() + "HG00096.manta.vcf.gz")
                .addVCF(getToolTestDataDir() + "HG00096.wham.vcf.gz")
                .addVCF(getToolTestDataDir() + "HG00129.manta.vcf.gz")
                .addVCF(getToolTestDataDir() + "HG00129.wham.vcf.gz")
                .addVCF(getToolTestDataDir() + "HG00140.manta.vcf.gz")
                .addVCF(getToolTestDataDir() + "HG00140.wham.vcf.gz")
                .add(SVCluster.ALGORITHM_LONG_NAME, SVCluster.CLUSTER_ALGORITHM.SINGLE_LINKAGE)
                .add(SVClusterEngineArgumentsCollection.DEPTH_SAMPLE_OVERLAP_FRACTION_NAME, 0)
                .add(SVClusterEngineArgumentsCollection.DEPTH_INTERVAL_OVERLAP_FRACTION_NAME, 1)
                .add(SVClusterEngineArgumentsCollection.DEPTH_BREAKEND_WINDOW_NAME, 0)
                .add(SVClusterEngineArgumentsCollection.MIXED_SAMPLE_OVERLAP_FRACTION_NAME, 0)
                .add(SVClusterEngineArgumentsCollection.MIXED_INTERVAL_OVERLAP_FRACTION_NAME, 1)
                .add(SVClusterEngineArgumentsCollection.MIXED_BREAKEND_WINDOW_NAME, 0)
                .add(SVClusterEngineArgumentsCollection.PESR_SAMPLE_OVERLAP_FRACTION_NAME, 0)
                .add(SVClusterEngineArgumentsCollection.PESR_INTERVAL_OVERLAP_FRACTION_NAME, 1)
                .add(SVClusterEngineArgumentsCollection.PESR_BREAKEND_WINDOW_NAME, 0);
        if (omitMembers) {
            args = args.addFlag(SVCluster.OMIT_MEMBERS_LONG_NAME);
        }

        runCommandLine(args, SVCluster.class.getSimpleName());

        final Pair<VCFHeader, List<VariantContext>> vcf = VariantContextTestUtils.readEntireVCFIntoMemory(output.getAbsolutePath());
        final VCFHeader header = vcf.getKey();
        final List<VariantContext> records = vcf.getValue();

        Assert.assertEquals(header.getSampleNamesInOrder(), Lists.newArrayList("HG00096", "HG00129", "HG00140", "NA18945", "NA18956"));

        Assert.assertEquals(records.size(), 1793);

        // Check for one record
        int expectedRecordsFound = 0;
        for (final VariantContext variant : records) {
            Assert.assertEquals(!variant.hasAttribute(GATKSVVCFConstants.CLUSTER_MEMBER_IDS_KEY), omitMembers);
            Assert.assertTrue(variant.hasAttribute(GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE));
            if (variant.getID().equals("manta_HG00096_8620")) {
                expectedRecordsFound++;
                final List<String> algorithms = variant.getAttributeAsStringList(GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE, null);
                Assert.assertEquals(algorithms.size(), 1);
                Assert.assertEquals(algorithms.get(0), "manta");
                if (!omitMembers) {
                    final List<String> members = variant.getAttributeAsStringList(GATKSVVCFConstants.CLUSTER_MEMBER_IDS_KEY, null);
                    Assert.assertEquals(members.size(), 2);
                }
                final List<Allele> alts = variant.getAlternateAlleles();
                Assert.assertEquals(alts.size(), 1);
                Assert.assertEquals(alts.get(0), Allele.SV_SIMPLE_INS);
                Assert.assertEquals(variant.getStructuralVariantType(), StructuralVariantType.INS);
                for (final Genotype g : variant.getGenotypes()) {
                    if (g.getSampleName().equals("HG00096") || g.getSampleName().equals("HG00140")) {
                        Assert.assertTrue(g.isHet());
                    } else {
                        Assert.assertTrue(g.isNoCall());
                    }
                }
            }
        }
        Assert.assertEquals(expectedRecordsFound, 1);
    }

    @Test
    public void testClusterSingleLinkage() {
        final File output = createTempFile("single_linkage_cluster", ".vcf");
        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addOutput(output)
                .addVCF(getToolTestDataDir() + "1kgp_test.merged.vcf.gz")
                .add(SVCluster.ALGORITHM_LONG_NAME, SVCluster.CLUSTER_ALGORITHM.SINGLE_LINKAGE)
                .add(SVClusterEngineArgumentsCollection.DEPTH_SAMPLE_OVERLAP_FRACTION_NAME, 0)
                .add(SVClusterEngineArgumentsCollection.DEPTH_INTERVAL_OVERLAP_FRACTION_NAME, 0.5)
                .add(SVClusterEngineArgumentsCollection.DEPTH_BREAKEND_WINDOW_NAME, 2000)
                .add(SVClusterEngineArgumentsCollection.MIXED_SAMPLE_OVERLAP_FRACTION_NAME, 0)
                .add(SVClusterEngineArgumentsCollection.MIXED_INTERVAL_OVERLAP_FRACTION_NAME, 0.1)
                .add(SVClusterEngineArgumentsCollection.MIXED_BREAKEND_WINDOW_NAME, 2000)
                .add(SVClusterEngineArgumentsCollection.PESR_SAMPLE_OVERLAP_FRACTION_NAME, 0)
                .add(SVClusterEngineArgumentsCollection.PESR_INTERVAL_OVERLAP_FRACTION_NAME, 0.1)
                .add(SVClusterEngineArgumentsCollection.PESR_BREAKEND_WINDOW_NAME, 500);

        runCommandLine(args, SVCluster.class.getSimpleName());

        final Pair<VCFHeader, List<VariantContext>> vcf = VariantContextTestUtils.readEntireVCFIntoMemory(output.getAbsolutePath());
        final VCFHeader header = vcf.getKey();
        final List<VariantContext> records = vcf.getValue();

        Assert.assertEquals(header.getSampleNamesInOrder(), Lists.newArrayList("HG00096", "HG00129", "HG00140", "NA18945", "NA18956"));

        Assert.assertEquals(records.size(), 1355);

        // Check for one record
        int expectedRecordsFound = 0;
        for (final VariantContext variant : records) {
            Assert.assertTrue(variant.hasAttribute(GATKSVVCFConstants.CLUSTER_MEMBER_IDS_KEY));
            Assert.assertTrue(variant.hasAttribute(GATKSVVCFConstants.SVLEN));
            Assert.assertTrue(variant.hasAttribute(GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE));
            if (variant.getID().equals("manta_HG00096_8505")) {
                expectedRecordsFound++;
                Assert.assertEquals(variant.getStart(), 1408502);
                Assert.assertEquals(variant.getEnd(), 1410171);
                final int length = variant.getAttributeAsInt(GATKSVVCFConstants.SVLEN, 0);
                Assert.assertEquals(length, 1670);
                final List<String> algorithms = variant.getAttributeAsStringList(GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE, null);
                Assert.assertEquals(algorithms.size(), 2);
                Assert.assertTrue(algorithms.contains("manta"));
                Assert.assertTrue(algorithms.contains("wham"));
                final List<String> members = variant.getAttributeAsStringList(GATKSVVCFConstants.CLUSTER_MEMBER_IDS_KEY, null);
                Assert.assertEquals(members.size(), 3);
                final List<Allele> alts = variant.getAlternateAlleles();
                Assert.assertEquals(alts.size(), 1);
                Assert.assertEquals(alts.get(0), Allele.SV_SIMPLE_DEL);
                Assert.assertEquals(variant.getStructuralVariantType(), StructuralVariantType.DEL);
                for (final Genotype g : variant.getGenotypes()) {
                    if (g.getSampleName().equals("HG00096") || g.getSampleName().equals("HG00140")) {
                        Assert.assertTrue(g.isHomVar());
                    } else if (g.getSampleName().equals("HG00129")) {
                        Assert.assertTrue(g.isHet());
                    } else {
                        Assert.assertTrue(g.isHomRef());
                    }
                }
            }
        }
        Assert.assertEquals(expectedRecordsFound, 1);
    }

    @DataProvider(name = "testClusterMaxCliqueData")
    public Object[][] testClusterMaxCliqueData() {
        return new Object[][]{
                {true},
                {false}
        };
    }

    @Test(dataProvider= "testClusterMaxCliqueData")
    public void testClusterMaxClique(final boolean fastMode) {
        final File output = createTempFile("max_clique_cluster", ".vcf");
        ArgumentsBuilder args = new ArgumentsBuilder()
                .addOutput(output)
                .addVCF(getToolTestDataDir() + "1kgp_test.merged.vcf.gz")
                .add(SVCluster.ALGORITHM_LONG_NAME, SVCluster.CLUSTER_ALGORITHM.MAX_CLIQUE)
                .add(SVClusterEngineArgumentsCollection.DEPTH_SAMPLE_OVERLAP_FRACTION_NAME, 0)
                .add(SVClusterEngineArgumentsCollection.DEPTH_INTERVAL_OVERLAP_FRACTION_NAME, 0.5)
                .add(SVClusterEngineArgumentsCollection.DEPTH_BREAKEND_WINDOW_NAME, 2000)
                .add(SVClusterEngineArgumentsCollection.MIXED_SAMPLE_OVERLAP_FRACTION_NAME, 0)
                .add(SVClusterEngineArgumentsCollection.MIXED_INTERVAL_OVERLAP_FRACTION_NAME, 0.1)
                .add(SVClusterEngineArgumentsCollection.MIXED_BREAKEND_WINDOW_NAME, 2000)
                .add(SVClusterEngineArgumentsCollection.PESR_SAMPLE_OVERLAP_FRACTION_NAME, 0)
                .add(SVClusterEngineArgumentsCollection.PESR_INTERVAL_OVERLAP_FRACTION_NAME, 0.1)
                .add(SVClusterEngineArgumentsCollection.PESR_BREAKEND_WINDOW_NAME, 500);
        if (fastMode) {
            args = args.addFlag(SVCluster.FAST_MODE_LONG_NAME);
        }

        runCommandLine(args, SVCluster.class.getSimpleName());

        final Pair<VCFHeader, List<VariantContext>> vcf = VariantContextTestUtils.readEntireVCFIntoMemory(output.getAbsolutePath());
        final VCFHeader header = vcf.getKey();
        final List<VariantContext> records = vcf.getValue();

        Assert.assertEquals(header.getSampleNamesInOrder(), Lists.newArrayList("HG00096", "HG00129", "HG00140", "NA18945", "NA18956"));

        Assert.assertEquals(records.size(), 1370);

        // Check for one record
        int expectedRecordsFound = 0;
        for (final VariantContext variant : records) {
            Assert.assertTrue(variant.hasAttribute(GATKSVVCFConstants.CLUSTER_MEMBER_IDS_KEY));
            Assert.assertTrue(variant.hasAttribute(GATKSVVCFConstants.SVLEN));
            Assert.assertTrue(variant.hasAttribute(GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE));
            if (variant.getID().equals("wham_HG00140_3361")) {
                expectedRecordsFound++;
                Assert.assertEquals(variant.getStart(), 10753974);
                Assert.assertEquals(variant.getEnd(), 10754620);
                final int length = variant.getAttributeAsInt(GATKSVVCFConstants.SVLEN, 0);
                Assert.assertEquals(length, 647);
                final List<String> algorithms = variant.getAttributeAsStringList(GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE, null);
                Assert.assertEquals(algorithms.size(), 1);
                Assert.assertTrue(algorithms.contains("wham"));
                final List<String> members = variant.getAttributeAsStringList(GATKSVVCFConstants.CLUSTER_MEMBER_IDS_KEY, null);
                Assert.assertEquals(members.size(), 2);
                final List<Allele> alts = variant.getAlternateAlleles();
                Assert.assertEquals(alts.size(), 1);
                Assert.assertEquals(alts.get(0), Allele.SV_SIMPLE_DEL);
                Assert.assertEquals(variant.getStructuralVariantType(), StructuralVariantType.DEL);
                for (final Genotype g : variant.getGenotypes()) {
                    if (g.getSampleName().equals("HG00140")) {
                        Assert.assertTrue(g.isHet());
                    } else if (!fastMode) {
                        Assert.assertTrue(g.isHomRef());
                    }
                }
            }
        }
        Assert.assertEquals(expectedRecordsFound, 1);
    }

    @Test
    public void testClusterSampleOverlap() {
        final File output = createTempFile("cluster_sample_overlap", ".vcf");
        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addOutput(output)
                .addVCF(getToolTestDataDir() + "1kgp_test.merged.vcf.gz")
                .add(SVCluster.ALGORITHM_LONG_NAME, SVCluster.CLUSTER_ALGORITHM.SINGLE_LINKAGE)
                .add(SVClusterEngineArgumentsCollection.DEPTH_SAMPLE_OVERLAP_FRACTION_NAME, 0.5)
                .add(SVClusterEngineArgumentsCollection.DEPTH_INTERVAL_OVERLAP_FRACTION_NAME, 0.5)
                .add(SVClusterEngineArgumentsCollection.DEPTH_BREAKEND_WINDOW_NAME, 2000)
                .add(SVClusterEngineArgumentsCollection.MIXED_SAMPLE_OVERLAP_FRACTION_NAME, 0.5)
                .add(SVClusterEngineArgumentsCollection.MIXED_INTERVAL_OVERLAP_FRACTION_NAME, 0.1)
                .add(SVClusterEngineArgumentsCollection.MIXED_BREAKEND_WINDOW_NAME, 2000)
                .add(SVClusterEngineArgumentsCollection.PESR_SAMPLE_OVERLAP_FRACTION_NAME, 0.5)
                .add(SVClusterEngineArgumentsCollection.PESR_INTERVAL_OVERLAP_FRACTION_NAME, 0.1)
                .add(SVClusterEngineArgumentsCollection.PESR_BREAKEND_WINDOW_NAME, 500);

        runCommandLine(args, SVCluster.class.getSimpleName());

        final Pair<VCFHeader, List<VariantContext>> vcf = VariantContextTestUtils.readEntireVCFIntoMemory(output.getAbsolutePath());
        final VCFHeader header = vcf.getKey();
        final List<VariantContext> records = vcf.getValue();

        Assert.assertEquals(header.getSampleNamesInOrder(), Lists.newArrayList("HG00096", "HG00129", "HG00140", "NA18945", "NA18956"));

        Assert.assertEquals(records.size(), 1618);

        // Check for one record
        int expectedRecordsFound = 0;
        for (final VariantContext variant : records) {
            Assert.assertTrue(variant.hasAttribute(GATKSVVCFConstants.CLUSTER_MEMBER_IDS_KEY));
            Assert.assertTrue(variant.hasAttribute(GATKSVVCFConstants.SVLEN));
            Assert.assertTrue(variant.hasAttribute(GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE));
            if (variant.getID().equals("batch1_DEL_549336")) {
                expectedRecordsFound++;
                Assert.assertEquals(variant.getStart(), 28754914);
                Assert.assertEquals(variant.getEnd(), 28757851);
                final int length = variant.getAttributeAsInt(GATKSVVCFConstants.SVLEN, 0);
                Assert.assertEquals(length, 2938);
                final List<String> algorithms = variant.getAttributeAsStringList(GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE, null);
                Assert.assertEquals(algorithms.size(), 2);
                Assert.assertTrue(algorithms.contains(GATKSVVCFConstants.DEPTH_ALGORITHM));
                Assert.assertTrue(algorithms.contains("manta"));
                final List<String> members = variant.getAttributeAsStringList(GATKSVVCFConstants.CLUSTER_MEMBER_IDS_KEY, null);
                Assert.assertEquals(members.size(), 2);
                final List<Allele> alts = variant.getAlternateAlleles();
                Assert.assertEquals(alts.size(), 1);
                Assert.assertEquals(alts.get(0), Allele.SV_SIMPLE_DEL);
                Assert.assertEquals(variant.getStructuralVariantType(), StructuralVariantType.DEL);
                for (final Genotype g : variant.getGenotypes()) {
                    if (g.getSampleName().equals("HG00129")) {
                        Assert.assertTrue(g.isHet());
                    } else if (g.getSampleName().equals("HG00140")) {
                        Assert.assertTrue(g.isHom());
                    } else {
                        Assert.assertTrue(g.isHomRef());
                    }
                }
            }
        }
        Assert.assertEquals(expectedRecordsFound, 1);
    }

}