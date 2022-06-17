package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVCallRecordUtils;
import org.broadinstitute.hellbender.tools.sv.SVTestUtils;
import org.broadinstitute.hellbender.tools.sv.cluster.*;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.reference.ReferenceUtils;
import org.broadinstitute.hellbender.utils.variant.VariantContextGetters;
import org.spark_project.guava.collect.Lists;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

public class SVClusterIntegrationTest extends CommandLineProgramTest {

    private static final String REFERENCE_PATH = GATKBaseTest.hg38Reference;

    @Test
    public void testDefragmentation() {
        final File output = createTempFile("defragmented", ".vcf");

        final String inputVcfPath = getToolTestDataDir() + "1kgp_test.cnvs.vcf.gz";
        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addOutput(output)
                .add(StandardArgumentDefinitions.REFERENCE_LONG_NAME, REFERENCE_PATH)
                .addVCF(inputVcfPath)
                .add(SVCluster.PLOIDY_TABLE_LONG_NAME, getToolTestDataDir() + "1kgp.batch1.ploidy.tsv")
                .add(SVCluster.VARIANT_PREFIX_LONG_NAME, "SVx")
                .add(SVCluster.ALGORITHM_LONG_NAME, SVCluster.CLUSTER_ALGORITHM.DEFRAGMENT_CNV)
                .add(SVCluster.DEFRAG_PADDING_FRACTION_LONG_NAME, 0.25)
                .add(SVClusterEngineArgumentsCollection.DEPTH_SAMPLE_OVERLAP_FRACTION_NAME, 0.5);

        runCommandLine(args, SVCluster.class.getSimpleName());

        final VCFHeader inputHeader =  VariantContextTestUtils.getVCFHeader(inputVcfPath);
        final Pair<VCFHeader, List<VariantContext>> vcf = VariantContextTestUtils.readEntireVCFIntoMemory(output.getAbsolutePath());
        final VCFHeader header = vcf.getKey();
        final List<VariantContext> records = vcf.getValue();

        Assert.assertEquals(header.getSampleNamesInOrder(), inputHeader.getSampleNamesInOrder());
        Assert.assertEquals(header.getSequenceDictionary().size(), inputHeader.getSequenceDictionary().size());

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
                        Assert.assertEquals(VariantContextGetters.getAttributeAsInt(g, GATKSVVCFConstants.COPY_NUMBER_FORMAT, -1), 3);
                        Assert.assertEquals(VariantContextGetters.getAttributeAsInt(g, GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, -1), 2);
                    } else {
                        Assert.assertEquals(VariantContextGetters.getAttributeAsInt(g, GATKSVVCFConstants.COPY_NUMBER_FORMAT, -1), 2);
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
                .add(StandardArgumentDefinitions.REFERENCE_LONG_NAME, REFERENCE_PATH)
                .addVCF(getToolTestDataDir() + "1kgp_test.cnvs.vcf.gz")
                .addVCF(getToolTestDataDir() + "HG00096.manta.vcf.gz")
                .addVCF(getToolTestDataDir() + "HG00096.wham.vcf.gz")
                .addVCF(getToolTestDataDir() + "HG00129.manta.vcf.gz")
                .addVCF(getToolTestDataDir() + "HG00129.wham.vcf.gz")
                .addVCF(getToolTestDataDir() + "HG00140.manta.vcf.gz")
                .addVCF(getToolTestDataDir() + "HG00140.wham.vcf.gz")
                .add(SVCluster.PLOIDY_TABLE_LONG_NAME, getToolTestDataDir() + "1kgp.batch1.ploidy.tsv")
                .add(SVCluster.VARIANT_PREFIX_LONG_NAME, "SVx")
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
            if (variant.getID().equals("SVx00000065")) {
                expectedRecordsFound++;
                Assert.assertEquals(variant.getStart(), 8835844);
                Assert.assertEquals(variant.getEnd(), 8835844);
                final List<String> algorithms = variant.getAttributeAsStringList(GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE, null);
                Assert.assertEquals(algorithms.size(), 1);
                Assert.assertEquals(algorithms.get(0), "manta");
                final int svlen = variant.getAttributeAsInt(GATKSVVCFConstants.SVLEN, -1);
                Assert.assertEquals(svlen, 301);
                Assert.assertEquals(variant.getReference(), Allele.REF_C);
                if (!omitMembers) {
                    final List<String> members = variant.getAttributeAsStringList(GATKSVVCFConstants.CLUSTER_MEMBER_IDS_KEY, null);
                    Assert.assertEquals(members.size(), 2);
                }
                final List<Allele> alts = variant.getAlternateAlleles();
                Assert.assertEquals(alts.size(), 1);
                Assert.assertEquals(alts.get(0), Allele.SV_SIMPLE_INS);
                Assert.assertEquals(variant.getStructuralVariantType(), StructuralVariantType.INS);
                for (final Genotype g : variant.getGenotypes()) {
                    if (g.getSampleName().equals("HG00096") || g.getSampleName().equals("HG00129")) {
                        Assert.assertTrue(g.isHet());
                        Assert.assertEquals(VariantContextGetters.getAttributeAsInt(g, GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, -1), 2);
                    } else {
                        Assert.assertTrue(g.isHomRef());
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
                .addVCF(getToolTestDataDir() + "1kgp_test.cnvs.vcf.gz")
                .addVCF(getToolTestDataDir() + "HG00096.manta.vcf.gz")
                .addVCF(getToolTestDataDir() + "HG00096.wham.vcf.gz")
                .addVCF(getToolTestDataDir() + "HG00129.manta.vcf.gz")
                .addVCF(getToolTestDataDir() + "HG00129.wham.vcf.gz")
                .addVCF(getToolTestDataDir() + "HG00140.manta.vcf.gz")
                .addVCF(getToolTestDataDir() + "HG00140.wham.vcf.gz")
                .add(SVCluster.PLOIDY_TABLE_LONG_NAME, getToolTestDataDir() + "1kgp.batch1.ploidy.tsv")
                .add(SVCluster.VARIANT_PREFIX_LONG_NAME, "SVx")
                .add(SVCluster.ALGORITHM_LONG_NAME, SVCluster.CLUSTER_ALGORITHM.SINGLE_LINKAGE)
                .add(StandardArgumentDefinitions.REFERENCE_LONG_NAME, REFERENCE_PATH)
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

        Assert.assertEquals(records.size(), 1338);

        // Check for one record
        int expectedRecordsFound = 0;
        for (final VariantContext variant : records) {
            Assert.assertTrue(variant.hasAttribute(GATKSVVCFConstants.CLUSTER_MEMBER_IDS_KEY));
            Assert.assertTrue(variant.hasAttribute(GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE));
            if (variant.getID().equals("SVx0000001c")) {
                expectedRecordsFound++;
                Assert.assertEquals(variant.getContig(), "chr20");
                Assert.assertEquals(variant.getStart(), 3067778);
                Assert.assertEquals(variant.getEnd(), 3067860);
                final List<String> algorithms = variant.getAttributeAsStringList(GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE, null);
                Assert.assertEquals(algorithms.size(), 2);
                Assert.assertTrue(algorithms.contains("manta"));
                Assert.assertTrue(algorithms.contains("wham"));
                final List<String> members = variant.getAttributeAsStringList(GATKSVVCFConstants.CLUSTER_MEMBER_IDS_KEY, null);
                Assert.assertEquals(members.size(), 4);
                final List<Allele> alts = variant.getAlternateAlleles();
                Assert.assertEquals(alts.size(), 1);
                Assert.assertEquals(alts.get(0), Allele.SV_SIMPLE_DUP);
                Assert.assertEquals(variant.getStructuralVariantType(), StructuralVariantType.DUP);
                for (final Genotype g : variant.getGenotypes()) {
                    if (g.getSampleName().equals("HG00096") || g.getSampleName().equals("HG00129")) {
                        Assert.assertTrue(g.isHet());
                        Assert.assertEquals(VariantContextGetters.getAttributeAsInt(g, GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, -1), 2);
                        Assert.assertEquals(VariantContextGetters.getAttributeAsInt(g, GATKSVVCFConstants.COPY_NUMBER_FORMAT, -1), 3);
                    } else {
                        Assert.assertTrue(g.isHomRef());
                    }
                }
            }
        }
        Assert.assertEquals(expectedRecordsFound, 1);
    }


    // Ensure the output buffer works correctly
    @Test
    public void testAgainstSimpleImplementation() {
        final File output = createTempFile("single_linkage_cluster", ".vcf");
        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addOutput(output)
                .add(SVCluster.PLOIDY_TABLE_LONG_NAME, getToolTestDataDir() + "1kgp.batch1.ploidy.tsv")
                .add(SVCluster.VARIANT_PREFIX_LONG_NAME, "SVx")
                .add(SVCluster.ALGORITHM_LONG_NAME, SVCluster.CLUSTER_ALGORITHM.SINGLE_LINKAGE)
                .add(JointGermlineCNVSegmentation.BREAKPOINT_SUMMARY_STRATEGY_LONG_NAME, CanonicalSVCollapser.BreakpointSummaryStrategy.MEDIAN_START_MEDIAN_END)
                .add(JointGermlineCNVSegmentation.ALT_ALLELE_SUMMARY_STRATEGY_LONG_NAME, CanonicalSVCollapser.AltAlleleSummaryStrategy.COMMON_SUBTYPE)
                .add(SVCluster.INSERTION_LENGTH_SUMMARY_STRATEGY_LONG_NAME, CanonicalSVCollapser.InsertionLengthSummaryStrategy.MEDIAN)
                .add(StandardArgumentDefinitions.REFERENCE_LONG_NAME, REFERENCE_PATH)
                .add(SVClusterEngineArgumentsCollection.DEPTH_SAMPLE_OVERLAP_FRACTION_NAME, 0)
                .add(SVClusterEngineArgumentsCollection.DEPTH_INTERVAL_OVERLAP_FRACTION_NAME, 0.5)
                .add(SVClusterEngineArgumentsCollection.DEPTH_BREAKEND_WINDOW_NAME, 2000)
                .add(SVClusterEngineArgumentsCollection.MIXED_SAMPLE_OVERLAP_FRACTION_NAME, 0)
                .add(SVClusterEngineArgumentsCollection.MIXED_INTERVAL_OVERLAP_FRACTION_NAME, 0.1)
                .add(SVClusterEngineArgumentsCollection.MIXED_BREAKEND_WINDOW_NAME, 2000)
                .add(SVClusterEngineArgumentsCollection.PESR_SAMPLE_OVERLAP_FRACTION_NAME, 0)
                .add(SVClusterEngineArgumentsCollection.PESR_INTERVAL_OVERLAP_FRACTION_NAME, 0.1)
                .add(SVClusterEngineArgumentsCollection.PESR_BREAKEND_WINDOW_NAME, 500);

        final List<String> vcfInputFilenames = Lists.newArrayList(
                "1kgp_test.cnvs.vcf.gz",
                "HG00096.manta.vcf.gz",
                "HG00096.wham.vcf.gz",
                "HG00129.manta.vcf.gz",
                "HG00129.wham.vcf.gz",
                "HG00140.manta.vcf.gz",
                "HG00140.wham.vcf.gz"
        );
        vcfInputFilenames.stream().forEach(v -> args.addVCF(getToolTestDataDir() + v));

        runCommandLine(args, SVCluster.class.getSimpleName());

        final Pair<VCFHeader, List<VariantContext>> testVcf = VariantContextTestUtils.readEntireVCFIntoMemory(output.getAbsolutePath());

        final ReferenceSequenceFile referenceSequenceFile = ReferenceUtils.createReferenceReader(new GATKPath(REFERENCE_PATH));
        final ClusteringParameters depthParameters = ClusteringParameters.createDepthParameters(0.5, 2000, 0);
        final ClusteringParameters mixedParameters = ClusteringParameters.createMixedParameters(0.1, 2000, 0);
        final ClusteringParameters pesrParameters = ClusteringParameters.createPesrParameters(0.1, 500, 0);
        final CanonicalSVClusterEngine<SVCallRecord> engine = SVClusterEngineFactory.createCanonical(
                SVClusterEngine.CLUSTERING_TYPE.SINGLE_LINKAGE,
                referenceSequenceFile.getSequenceDictionary(),
                false,
                depthParameters,
                mixedParameters,
                pesrParameters);

        vcfInputFilenames.stream()
                .flatMap(vcfFilename -> VariantContextTestUtils.readEntireVCFIntoMemory(getToolTestDataDir() + vcfFilename).getValue().stream())
                .sorted(IntervalUtils.getDictionaryOrderComparator(referenceSequenceFile.getSequenceDictionary()))
                .map(SVCallRecordUtils::create)
                .forEach(engine::add);

        final Comparator<SVCallRecord> recordComparator = SVCallRecordUtils.getCallComparator(referenceSequenceFile.getSequenceDictionary());
        final SVCollapser<SVCallRecord> collapser = SVTestUtils.defaultCollapser;
        final List<VariantContext> expectedVariants = engine.forceFlush().stream()
                .map(BasicOutputCluster::getMembers)
                .map(collapser::collapse)
                .sorted(recordComparator)
                .map(SVCallRecordUtils::getVariantBuilder)
                .map(VariantContextBuilder::make)
                .collect(Collectors.toList());
        final List<VariantContext> testVariants = testVcf.getValue();

        Assert.assertEquals(testVariants.size(), expectedVariants.size());
        for (int i = 0; i < testVariants.size(); i++) {
            final VariantContext testVariant = testVariants.get(i);
            final VariantContext expectedVariant = expectedVariants.get(i);
            Assert.assertEquals(testVariant.getContig(), expectedVariant.getContig());
            Assert.assertEquals(testVariant.getStart(), expectedVariant.getStart());
            Assert.assertEquals(testVariant.getEnd(), expectedVariant.getEnd());
            Assert.assertEquals(testVariant.getAlleles(), expectedVariant.getAlleles());
            Assert.assertEquals(testVariant.getAttributeAsStringList(GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE, "test_default"), expectedVariant.getAttributeAsStringList(GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE, "expected_default"));
            Assert.assertEquals(testVariant.getAttributeAsStringList(GATKSVVCFConstants.CLUSTER_MEMBER_IDS_KEY, "test_default"), expectedVariant.getAttributeAsStringList(GATKSVVCFConstants.CLUSTER_MEMBER_IDS_KEY, "expected_default"));
            Assert.assertEquals(testVariant.getAttributeAsString(GATKSVVCFConstants.SVTYPE, "test_default"), expectedVariant.getAttributeAsString(GATKSVVCFConstants.SVTYPE, "expected_default"));
        }
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
                .addVCF(getToolTestDataDir() + "1kgp_test.cnvs.vcf.gz")
                .addVCF(getToolTestDataDir() + "HG00096.manta.vcf.gz")
                .addVCF(getToolTestDataDir() + "HG00096.wham.vcf.gz")
                .addVCF(getToolTestDataDir() + "HG00129.manta.vcf.gz")
                .addVCF(getToolTestDataDir() + "HG00129.wham.vcf.gz")
                .addVCF(getToolTestDataDir() + "HG00140.manta.vcf.gz")
                .addVCF(getToolTestDataDir() + "HG00140.wham.vcf.gz")
                .add(SVCluster.PLOIDY_TABLE_LONG_NAME, getToolTestDataDir() + "1kgp.batch1.ploidy.tsv")
                .add(SVCluster.VARIANT_PREFIX_LONG_NAME, "SVx")
                .add(StandardArgumentDefinitions.REFERENCE_LONG_NAME, REFERENCE_PATH)
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

        //Assert.assertEquals(records.size(), 1353);

        // Check for one record
        int expectedRecordsFound = 0;
        for (final VariantContext variant : records) {
            Assert.assertTrue(variant.hasAttribute(GATKSVVCFConstants.CLUSTER_MEMBER_IDS_KEY));
            Assert.assertTrue(variant.hasAttribute(GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE));
            if (variant.getContig().equals("chr20") && variant.getStart() == 28654436) {
                expectedRecordsFound++;
                Assert.assertEquals(variant.getEnd(), 28719092);
                Assert.assertFalse(variant.hasAttribute(GATKSVVCFConstants.SVLEN));
                final List<String> algorithms = variant.getAttributeAsStringList(GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE, null);
                Assert.assertEquals(algorithms.size(), 1);
                Assert.assertTrue(algorithms.contains("manta"));
                final List<String> members = variant.getAttributeAsStringList(GATKSVVCFConstants.CLUSTER_MEMBER_IDS_KEY, null);
                Assert.assertEquals(members.size(), 2);
                final List<Allele> alts = variant.getAlternateAlleles();
                Assert.assertEquals(alts.size(), 1);
                Assert.assertEquals(alts.get(0), Allele.SV_SIMPLE_INV);
                Assert.assertEquals(variant.getStructuralVariantType(), StructuralVariantType.INV);
                final String strands = variant.getAttributeAsString(GATKSVVCFConstants.STRANDS_ATTRIBUTE, null);
                Assert.assertEquals(strands, "--");
                for (final Genotype g : variant.getGenotypes()) {
                    if (g.getSampleName().equals("HG00140") || g.getSampleName().equals("HG00129")) {
                        Assert.assertTrue(g.isHom());
                        Assert.assertEquals(VariantContextGetters.getAttributeAsInt(g, GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, -1), 2);
                    } else {
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
                .addVCF(getToolTestDataDir() + "1kgp_test.batch1.pesr.chr22.vcf.gz")
                .addVCF(getToolTestDataDir() + "1kgp_test.batch1.depth.chr22.vcf.gz")
                .add(SVCluster.PLOIDY_TABLE_LONG_NAME, getToolTestDataDir() + "1kgp.batch1.ploidy.tsv")
                .add(SVCluster.VARIANT_PREFIX_LONG_NAME, "SVx")
                .add(StandardArgumentDefinitions.REFERENCE_LONG_NAME, REFERENCE_PATH)
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

        Assert.assertEquals(header.getSampleNamesInOrder().size(), 156);

        Assert.assertEquals(records.size(), 1705);

        // Check for one record
        int expectedRecordsFound = 0;
        for (final VariantContext variant : records) {
            Assert.assertTrue(variant.hasAttribute(GATKSVVCFConstants.CLUSTER_MEMBER_IDS_KEY));
            Assert.assertTrue(variant.hasAttribute(GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE));
            if (variant.getID().equals("SVx0000041c")) {
                expectedRecordsFound++;
                Assert.assertEquals(variant.getContig(), "chr22");
                Assert.assertEquals(variant.getStart(), 38898103);
                Assert.assertEquals(variant.getEnd(), 38902679);
                final List<String> algorithms = variant.getAttributeAsStringList(GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE, null);
                Assert.assertEquals(algorithms.size(), 2);
                Assert.assertTrue(algorithms.contains(GATKSVVCFConstants.DEPTH_ALGORITHM));
                Assert.assertTrue(algorithms.contains("manta"));
                final List<String> members = variant.getAttributeAsStringList(GATKSVVCFConstants.CLUSTER_MEMBER_IDS_KEY, null);
                Assert.assertEquals(members.size(), 2);
                Assert.assertEquals(variant.getReference(), Allele.REF_A);
                final List<Allele> alts = variant.getAlternateAlleles();
                Assert.assertEquals(alts.size(), 1);
                Assert.assertEquals(alts.get(0), Allele.SV_SIMPLE_DEL);
                Assert.assertEquals(variant.getStructuralVariantType(), StructuralVariantType.DEL);
                final int nonRefGenotypeCount = (int) variant.getGenotypes().stream().filter(g -> SVCallRecordUtils.isAltGenotype(g)).count();
                Assert.assertEquals(nonRefGenotypeCount, 71);
                final int alleleCount = (int) variant.getGenotypes().stream().flatMap(g -> g.getAlleles().stream()).filter(SVCallRecordUtils::isAltAllele).count();
                Assert.assertEquals(alleleCount, 87);
                final Genotype g = variant.getGenotype("HG00129");
                Assert.assertTrue(g.isHet());
                Assert.assertEquals(VariantContextGetters.getAttributeAsInt(g, GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, -1), 2);
                Assert.assertEquals(VariantContextGetters.getAttributeAsInt(g, GATKSVVCFConstants.COPY_NUMBER_FORMAT, -1), 1);
            }
        }
        Assert.assertEquals(expectedRecordsFound, 1);
    }


    @Test
    public void testAllosome() {
        final File output = createTempFile("allosome", ".vcf");
        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addOutput(output)
                .addVCF(getToolTestDataDir() + "HG00096.manta.allo.vcf.gz")
                .addVCF(getToolTestDataDir() + "HG00129.manta.allo.vcf.gz")
                .addVCF(getToolTestDataDir() + "HG00150.manta.allo.vcf.gz")
                .add(SVCluster.PLOIDY_TABLE_LONG_NAME, getToolTestDataDir() + "1kgp.batch1.ploidy.tsv")
                .add(SVCluster.VARIANT_PREFIX_LONG_NAME, "SVx")
                .add(SVCluster.ALGORITHM_LONG_NAME, SVCluster.CLUSTER_ALGORITHM.SINGLE_LINKAGE)
                .add(StandardArgumentDefinitions.REFERENCE_LONG_NAME, REFERENCE_PATH)
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

        Assert.assertEquals(header.getSampleNamesInOrder(), Lists.newArrayList("HG00096", "HG00129", "HG00150"));

        Assert.assertEquals(records.size(), 536);

        // Check for one record
        int expectedRecordsFound = 0;
        for (final VariantContext variant : records) {
            Assert.assertTrue(variant.hasAttribute(GATKSVVCFConstants.CLUSTER_MEMBER_IDS_KEY));
            Assert.assertTrue(variant.hasAttribute(GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE));
            if (variant.getID().equals("SVx00000203")) {
                expectedRecordsFound++;
                Assert.assertEquals(variant.getContig(), "chrY");
                Assert.assertEquals(variant.getStart(), 10676436);
                Assert.assertEquals(variant.getEnd(), 10694243);
                final List<String> algorithms = variant.getAttributeAsStringList(GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE, null);
                Assert.assertEquals(algorithms.size(), 1);
                Assert.assertTrue(algorithms.contains("manta"));
                final List<String> members = variant.getAttributeAsStringList(GATKSVVCFConstants.CLUSTER_MEMBER_IDS_KEY, null);
                Assert.assertEquals(members.size(), 1);
                final List<Allele> alts = variant.getAlternateAlleles();
                Assert.assertEquals(alts.size(), 1);
                Assert.assertEquals(alts.get(0), Allele.SV_SIMPLE_DUP);
                Assert.assertEquals(variant.getStructuralVariantType(), StructuralVariantType.DUP);
                for (final Genotype g : variant.getGenotypes()) {
                    if (g.getSampleName().equals("HG00096")) {
                        Assert.assertTrue(g.isHomVar());
                        Assert.assertEquals(g.getAlleles().size(), 1);
                        Assert.assertEquals(VariantContextGetters.getAttributeAsInt(g, GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, -1), 1);
                        Assert.assertEquals(VariantContextGetters.getAttributeAsInt(g, GATKSVVCFConstants.COPY_NUMBER_FORMAT, -1), 3);
                    } else if (g.getSampleName().equals("HG00129")) {
                        Assert.assertTrue(g.isHomRef());
                        Assert.assertEquals(g.getAlleles().size(), 1);
                        Assert.assertEquals(VariantContextGetters.getAttributeAsInt(g, GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, -1), 1);
                        Assert.assertEquals(VariantContextGetters.getAttributeAsInt(g, GATKSVVCFConstants.COPY_NUMBER_FORMAT, -1), 1);
                    } else if (g.getSampleName().equals("HG00150")) {
                        Assert.assertTrue(g.isNoCall());
                        Assert.assertEquals(VariantContextGetters.getAttributeAsInt(g, GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, -1), 0);
                        Assert.assertEquals(VariantContextGetters.getAttributeAsInt(g, GATKSVVCFConstants.COPY_NUMBER_FORMAT, -1), 0);
                    }
                }
            }
        }
        Assert.assertEquals(expectedRecordsFound, 1);
    }

}