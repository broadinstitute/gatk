package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.variant.variantcontext.Allele;
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
import org.broadinstitute.hellbender.tools.sv.stratify.SVStratificationEngineArgumentsCollection;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.List;

public class GroupedSVClusterIntegrationTest extends CommandLineProgramTest {

    @Test
    public void testClusterStratified() {
        final File output = createTempFile("single_linkage_cluster", ".vcf");

        final String clusteringConfigFile = getToolTestDataDir() + "stratified_cluster_params.tsv";
        final String stratifyConfigFile = getToolTestDataDir() + "stratified_cluster_strata.tsv";
        final String segdupFile = getToolTestDataDir() + "../SVStratify/hg38.SegDup.chr22.bed";
        final String segdupName = "SD";
        final String repeatmaskerFile = getToolTestDataDir() + "../SVStratify/hg38.RM.chr22_subsampled.bed";
        final String repeatmaskerName = "RM";

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addOutput(output)
                .addVCF(getToolTestDataDir() + "../SVStratify/bwa_melt.chr22.vcf.gz")
                .add(SVCluster.PLOIDY_TABLE_LONG_NAME, getToolTestDataDir() + "../SVCluster/1kgp.batch1.ploidy.tsv")
                .add(SVCluster.VARIANT_PREFIX_LONG_NAME, "SVx")
                .add(SVCluster.ALGORITHM_LONG_NAME, SVCluster.CLUSTER_ALGORITHM.SINGLE_LINKAGE)
                .add(GroupedSVCluster.CLUSTERING_CONFIG_FILE_LONG_NAME, clusteringConfigFile)
                .add(SVStratificationEngineArgumentsCollection.STRATIFY_CONFIG_FILE_LONG_NAME, stratifyConfigFile)
                .add(SVStratificationEngineArgumentsCollection.CONTEXT_NAME_FILE_LONG_NAME, segdupName)
                .add(SVStratificationEngineArgumentsCollection.CONTEXT_INTERVAL_FILE_LONG_NAME, segdupFile)
                .add(SVStratificationEngineArgumentsCollection.CONTEXT_NAME_FILE_LONG_NAME, repeatmaskerName)
                .add(SVStratificationEngineArgumentsCollection.CONTEXT_INTERVAL_FILE_LONG_NAME, repeatmaskerFile)
                .add(SVStratificationEngineArgumentsCollection.OVERLAP_FRACTION_LONG_NAME, 0.5)
                .add(StandardArgumentDefinitions.REFERENCE_LONG_NAME, GATKBaseTest.hg38Reference);

        runCommandLine(args, GroupedSVCluster.class.getSimpleName());

        final Pair<VCFHeader, List<VariantContext>> vcf = VariantContextTestUtils.readEntireVCFIntoMemory(output.getAbsolutePath());
        final List<VariantContext> records = vcf.getValue();

        Assert.assertEquals(records.size(), 1437);

        // Check for specific records
        int expectedRecordsFound = 0;
        for (final VariantContext variant : records) {
            Assert.assertTrue(variant.hasAttribute(GATKSVVCFConstants.CLUSTER_MEMBER_IDS_KEY));
            Assert.assertTrue(variant.hasAttribute(GATKSVVCFConstants.STRATUM_INFO_KEY));
            Assert.assertTrue(variant.hasAttribute(GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE));
            if (variant.getID().equals("SVx00000032")) {
                expectedRecordsFound++;
                Assert.assertEquals(variant.getContig(), "chr22");
                Assert.assertEquals(variant.getStart(), 11628747);
                Assert.assertEquals(variant.getEnd(), 11629803);
                final List<String> algorithms = variant.getAttributeAsStringList(GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE, null);
                Assert.assertEquals(algorithms.size(), 2);
                Assert.assertTrue(algorithms.contains("manta"));
                Assert.assertTrue(algorithms.contains("wham"));
                final List<String> members = variant.getAttributeAsStringList(GATKSVVCFConstants.CLUSTER_MEMBER_IDS_KEY, null);
                Assert.assertEquals(members.size(), 2);
                final List<Allele> alts = variant.getAlternateAlleles();
                Assert.assertEquals(alts.size(), 1);
                Assert.assertEquals(alts.get(0), Allele.SV_SIMPLE_DEL);
                Assert.assertEquals(variant.getStructuralVariantType(), StructuralVariantType.DEL);
                Assert.assertEquals(variant.getAttribute(GATKSVVCFConstants.STRATUM_INFO_KEY), "DEL_50_5k_SD_RM");
            } else if (variant.getID().equals("SVx00000125")) {
                expectedRecordsFound++;
                Assert.assertEquals(variant.getContig(), "chr22");
                Assert.assertEquals(variant.getStart(), 22563654);
                Assert.assertEquals(variant.getEnd(), 22567049);
                final List<String> algorithms = variant.getAttributeAsStringList(GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE, null);
                Assert.assertEquals(algorithms.size(), 1);
                Assert.assertTrue(algorithms.contains("manta"));
                final List<String> members = variant.getAttributeAsStringList(GATKSVVCFConstants.CLUSTER_MEMBER_IDS_KEY, null);
                Assert.assertEquals(members.size(), 1);
                final List<Allele> alts = variant.getAlternateAlleles();
                Assert.assertEquals(alts.size(), 1);
                Assert.assertEquals(alts.get(0), Allele.SV_SIMPLE_DEL);
                Assert.assertEquals(variant.getStructuralVariantType(), StructuralVariantType.DEL);
                Assert.assertEquals(variant.getAttribute(GATKSVVCFConstants.STRATUM_INFO_KEY), SVStratify.DEFAULT_STRATUM);
            } else if (variant.getID().equals("SVx000001dc")) {
                expectedRecordsFound++;
                Assert.assertEquals(variant.getContig(), "chr22");
                Assert.assertEquals(variant.getStart(), 26060912);
                Assert.assertEquals(variant.getEnd(), 26060989);
                final List<String> algorithms = variant.getAttributeAsStringList(GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE, null);
                Assert.assertEquals(algorithms.size(), 1);
                Assert.assertTrue(algorithms.contains("manta"));
                final List<String> members = variant.getAttributeAsStringList(GATKSVVCFConstants.CLUSTER_MEMBER_IDS_KEY, null);
                Assert.assertEquals(members.size(), 1);
                final List<Allele> alts = variant.getAlternateAlleles();
                Assert.assertEquals(alts.size(), 1);
                Assert.assertEquals(alts.get(0), Allele.SV_SIMPLE_DUP);
                Assert.assertEquals(variant.getStructuralVariantType(), StructuralVariantType.DUP);
                Assert.assertEquals(variant.getAttribute(GATKSVVCFConstants.STRATUM_INFO_KEY), SVStratify.DEFAULT_STRATUM);
            }
        }
        Assert.assertEquals(expectedRecordsFound, 3);
    }
}
