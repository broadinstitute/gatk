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
import org.broadinstitute.hellbender.tools.sv.cluster.CanonicalSVCollapser;
import org.broadinstitute.hellbender.tools.sv.cluster.SVClusterWalker;
import org.broadinstitute.hellbender.tools.sv.stratify.SVStratificationEngineArgumentsCollection;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
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
                .add(SVStratificationEngineArgumentsCollection.TRACK_NAME_LONG_NAME, segdupName)
                .add(SVStratificationEngineArgumentsCollection.TRACK_INTERVAL_FILE_LONG_NAME, segdupFile)
                .add(SVStratificationEngineArgumentsCollection.TRACK_NAME_LONG_NAME, repeatmaskerName)
                .add(SVStratificationEngineArgumentsCollection.TRACK_INTERVAL_FILE_LONG_NAME, repeatmaskerFile)
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

    /**
     * Data provider for the mandatory low-mem parity matrix for {@link GroupedSVCluster}. Each row is:
     * {@code [String label, ArgumentsBuilder commonArgs]}.
     *
     * <p>Cases:
     * <ol>
     *   <li>(a) representative — SINGLE_LINKAGE, default REPRESENTATIVE breakpoint strategy,
     *       default clustering config (sample overlap 0)</li>
     *   <li>(b) medianStrategy — SINGLE_LINKAGE with MEDIAN_START_MEDIAN_END breakpoint strategy</li>
     *   <li>(c) maxClique — MAX_CLIQUE algorithm, default config, default strategy</li>
     *   <li>(d) sampleOverlap — SINGLE_LINKAGE, REPRESENTATIVE strategy, clustering config with
     *       SAMPLE_OVERLAP=0.5 for both strata</li>
     *   <li>(e) intervals — SINGLE_LINKAGE, default config, restricted to chr22:1-25000000</li>
     * </ol>
     */
    @DataProvider(name = "groupedLowMemMatrixData")
    public Object[][] groupedLowMemMatrixData() {
        final String defaultClusteringConfig = getToolTestDataDir() + "stratified_cluster_params.tsv";
        final String overlapClusteringConfig = getToolTestDataDir() + "stratified_cluster_params_overlap.tsv";
        final String stratifyConfigFile = getToolTestDataDir() + "stratified_cluster_strata.tsv";
        final String segdupFile = getToolTestDataDir() + "../SVStratify/hg38.SegDup.chr22.bed";
        final String segdupName = "SD";
        final String repeatmaskerFile = getToolTestDataDir() + "../SVStratify/hg38.RM.chr22_subsampled.bed";
        final String repeatmaskerName = "RM";
        final String inputVcf = getToolTestDataDir() + "../SVStratify/bwa_melt.chr22.vcf.gz";
        final String ploidyTable = getToolTestDataDir() + "../SVCluster/1kgp.batch1.ploidy.tsv";

        // Shared base args used by all cases
        // (tracks, stratify config, VCF, ploidy, variant prefix, overlap fraction, reference)
        // Caller adds algorithm, breakpoint strategy, clustering config, and optional interval.

        // (a) representative — default REPRESENTATIVE strategy (not added explicitly), SINGLE_LINKAGE, default config
        final ArgumentsBuilder argsRepresentative = new ArgumentsBuilder()
                .addVCF(inputVcf)
                .add(SVCluster.PLOIDY_TABLE_LONG_NAME, ploidyTable)
                .add(SVCluster.VARIANT_PREFIX_LONG_NAME, "SVx")
                .add(SVCluster.ALGORITHM_LONG_NAME, SVCluster.CLUSTER_ALGORITHM.SINGLE_LINKAGE)
                .add(GroupedSVCluster.CLUSTERING_CONFIG_FILE_LONG_NAME, defaultClusteringConfig)
                .add(SVStratificationEngineArgumentsCollection.STRATIFY_CONFIG_FILE_LONG_NAME, stratifyConfigFile)
                .add(SVStratificationEngineArgumentsCollection.TRACK_NAME_LONG_NAME, segdupName)
                .add(SVStratificationEngineArgumentsCollection.TRACK_INTERVAL_FILE_LONG_NAME, segdupFile)
                .add(SVStratificationEngineArgumentsCollection.TRACK_NAME_LONG_NAME, repeatmaskerName)
                .add(SVStratificationEngineArgumentsCollection.TRACK_INTERVAL_FILE_LONG_NAME, repeatmaskerFile)
                .add(SVStratificationEngineArgumentsCollection.OVERLAP_FRACTION_LONG_NAME, 0.5)
                .add(StandardArgumentDefinitions.REFERENCE_LONG_NAME, GATKBaseTest.hg38Reference);

        // (b) medianStrategy — SINGLE_LINKAGE with MEDIAN_START_MEDIAN_END
        final ArgumentsBuilder argsMedian = new ArgumentsBuilder()
                .addVCF(inputVcf)
                .add(SVCluster.PLOIDY_TABLE_LONG_NAME, ploidyTable)
                .add(SVCluster.VARIANT_PREFIX_LONG_NAME, "SVx")
                .add(SVCluster.ALGORITHM_LONG_NAME, SVCluster.CLUSTER_ALGORITHM.SINGLE_LINKAGE)
                .add(JointGermlineCNVSegmentation.BREAKPOINT_SUMMARY_STRATEGY_LONG_NAME,
                        CanonicalSVCollapser.BreakpointSummaryStrategy.MEDIAN_START_MEDIAN_END)
                .add(GroupedSVCluster.CLUSTERING_CONFIG_FILE_LONG_NAME, defaultClusteringConfig)
                .add(SVStratificationEngineArgumentsCollection.STRATIFY_CONFIG_FILE_LONG_NAME, stratifyConfigFile)
                .add(SVStratificationEngineArgumentsCollection.TRACK_NAME_LONG_NAME, segdupName)
                .add(SVStratificationEngineArgumentsCollection.TRACK_INTERVAL_FILE_LONG_NAME, segdupFile)
                .add(SVStratificationEngineArgumentsCollection.TRACK_NAME_LONG_NAME, repeatmaskerName)
                .add(SVStratificationEngineArgumentsCollection.TRACK_INTERVAL_FILE_LONG_NAME, repeatmaskerFile)
                .add(SVStratificationEngineArgumentsCollection.OVERLAP_FRACTION_LONG_NAME, 0.5)
                .add(StandardArgumentDefinitions.REFERENCE_LONG_NAME, GATKBaseTest.hg38Reference);

        // (c) maxClique — MAX_CLIQUE algorithm, default config, default (REPRESENTATIVE) strategy
        final ArgumentsBuilder argsMaxClique = new ArgumentsBuilder()
                .addVCF(inputVcf)
                .add(SVCluster.PLOIDY_TABLE_LONG_NAME, ploidyTable)
                .add(SVCluster.VARIANT_PREFIX_LONG_NAME, "SVx")
                .add(SVCluster.ALGORITHM_LONG_NAME, SVCluster.CLUSTER_ALGORITHM.MAX_CLIQUE)
                .add(GroupedSVCluster.CLUSTERING_CONFIG_FILE_LONG_NAME, defaultClusteringConfig)
                .add(SVStratificationEngineArgumentsCollection.STRATIFY_CONFIG_FILE_LONG_NAME, stratifyConfigFile)
                .add(SVStratificationEngineArgumentsCollection.TRACK_NAME_LONG_NAME, segdupName)
                .add(SVStratificationEngineArgumentsCollection.TRACK_INTERVAL_FILE_LONG_NAME, segdupFile)
                .add(SVStratificationEngineArgumentsCollection.TRACK_NAME_LONG_NAME, repeatmaskerName)
                .add(SVStratificationEngineArgumentsCollection.TRACK_INTERVAL_FILE_LONG_NAME, repeatmaskerFile)
                .add(SVStratificationEngineArgumentsCollection.OVERLAP_FRACTION_LONG_NAME, 0.5)
                .add(StandardArgumentDefinitions.REFERENCE_LONG_NAME, GATKBaseTest.hg38Reference);

        // (d) sampleOverlap — SINGLE_LINKAGE, REPRESENTATIVE strategy, clustering config with SAMPLE_OVERLAP=0.5
        final ArgumentsBuilder argsSampleOverlap = new ArgumentsBuilder()
                .addVCF(inputVcf)
                .add(SVCluster.PLOIDY_TABLE_LONG_NAME, ploidyTable)
                .add(SVCluster.VARIANT_PREFIX_LONG_NAME, "SVx")
                .add(SVCluster.ALGORITHM_LONG_NAME, SVCluster.CLUSTER_ALGORITHM.SINGLE_LINKAGE)
                .add(GroupedSVCluster.CLUSTERING_CONFIG_FILE_LONG_NAME, overlapClusteringConfig)
                .add(SVStratificationEngineArgumentsCollection.STRATIFY_CONFIG_FILE_LONG_NAME, stratifyConfigFile)
                .add(SVStratificationEngineArgumentsCollection.TRACK_NAME_LONG_NAME, segdupName)
                .add(SVStratificationEngineArgumentsCollection.TRACK_INTERVAL_FILE_LONG_NAME, segdupFile)
                .add(SVStratificationEngineArgumentsCollection.TRACK_NAME_LONG_NAME, repeatmaskerName)
                .add(SVStratificationEngineArgumentsCollection.TRACK_INTERVAL_FILE_LONG_NAME, repeatmaskerFile)
                .add(SVStratificationEngineArgumentsCollection.OVERLAP_FRACTION_LONG_NAME, 0.5)
                .add(StandardArgumentDefinitions.REFERENCE_LONG_NAME, GATKBaseTest.hg38Reference);

        // (e) intervals — SINGLE_LINKAGE, default config, restricted to chr22:1-25000000
        final ArgumentsBuilder argsIntervals = new ArgumentsBuilder()
                .addVCF(inputVcf)
                .add(SVCluster.PLOIDY_TABLE_LONG_NAME, ploidyTable)
                .add(SVCluster.VARIANT_PREFIX_LONG_NAME, "SVx")
                .add(SVCluster.ALGORITHM_LONG_NAME, SVCluster.CLUSTER_ALGORITHM.SINGLE_LINKAGE)
                .add(StandardArgumentDefinitions.INTERVALS_LONG_NAME, "chr22:1-25000000")
                .add(GroupedSVCluster.CLUSTERING_CONFIG_FILE_LONG_NAME, defaultClusteringConfig)
                .add(SVStratificationEngineArgumentsCollection.STRATIFY_CONFIG_FILE_LONG_NAME, stratifyConfigFile)
                .add(SVStratificationEngineArgumentsCollection.TRACK_NAME_LONG_NAME, segdupName)
                .add(SVStratificationEngineArgumentsCollection.TRACK_INTERVAL_FILE_LONG_NAME, segdupFile)
                .add(SVStratificationEngineArgumentsCollection.TRACK_NAME_LONG_NAME, repeatmaskerName)
                .add(SVStratificationEngineArgumentsCollection.TRACK_INTERVAL_FILE_LONG_NAME, repeatmaskerFile)
                .add(SVStratificationEngineArgumentsCollection.OVERLAP_FRACTION_LONG_NAME, 0.5)
                .add(StandardArgumentDefinitions.REFERENCE_LONG_NAME, GATKBaseTest.hg38Reference);

        return new Object[][]{
                {"representative",  argsRepresentative},
                {"medianStrategy",  argsMedian},
                {"maxClique",       argsMaxClique},
                {"sampleOverlap",   argsSampleOverlap},
                {"intervals",       argsIntervals},
        };
    }

    /**
     * Mandatory low-mem parity matrix for {@link GroupedSVCluster}. For each scenario, runs the tool
     * single-pass (baseline) and with {@code --low-mem} and asserts byte-identical output: same record
     * count, contig, start, end, ID, all INFO attributes, filters, and all per-sample genotype alleles
     * and extended attributes.
     */
    @Test(dataProvider = "groupedLowMemMatrixData")
    public void testGroupedLowMemMatrix(final String label, final ArgumentsBuilder commonArgs) {
        // Baseline (single-pass)
        final File outputBaseline = createTempFile("grouped_lowmem_matrix_baseline_" + label, ".vcf");
        runCommandLine(commonArgs.copy().addOutput(outputBaseline), GroupedSVCluster.class.getSimpleName());

        // Low-mem pass
        final File outputLowMem = createTempFile("grouped_lowmem_matrix_lowmem_" + label, ".vcf");
        runCommandLine(commonArgs.copy().addOutput(outputLowMem)
                .add(SVClusterWalker.LOW_MEM_LONG_NAME, true), GroupedSVCluster.class.getSimpleName());

        final List<VariantContext> baselineRecords =
                VariantContextTestUtils.readEntireVCFIntoMemory(outputBaseline.getAbsolutePath()).getValue();
        final List<VariantContext> lowMemRecords =
                VariantContextTestUtils.readEntireVCFIntoMemory(outputLowMem.getAbsolutePath()).getValue();

        Assert.assertEquals(lowMemRecords.size(), baselineRecords.size(),
                label + ": --low-mem produced a different number of records than baseline");

        for (int i = 0; i < baselineRecords.size(); i++) {
            final VariantContext base = baselineRecords.get(i);
            final VariantContext lowMem = lowMemRecords.get(i);
            final String ctx = label + " record[" + i + "] id=" + base.getID();
            Assert.assertEquals(lowMem.getContig(), base.getContig(), ctx + ": contig mismatch");
            Assert.assertEquals(lowMem.getStart(), base.getStart(), ctx + ": start mismatch");
            Assert.assertEquals(lowMem.getEnd(), base.getEnd(), ctx + ": end mismatch");
            Assert.assertEquals(lowMem.getID(), base.getID(), ctx + ": ID mismatch");
            Assert.assertEquals(lowMem.getAlleles(), base.getAlleles(), ctx + ": alleles mismatch");
            Assert.assertEquals(lowMem.getAttributes(), base.getAttributes(), ctx + ": attributes mismatch");
            Assert.assertEquals(lowMem.getFilters(), base.getFilters(), ctx + ": filters mismatch");
            Assert.assertEquals(lowMem.getGenotypes().size(), base.getGenotypes().size(),
                    ctx + ": genotype count mismatch");
            for (final Genotype baseGt : base.getGenotypes()) {
                final Genotype lowMemGt = lowMem.getGenotype(baseGt.getSampleName());
                Assert.assertNotNull(lowMemGt, ctx + ": missing genotype for sample " + baseGt.getSampleName());
                Assert.assertEquals(lowMemGt.getAlleles(), baseGt.getAlleles(),
                        ctx + ": alleles mismatch for sample " + baseGt.getSampleName());
                Assert.assertEquals(lowMemGt.getExtendedAttributes(), baseGt.getExtendedAttributes(),
                        ctx + ": extended attributes mismatch for sample " + baseGt.getSampleName());
            }
        }
    }
}
