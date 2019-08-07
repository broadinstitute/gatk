package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamFiles;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import java.nio.file.Path;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.Main;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.AssemblyRegionWalker;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.testutils.CommandLineProgramTester;
import org.broadinstitute.hellbender.tools.exome.orientationbiasvariantfilter.OrientationBiasUtils;
import org.broadinstitute.hellbender.tools.walkers.annotator.StrandBiasBySample;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyBasedCallerArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.FilterMutectCalls;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.M2FiltersArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.readorientation.LearnReadOrientationModel;
import org.broadinstitute.hellbender.tools.walkers.validation.ConcordanceSummaryRecord;
import org.broadinstitute.hellbender.tools.walkers.variantutils.ValidateVariants;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.GATKProtectedVariantContextUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMFileGATKReadWriter;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

/**
 * Created by davidben on 9/1/16.
 */
@Test(groups = {"variantcalling"})
public class Mutect2IntegrationTest extends CommandLineProgramTest {
    // positions 10,000,000 - 11,000,000 of chr 20 and with most annotations removed
    private static final File GNOMAD = new File(largeFileTestDir, "very-small-gnomad.vcf");

    private static final String DREAM_BAMS_DIR = largeFileTestDir + "mutect/dream_synthetic_bams/";
    private static final File DREAM_4_NORMAL = new File(DREAM_BAMS_DIR, "normal_4.bam");
    private static final File DREAM_3_NORMAL = new File(DREAM_BAMS_DIR, "normal_3.bam");
    private static final File DREAM_4_TUMOR = new File(DREAM_BAMS_DIR, "tumor_4.bam");
    private static final File DREAM_3_TUMOR = new File(DREAM_BAMS_DIR, "tumor_3.bam");
    private static final File DREAM_2_NORMAL = new File(DREAM_BAMS_DIR, "normal_2.bam");
    private static final File DREAM_1_NORMAL = new File(DREAM_BAMS_DIR, "normal_1.bam");
    private static final File DREAM_2_TUMOR = new File(DREAM_BAMS_DIR, "tumor_2.bam");
    private static final File DREAM_1_TUMOR = new File(DREAM_BAMS_DIR, "tumor_1.bam");
    ;

    private static final String DREAM_VCFS_DIR = toolsTestDir + "mutect/dream/vcfs/";
    private static final File DREAM_4_TRUTH = new File(DREAM_VCFS_DIR, "sample_4.vcf");
    private static final File DREAM_3_TRUTH = new File(DREAM_VCFS_DIR, "sample_3.vcf");
    private static final File DREAM_2_TRUTH = new File(DREAM_VCFS_DIR, "sample_2.vcf");
    private static final File DREAM_1_TRUTH = new File(DREAM_VCFS_DIR, "sample_1.vcf");

    private static final String DREAM_MASKS_DIR = toolsTestDir + "mutect/dream/masks/";
    private static final File DREAM_4_MASK = new File(DREAM_MASKS_DIR, "mask4.list");
    private static final File DREAM_3_MASK = new File(DREAM_MASKS_DIR, "mask3.list");
    private static final File DREAM_2_MASK = new File(DREAM_MASKS_DIR, "mask2.list");
    private static final File DREAM_1_MASK = new File(DREAM_MASKS_DIR, "mask1.list");

    private static final File DREAM_4_FALSE_POSITIVES = new File(DREAM_VCFS_DIR, "sample_4.false_positives.vcf");

    private static final File NO_CONTAMINATION_TABLE = new File(toolsTestDir, "mutect/no-contamination.table");
    private static final File FIVE_PCT_CONTAMINATION_TABLE = new File(toolsTestDir, "mutect/five-pct-contamination.table");
    private static final File TEN_PCT_CONTAMINATION_TABLE = new File(toolsTestDir, "mutect/ten-pct-contamination.table");

    private static final File NA12878_MITO_BAM = new File(toolsTestDir, "mutect/mito/NA12878.bam");
    private static final File NA12878_MITO_VCF = new File(toolsTestDir, "mutect/mito/unfiltered.vcf");
    private static final File NA12878_MITO_GVCF = new File(toolsTestDir, "mitochondria/NA12878.MT.g.vcf");
    private static final File MITO_REF = new File(toolsTestDir, "mutect/mito/Homo_sapiens_assembly38.mt_only.fasta");
    private static final File DEEP_MITO_BAM = new File(largeFileTestDir, "mutect/highDPMTsnippet.bam");
    private static final String DEEP_MITO_SAMPLE_NAME = "mixture";

    private static final File FILTERING_DIR = new File(toolsTestDir, "mutect/filtering");

    private static final File GNOMAD_WITHOUT_AF_SNIPPET = new File(toolsTestDir, "mutect/gnomad-without-af.vcf");

    private static final double TLOD_MATCH_EPSILON = 0.05;
    private static final double VARIANT_TLOD_MATCH_PCT = 0.01;

    /**
     * Several DREAM challenge bams with synthetic truth data.  In order to keep file sizes manageable, bams are restricted
     * to chromosome 20, leaving ~100-200 variants, and then further restricted to 400-bp intervals centered around
     * these variants.
     * <p>
     * Because this truth data is synthetic, ideal performance is perfect concordance with the truth vcfs.
     * <p>
     * The characteristics of the data are as follows (note that we removed structural variants from the original
     * DREAM challenge vcfs):
     * <p>
     * Sample 1: pure monoclonal sample, SNVs only
     * Sample 2: 80% pure monoclonal sample, SNVs only
     * Sample 3: pure triclonal sample, subclone minor allele frequencies are 1/2, 1/3, and 1/5, SNVs and indels
     * Sample 4: 80% biclonal sample, subclone minor allele fractions are 50% and 35%, SNVs and indels
     *
     * @throws Exception
     */
    @Test(dataProvider = "dreamSyntheticData")
    public void testDreamTumorNormal(final File tumorBam, final File normalBam, final File truthVcf, final File mask,
                                     final double requiredSensitivity, final boolean tumorOnly) throws Exception {
        Utils.resetRandomGenerator();
        final File unfilteredVcf = createTempFile("unfiltered", ".vcf");
        final File filteredVcf = createTempFile("filtered", ".vcf");
        final File f1r2Counts = createTempFile("f1r2", ".tar.gz");
        final File orientationModel = createTempFile("orientation", ".tar.gz");

        final List<String> args = Arrays.asList(
                "-I", tumorBam.getAbsolutePath(),
                "-R", b37_reference_20_21,
                "-L", "20",
                "--" + M2ArgumentCollection.GERMLINE_RESOURCE_LONG_NAME, GNOMAD.getAbsolutePath(),
                "-XL", mask.getAbsolutePath(),
                "-O", unfilteredVcf.getAbsolutePath(),
                "--" + M2ArgumentCollection.F1R2_TAR_GZ_NAME, f1r2Counts.getAbsolutePath(),
                "--" + M2ArgumentCollection.DOWNSAMPLING_STRIDE_LONG_NAME, "20",
                "--max-reads-per-alignment-start", "4",
                "--" + M2ArgumentCollection.MAX_SUSPICIOUS_READS_PER_ALIGNMENT_START_LONG_NAME, "4").stream().collect(Collectors.toList());

        // tumor-only calling with gnomAD
        if (!tumorOnly) {
            final String normal = getSampleName(normalBam);
            args.addAll(Arrays.asList("-I", normalBam.getAbsolutePath(), "-" + M2ArgumentCollection.NORMAL_SAMPLE_SHORT_NAME, normal));
        }

        runCommandLine(args);

        // verify that alleles contained in likelihoods matrix but dropped from somatic calls do not show up in annotations
        // also check that alleles have been properly clipped after dropping any non-called alleles, i.e. if we had AAA AA A
        // and A got dropped, we need AAA AA -> AA A.  The condition we don't want is that all alleles share a common first base
        // and no allele has length 1.
        VariantContextTestUtils.streamVcf(unfilteredVcf)
                .forEach(vc -> {
                    for (final Genotype genotype : vc.getGenotypes()) {
                        final int[] f1r2 = OrientationBiasUtils.getF1R2(genotype);
                        Assert.assertEquals(f1r2.length, vc.getNAlleles());
                        if (vc.getAlleles().stream().filter(a -> !a.isSymbolic()).map(a -> a.getBases()[0]).distinct().count() == 1) {
                            Assert.assertTrue(vc.getAlleles().stream().anyMatch(a -> a.getBases().length == 1));
                        }
                    }
                });

        // learn orientation bias model
        new Main().instanceMain(makeCommandLineArgs(Arrays.asList("-I", f1r2Counts.getAbsolutePath(), "-O", orientationModel.getAbsolutePath()), LearnReadOrientationModel.class.getSimpleName()));

        for (final boolean runOrientationFilter : new boolean[] { true, false}) {
            // run FilterMutectCalls
            final List<String> filterArgs = new ArrayList<>();
            filterArgs.addAll(Arrays.asList("-R", b37Reference, "-V", unfilteredVcf.getAbsolutePath(), "-O", filteredVcf.getAbsolutePath()));
            if (runOrientationFilter) {
                filterArgs.addAll(Arrays.asList("--" + M2FiltersArgumentCollection.ARTIFACT_PRIOR_TABLE_NAME, orientationModel.getAbsolutePath()));
            }

            new Main().instanceMain(makeCommandLineArgs(filterArgs, FilterMutectCalls.class.getSimpleName()));


            // run Concordance
            final Path concordanceSummary = createTempPath("concordance", ".txt");
            new Main().instanceMain(makeCommandLineArgs(Arrays.asList("-truth", truthVcf.getAbsolutePath(), "-eval", filteredVcf.getAbsolutePath(), "-L", "20", "-XL", mask.getAbsolutePath(), "-summary", concordanceSummary.toAbsolutePath().toString()), "Concordance"));

            final List<ConcordanceSummaryRecord> summaryRecords = new ConcordanceSummaryRecord.Reader(concordanceSummary).toList();
            summaryRecords.forEach(rec -> {
                if (rec.getTruePositives() + rec.getFalseNegatives() > 0) {
                    Assert.assertTrue(rec.getSensitivity() > requiredSensitivity);
                    // tumor-only will have germline variants sneak in
                    if (!tumorOnly) {
                        Assert.assertTrue(rec.getPrecision() > 0.5);
                    }
                }
            });
        }
    }

    @Test
    public void testNA12878NormalNormalFiltering() throws Exception {
        Utils.resetRandomGenerator();
        final File unfilteredVcf = new File(FILTERING_DIR, "NA12878.vcf");
        final File contamination = new File(FILTERING_DIR, "contamination.table");
        final File segments = new File(FILTERING_DIR, "segments.table");
        final File stats = new File(FILTERING_DIR, "merged.stats");

        final File filteredVcf = createTempFile("filtered", ".vcf");

        // run FilterMutectCalls
        new Main().instanceMain(makeCommandLineArgs(Arrays.asList(
                "-R", b37Reference,
                "-V", unfilteredVcf.getAbsolutePath(), "-O", filteredVcf.getAbsolutePath(),
                "--" + M2FiltersArgumentCollection.TUMOR_SEGMENTATION_LONG_NAME, segments.getAbsolutePath(),
                "--" + M2FiltersArgumentCollection.CONTAMINATION_TABLE_LONG_NAME, contamination.getAbsolutePath(),
                "--" + FilterMutectCalls.FILTERING_STATS_LONG_NAME, stats.getAbsolutePath()),
                FilterMutectCalls.class.getSimpleName()));

        final long numPassVariants = VariantContextTestUtils.streamVcf(filteredVcf)
                .filter(vc -> vc.getFilters().isEmpty()).count();

        Assert.assertTrue(numPassVariants < 10);
    }

    private String getSampleName(File bam) throws IOException {
        final File nameFile = createTempFile("sample_name", ".txt");
        new Main().instanceMain(makeCommandLineArgs(Arrays.asList("-I", bam.getAbsolutePath(), "-O", nameFile.getAbsolutePath(), "-encode"), "GetSampleName"));
        return Files.readAllLines(nameFile.toPath()).get(0);
    }

    @Test(dataProvider = "twoTumorData")
    public void testTwoDreamTumorSamples(final File tumorBam1, final File tumorBam2, final Optional<File> normalBam,
                                         final File truthVcf, final File mask, final double requiredSensitivity) throws Exception {
        Utils.resetRandomGenerator();
        final File unfilteredVcf = createTempFile("unfiltered", ".vcf");
        final File filteredVcf = createTempFile("filtered", ".vcf");

        String normal = "";
        if (normalBam.isPresent()) {
            final File normalNameFile = createTempFile("normal_name", ".txt");
            new Main().instanceMain(makeCommandLineArgs(Arrays.asList("-I", normalBam.get().getAbsolutePath(), "-O", normalNameFile.getAbsolutePath(), "-encode"), "GetSampleName"));
            normal = Files.readAllLines(normalNameFile.toPath()).get(0);
        }

        final List<String> args = Arrays.asList(
                "-I", tumorBam1.getAbsolutePath(),
                "-I", tumorBam2.getAbsolutePath(),
                "-R", b37_reference_20_21,
                "-L", "20",
                "--" + M2ArgumentCollection.GERMLINE_RESOURCE_LONG_NAME, GNOMAD.getAbsolutePath(),
                "-XL", mask.getAbsolutePath(),
                "-A", "StrandBiasBySample",
                "-O", unfilteredVcf.getAbsolutePath(),
                "--" + M2ArgumentCollection.DOWNSAMPLING_STRIDE_LONG_NAME, "20",
                "--max-reads-per-alignment-start", "4",
                "--" + M2ArgumentCollection.MAX_SUSPICIOUS_READS_PER_ALIGNMENT_START_LONG_NAME, "4").stream().collect(Collectors.toList());
        ;

        if (normalBam.isPresent()) {
            args.addAll(Arrays.asList("-I", normalBam.get().getAbsolutePath(), "-" + M2ArgumentCollection.NORMAL_SAMPLE_SHORT_NAME, normal));
        }

        runCommandLine(args);

        // run FilterMutectCalls
        new Main().instanceMain(makeCommandLineArgs(Arrays.asList("-R", b37Reference, "-V", unfilteredVcf.getAbsolutePath(), "-O", filteredVcf.getAbsolutePath()), FilterMutectCalls.class.getSimpleName()));

        // run Concordance
        final Path concordanceSummary = createTempPath("concordance", ".txt");
        new Main().instanceMain(makeCommandLineArgs(Arrays.asList("-truth", truthVcf.getAbsolutePath(), "-eval", filteredVcf.getAbsolutePath(), "-L", "20", "-XL", mask.getAbsolutePath(), "-summary", concordanceSummary.toAbsolutePath().toString()), "Concordance"));

        final List<ConcordanceSummaryRecord> summaryRecords = new ConcordanceSummaryRecord.Reader(concordanceSummary).toList();
        summaryRecords.forEach(rec -> {
            if (rec.getTruePositives() + rec.getFalseNegatives() > 0) {
                Assert.assertTrue(rec.getSensitivity() > requiredSensitivity);
                // tumor-only will have germline variants sneak in
                if (normalBam.isPresent()) {
                    //Assert.assertTrue(rec.getPrecision() > 0.5);
                }
            }
        });
    }

    // make a pon with a tumor and then use this pon to call somatic variants on the same tumor
    // if the pon is doing its job all calls should be filtered by this pon
    @Test(dataProvider = "dreamSyntheticDataSample1")
    public void testPon(final File tumorBam, final File normalBam) throws Exception {
        Utils.resetRandomGenerator();
        final String normalSample = getSampleName(normalBam);

        final File ponVcf = createTempFile("pon", ".vcf");
        final String[] createPonArgs = {
                "-I", tumorBam.getAbsolutePath(),
                "-I", normalBam.getAbsolutePath(),
                "-" + M2ArgumentCollection.NORMAL_SAMPLE_SHORT_NAME, normalSample,
                "-R", b37_reference_20_21,
                "-L", "20",
                "-O", ponVcf.getAbsolutePath()
        };

        runCommandLine(createPonArgs);

        final File unfilteredVcf = createTempFile("unfiltered", ".vcf");
        final File filteredVcf = createTempFile("filtered", ".vcf");
        final String[] callWithPonArgs = {
                "-I", tumorBam.getAbsolutePath(),
                "-I", normalBam.getAbsolutePath(),
                "-" + M2ArgumentCollection.NORMAL_SAMPLE_SHORT_NAME, normalSample,
                "-" + M2ArgumentCollection.PANEL_OF_NORMALS_SHORT_NAME, ponVcf.getAbsolutePath(),
                "-R", b37_reference_20_21,
                "-L", "20",
                "-O", unfilteredVcf.getAbsolutePath()
        };

        runCommandLine(callWithPonArgs);

        // run FilterMutectCalls
        new Main().instanceMain(makeCommandLineArgs(Arrays.asList("-R", b37Reference, "-V", unfilteredVcf.getAbsolutePath(), "-O", filteredVcf.getAbsolutePath()), FilterMutectCalls.class.getSimpleName()));

        final long numVariants = VariantContextTestUtils.streamVcf(filteredVcf)
                .filter(vc -> vc.getFilters().isEmpty()).count();

        Assert.assertEquals(numVariants, 0);
    }

    // run tumor-normal mode using the original DREAM synthetic sample 1 tumor and normal restricted to
    // 1/3 of our dbSNP interval, in which there is only one true positive.
    // we want to see that the number of false positives is small
    @Test
    public void testTumorNormal() throws Exception {
        Utils.resetRandomGenerator();
        final File outputVcf = createTempFile("output", ".vcf");

        final File filteredVcf = createTempFile("filtered", ".vcf");

        final File tumorBam = new File(DREAM_BAMS_DIR, "tumor.bam");
        final File normalBam = new File(DREAM_BAMS_DIR, "normal.bam");
        final File normal2 = DREAM_2_NORMAL;
        final String normalName = getSampleName(normalBam);
        final String normal2Name = getSampleName(normal2);

        final String[] args = {
                "-I", tumorBam.getAbsolutePath(),
                "-I", normalBam.getAbsolutePath(),
                "-I", normal2.getAbsolutePath(),
                "-normal", normal2Name,
                "-" + M2ArgumentCollection.NORMAL_SAMPLE_SHORT_NAME, normalName,
                "-R", b37_reference_20_21,
                "-L", "20:10000000-10100000", // this is 1/3 of the chr 20 interval of our mini-dbSNP
                "-O", outputVcf.getAbsolutePath()
        };

        runCommandLine(args);
        new Main().instanceMain(makeCommandLineArgs(Arrays.asList("-R", b37Reference, "-V", outputVcf.getAbsolutePath(), "-O", filteredVcf.getAbsolutePath()), FilterMutectCalls.class.getSimpleName()));

        VariantContextTestUtils.streamVcf(outputVcf).flatMap(vc -> vc.getGenotypes().stream()).forEach(g -> Assert.assertTrue(g.hasAD()));
        final long numVariants = VariantContextTestUtils.streamVcf(outputVcf).count();
        Assert.assertTrue(numVariants < 4);
    }

    // run tumor-only using our mini gnomAD on NA12878, which is not a tumor
    @Test
    public void testTumorOnly() {
        Utils.resetRandomGenerator();
        final File unfilteredVcf = createTempFile("unfiltered", ".vcf");
        final File filteredVcf = createTempFile("filtered", ".vcf");

        final List<String> args = Arrays.asList("-I", NA12878_20_21_WGS_bam,
                "-R", b37_reference_20_21,
                "-L", "20:10000000-10010000",
                "-O", unfilteredVcf.getAbsolutePath(),
                "--" + M2ArgumentCollection.GERMLINE_RESOURCE_LONG_NAME, GNOMAD.getAbsolutePath());
        runCommandLine(args);

        // run FilterMutectCalls
        new Main().instanceMain(makeCommandLineArgs(Arrays.asList("-R", b37Reference, "-V", unfilteredVcf.getAbsolutePath(), "-O", filteredVcf.getAbsolutePath()), "FilterMutectCalls"));

        final long numVariantsBeforeFiltering = VariantContextTestUtils.streamVcf(unfilteredVcf).count();

        final long numVariantsPassingFilters = VariantContextTestUtils.streamVcf(filteredVcf)
                .filter(vc -> vc.getFilters().isEmpty()).count();

        // just a sanity check that this bam has some germline variants on this interval so that our test doesn't pass trivially!
        Assert.assertTrue(numVariantsBeforeFiltering > 15);

        // every variant on this interval in this sample is in gnomAD
        Assert.assertTrue(numVariantsPassingFilters < 2);
    }

    // test on an artificial bam with several contrived MNPs
    @Test
    public void testMnps() throws Exception {
        Utils.resetRandomGenerator();
        final File bam = new File(toolsTestDir, "mnp.bam");

        for (final int maxMnpDistance : new int[]{0, 1, 2, 3, 5}) {
            final File outputVcf = createTempFile("unfiltered", ".vcf");

            final List<String> args = Arrays.asList("-I", bam.getAbsolutePath(),
                    "-R", b37_reference_20_21,
                    "-L", "20:10019000-10022000",
                    "-O", outputVcf.getAbsolutePath(),
                    "-" + M2ArgumentCollection.EMISSION_LOG_SHORT_NAME, "15",
                    "-" + AssemblyBasedCallerArgumentCollection.MAX_MNP_DISTANCE_SHORT_NAME, Integer.toString(maxMnpDistance));
            runCommandLine(args);

            checkMnpOutput(maxMnpDistance, outputVcf);
        }
    }

    // this is particular to our particular artificial MNP bam -- we extract a method in order to use it for HaplotypeCaller
    private static void checkMnpOutput(int maxMnpDistance, File outputVcf) {
        // note that for testing HaplotypeCaller GVCF mode we will always have the symbolic <NON REF> allele
        final Map<Integer, List<String>> alleles = VariantContextTestUtils.streamVcf(outputVcf)
                .collect(Collectors.toMap(VariantContext::getStart, vc -> vc.getAlternateAlleles().stream().filter(a -> !a.isSymbolic()).map(Allele::getBaseString).collect(Collectors.toList())));

        // phased, two bases apart
        if (maxMnpDistance < 2) {
            Assert.assertEquals(alleles.get(10019968), Arrays.asList("G"));
            Assert.assertEquals(alleles.get(10019970), Arrays.asList("G"));
        } else {
            Assert.assertEquals(alleles.get(10019968), Arrays.asList("GAG"));
            Assert.assertTrue(!alleles.containsKey(10019970));
        }

        // adjacent and out of phase
        Assert.assertEquals(alleles.get(10020229), Arrays.asList("A"));
        Assert.assertEquals(alleles.get(10020230), Arrays.asList("G"));

        // 4-substitution MNP w/ spacings 2, 3, 4
        if (maxMnpDistance < 2) {
            Assert.assertEquals(alleles.get(10020430), Arrays.asList("G"));
            Assert.assertEquals(alleles.get(10020432), Arrays.asList("G"));
            Assert.assertEquals(alleles.get(10020435), Arrays.asList("G"));
            Assert.assertEquals(alleles.get(10020439), Arrays.asList("G"));
        } else if (maxMnpDistance < 3) {
            Assert.assertEquals(alleles.get(10020430), Arrays.asList("GAG"));
            Assert.assertEquals(alleles.get(10020435), Arrays.asList("G"));
            Assert.assertEquals(alleles.get(10020439), Arrays.asList("G"));
        } else if (maxMnpDistance < 4) {
            Assert.assertEquals(alleles.get(10020430), Arrays.asList("GAGTTG"));
            Assert.assertEquals(alleles.get(10020439), Arrays.asList("G"));
        } else {
            Assert.assertEquals(alleles.get(10020430), Arrays.asList("GAGTTGTCTG"));
        }

        // two out of phase DNPs that overlap and have a base in common
        if (maxMnpDistance > 0) {
            Assert.assertEquals(alleles.get(10020680), Arrays.asList("TA"));
            Assert.assertEquals(alleles.get(10020681), Arrays.asList("AT"));
        }
    }

    @Test
    public void testForceCalling() throws Exception {
        Utils.resetRandomGenerator();
        final File unfilteredVcf = createTempFile("unfiltered", ".vcf");

        final File givenAllelesVcf = new File(toolsTestDir, "mutect/gga_mode.vcf");
        final List<String> args = Arrays.asList("-I", NA12878_20_21_WGS_bam,
                "-R", b37_reference_20_21,
                "-L", "20:9998500-10010000",
                "-O", unfilteredVcf.getAbsolutePath(),
                "--alleles", givenAllelesVcf.getAbsolutePath());
        runCommandLine(args);

        final Map<Integer, List<Allele>> altAllelesByPosition = VariantContextTestUtils.streamVcf(unfilteredVcf)
                .collect(Collectors.toMap(vc -> vc.getStart(), vc -> vc.getAlternateAlleles()));
        for (final VariantContext vc : new FeatureDataSource<VariantContext>(givenAllelesVcf)) {
            final List<Allele> altAllelesAtThisLocus = altAllelesByPosition.get(vc.getStart());
            vc.getAlternateAlleles().forEach(a -> Assert.assertTrue(altAllelesAtThisLocus.contains(a)));
        }
    }

    /**
     * Here we give Mutect2 ridiculous kmer settings in order to force assembly to fail.
     * @throws Exception
     */
    @Test
    public void testGivenAllelesModeWithCycles() throws Exception {
        Utils.resetRandomGenerator();
        final File unfilteredVcf = createTempFile("unfiltered", ".vcf");

        final File givenAllelesVcf = new File(toolsTestDir, "mutect/gga_mode.vcf");
        final List<String> args = Arrays.asList("-I", NA12878_20_21_WGS_bam,
                "-R", b37_reference_20_21,
                "-L", "20:9998500-10010000",
                "-O", unfilteredVcf.getAbsolutePath(),
                "--alleles", givenAllelesVcf.getAbsolutePath(),
                "--kmer-size", "1",
                "--dont-increase-kmer-sizes-for-cycles");
        runCommandLine(args);

        final Map<Integer, List<Allele>> altAllelesByPosition = VariantContextTestUtils.streamVcf(unfilteredVcf)
                .collect(Collectors.toMap(vc -> vc.getStart(), vc -> vc.getAlternateAlleles()));
        for (final VariantContext vc : new FeatureDataSource<VariantContext>(givenAllelesVcf)) {
            final List<Allele> altAllelesAtThisLocus = altAllelesByPosition.get(vc.getStart());
            vc.getAlternateAlleles().forEach(a -> Assert.assertTrue(altAllelesAtThisLocus.contains(a)));
        }
    }

    // make sure that GGA mode with given alleles that normally wouldn't be called due to complete lack of coverage
    // doesn't run into any edge case bug involving empty likelihoods matrices
    @Test
    public void testGivenAllelesZeroCoverage() throws Exception {
        Utils.resetRandomGenerator();
        final File bam = DREAM_3_TUMOR;
        final File unfilteredVcf = createTempFile("unfiltered", ".vcf");
        final File givenAllelesVcf = new File(toolsTestDir, "mutect/gga_mode_2.vcf");
        final List<String> args = Arrays.asList("-I", bam.getAbsolutePath(),
                "-R", b37_reference_20_21,
                "-L", "20:1119000-1120000",
                "-O", unfilteredVcf.getAbsolutePath(),
                "--alleles", givenAllelesVcf.getAbsolutePath());
        runCommandLine(args);
    }

    // make sure we have fixed a bug where germline resources with AF=. throw errors
    @Test
    public void testMissingAF() {
        final File bam = DREAM_4_TUMOR;
        final String sample = "synthetic.challenge.set4.tumour";
        final File unfilteredVcf = createTempFile("unfiltered", ".vcf");
        final List<String> args = Arrays.asList("-I", bam.getAbsolutePath(),
                "-R", b37_reference_20_21,
                "--" + M2ArgumentCollection.GERMLINE_RESOURCE_LONG_NAME, GNOMAD_WITHOUT_AF_SNIPPET.getAbsolutePath(),
                "-L", "20:10086110",
                "-L", "20:10837425-10837426",
                "-O", unfilteredVcf.getAbsolutePath());
        runCommandLine(args);
    }

    @Test
    public void testContaminationFilter() throws Exception {
        Utils.resetRandomGenerator();
        final File unfilteredVcf = createTempFile("unfiltered", ".vcf");
        final File filteredVcfNoContamination = createTempFile("filtered-zero", ".vcf");
        final File filteredVcfFivePctContamination = createTempFile("filtered-five", ".vcf");
        final File filteredVcfTenPctContamination = createTempFile("filtered-ten", ".vcf");

        final String[] args = {
                "-I", NA12878_20_21_WGS_bam,
                "-R", b37_reference_20_21,
                "-L", "20:10000000-20010000",
                "--" + M2ArgumentCollection.GERMLINE_RESOURCE_LONG_NAME, GNOMAD.getAbsolutePath(),
                "-O", unfilteredVcf.getAbsolutePath()
        };

        runCommandLine(args);

        // run FilterMutectCalls
        new Main().instanceMain(makeCommandLineArgs(Arrays.asList("-R", b37Reference, "-V", unfilteredVcf.getAbsolutePath(), "-O", filteredVcfNoContamination.getAbsolutePath(), "--contamination-table", NO_CONTAMINATION_TABLE.getAbsolutePath()), FilterMutectCalls.class.getSimpleName()));
        new Main().instanceMain(makeCommandLineArgs(Arrays.asList("-R", b37Reference, "-V", unfilteredVcf.getAbsolutePath(), "-O", filteredVcfFivePctContamination.getAbsolutePath(), "--contamination-table", FIVE_PCT_CONTAMINATION_TABLE.getAbsolutePath()), FilterMutectCalls.class.getSimpleName()));
        new Main().instanceMain(makeCommandLineArgs(Arrays.asList("-R", b37Reference, "-V", unfilteredVcf.getAbsolutePath(), "-O", filteredVcfTenPctContamination.getAbsolutePath(), "--contamination-table", TEN_PCT_CONTAMINATION_TABLE.getAbsolutePath()), FilterMutectCalls.class.getSimpleName()));


        final Set<VariantContext> variantsFilteredAtZeroPct =
                VariantContextTestUtils.streamVcf(filteredVcfNoContamination)
                        .filter(vc -> vc.getFilters().contains(GATKVCFConstants.CONTAMINATION_FILTER_NAME))
                        .collect(Collectors.toSet());

        final Set<VariantContext> variantsFilteredAtFivePct =
                VariantContextTestUtils.streamVcf(filteredVcfFivePctContamination)
                        .filter(vc -> vc.getFilters().contains(GATKVCFConstants.CONTAMINATION_FILTER_NAME))
                        .collect(Collectors.toSet());

        final Set<VariantContext> variantsFilteredAtTenPct =
                VariantContextTestUtils.streamVcf(filteredVcfTenPctContamination)
                        .filter(vc -> vc.getFilters().contains(GATKVCFConstants.CONTAMINATION_FILTER_NAME))
                        .collect(Collectors.toSet());

        Assert.assertTrue(variantsFilteredAtZeroPct.isEmpty());
        Assert.assertTrue(variantsFilteredAtFivePct.size() < variantsFilteredAtTenPct.size());

        final List<VariantContext> missedObviousVariantsAtTenPercent = VariantContextTestUtils.streamVcf(filteredVcfTenPctContamination)
                .filter(vc -> !vc.getFilters().contains(GATKVCFConstants.CONTAMINATION_FILTER_NAME))
                .filter(VariantContext::isBiallelic)
                .filter(vc -> {
                    final int[] AD = vc.getGenotype(0).getAD();
                    return AD[1] < 0.15 * AD[0];
                }).collect(Collectors.toList());

        Assert.assertTrue(missedObviousVariantsAtTenPercent.isEmpty());

        // If the filter is smart, it won't filter variants with allele fraction much higher than the contamination
        final List<VariantContext> highAlleleFractionFilteredVariantsAtFivePercent = VariantContextTestUtils.streamVcf(filteredVcfFivePctContamination)
                .filter(vc -> vc.getFilters().contains(GATKVCFConstants.CONTAMINATION_FILTER_NAME))
                .filter(VariantContext::isBiallelic)
                .filter(vc -> {
                    final int[] AD = vc.getGenotype(0).getAD();
                    return MathUtils.sum(AD) > 30 && AD[1] > AD[0];
                }).collect(Collectors.toList());

        Assert.assertTrue(highAlleleFractionFilteredVariantsAtFivePercent.isEmpty());
    }

    // test that ReadFilterLibrary.NON_ZERO_REFERENCE_LENGTH_ALIGNMENT removes reads that consume zero reference bases
    // e.g. read name HAVCYADXX150109:1:2102:20528:2129 with cigar 23S53I
    @Test
    public void testReadsThatConsumeZeroReferenceReads() throws Exception {
        final String CONSUMES_ZERO_REFERENCE_BASES = publicTestDir + "org/broadinstitute/hellbender/tools/mutect/na12878-chr20-consumes-zero-reference-bases.bam";
        final File outputVcf = createTempFile("output", ".vcf");
        final String[] args = {
                "-I", CONSUMES_ZERO_REFERENCE_BASES,
                "-R", b37_reference_20_21,
                "-O", outputVcf.getAbsolutePath()
        };
        runCommandLine(args);
    }

    // make sure that unpaired reads that pass filtering do not cause errors
    // in particular, the read HAVCYADXX150109:1:1109:11610:46575 with SAM flag 16 fails without the patch
    @Test
    public void testUnpairedReads() throws Exception {
        final String bamWithUnpairedReads = toolsTestDir + "unpaired.bam";
        final File outputVcf = createTempFile("output", ".vcf");
        final String[] args = {
                "-I", bamWithUnpairedReads,
                "-R", b37_reference_20_21,
                "-O", outputVcf.getAbsolutePath()
        };
        runCommandLine(args);
    }

    // some bams from external pipelines use faulty adapter trimming programs that introduce identical repeated reads
    // into bams.  Although these bams fail the Picard tool ValidateSamFile, one can run HaplotypeCaller and Mutect on them
    // and get fine results.  This test ensures that this remains the case.  The test bam is a small chunk of reads surrounding
    // a germline SNP in NA12878, where we have duplicated 40 of the reads. (In practice bams of this nature might have one bad read
    // per megabase).
    @Test
    public void testBamWithRepeatedReads() {
        doMutect2Test(
                toolsTestDir + "mutect/repeated_reads.bam",
                "SM-612V3",
                "20:10018000-10020000",
                false,
                false,
                false
        );
    }

    /*
     * Test that the min_base_quality_score parameter works
     */
    @Test
    public void testMinBaseQualityScore() throws Exception {
        Utils.resetRandomGenerator();

        final File tumor = DREAM_1_TUMOR;

        final File outputAtLowThreshold = createTempFile("output", ".vcf");
        final File outputAtHighThreshold = createTempFile("output", ".vcf");

        final String[] lowThresholdArgs = {
                "-I", tumor.getAbsolutePath(),
                "-R", b37_reference_20_21,
                "-L", "20:10000000-13000000",
                "-O", outputAtLowThreshold.getAbsolutePath(),
                "--" + AssemblyBasedCallerArgumentCollection.MIN_BASE_QUALITY_SCORE_LONG_NAME, "20"
        };

        runCommandLine(lowThresholdArgs);

        final String[] highThresholdArgs = {
                "-I", tumor.getAbsolutePath(),
                "-R", b37_reference_20_21,
                "-L", "20:10000000-13000000",
                "-O", outputAtHighThreshold.getAbsolutePath(),
                "--" + AssemblyBasedCallerArgumentCollection.MIN_BASE_QUALITY_SCORE_LONG_NAME, "30"
        };

        runCommandLine(highThresholdArgs);

        try (final FeatureDataSource<VariantContext> lowThresholdSource = new FeatureDataSource<>(outputAtLowThreshold);
             final FeatureDataSource<VariantContext> highThresholdSource = new FeatureDataSource<>(outputAtHighThreshold)) {
            final List<VariantContext> variantsWithLowThreshold =
                    StreamSupport.stream(lowThresholdSource.spliterator(), false).collect(Collectors.toList());


            final List<VariantContext> variantsWithHighThreshold =
                    StreamSupport.stream(highThresholdSource.spliterator(), false).collect(Collectors.toList());

            final Set<Integer> lowStarts = variantsWithLowThreshold.stream().map(VariantContext::getStart).collect(Collectors.toSet());
            final Set<Integer> highStarts = variantsWithHighThreshold.stream().map(VariantContext::getStart).collect(Collectors.toSet());
            lowStarts.removeAll(highStarts);
            final List<Integer> diff = lowStarts.stream().sorted().collect(Collectors.toList());
            Assert.assertEquals(diff, Arrays.asList(11000090, 11000515, 12753594));
        }
    }

    // basic test on a small chunk of NA12878 mitochondria.  This is not a validation, but rather a sanity check
    // that M2 makes obvious calls, doesn't trip up on the beginning of the circular chromosome, and can handle high depth
    @Test
    public void testMitochondria() throws Exception {
        Utils.resetRandomGenerator();
        final File unfilteredVcf = createTempFile("unfiltered", ".vcf");

        final List<String> args = Arrays.asList("-I", NA12878_MITO_BAM.getAbsolutePath(),
                "-R", MITO_REF.getAbsolutePath(),
                "-L", "chrM:1-1000",
                "--" + M2ArgumentCollection.MITOCHONDRIA_MODE_LONG_NAME,
                "-O", unfilteredVcf.getAbsolutePath());
        runCommandLine(args);


        final List<VariantContext> variants = VariantContextTestUtils.streamVcf(unfilteredVcf).collect(Collectors.toList());
        final List<String> variantKeys = variants.stream().map(vc -> keyForVariant(vc)).collect(Collectors.toList());

        final List<String> expectedKeys = Arrays.asList(
                "chrM:152-152 T*, [C]",
                "chrM:263-263 A*, [G]",
                "chrM:302-302 A*, [AC, ACC, C]",
                "chrM:310-310 T*, [C, TC]",
                "chrM:750-750 A*, [G]");
        Assert.assertTrue(expectedKeys.stream().allMatch(variantKeys::contains));

        Assert.assertEquals(variants.get(0).getAttributeAsInt(GATKVCFConstants.ORIGINAL_CONTIG_MISMATCH_KEY, 0), 1671);
    }

    @Test(dataProvider = "vcfsForFiltering")
    public void testFilterMitochondria(String unfiltered, List<String> filterList, String[] extraArgs) throws Exception {
        final File filteredVcf = createTempFile("filtered", ".vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addArgument("R", MITO_REF.getAbsolutePath());
        args.addArgument("V", unfiltered);
        args.addArgument("O", filteredVcf.getPath());
        args.addBooleanArgument(M2ArgumentCollection.MITOCHONDRIA_MODE_LONG_NAME, true);
        args.addBooleanArgument(StandardArgumentDefinitions.DISABLE_SEQUENCE_DICT_VALIDATION_NAME, true);   // vcf sequence dicts don't match ref
        Arrays.stream(extraArgs).forEach(args::add);

        new Main().instanceMain(makeCommandLineArgs(args.getArgsList(), FilterMutectCalls.class.getSimpleName()));

        final List<VariantContext> variants = VariantContextTestUtils.streamVcf(filteredVcf).collect(Collectors.toList());
        final Iterator<String> expectedFilters = filterList.iterator();

        for (VariantContext v : variants) {
            final List<String> sortedFilters = new ArrayList<>(v.getFilters());
            Collections.sort(sortedFilters);
            Assert.assertEquals(sortedFilters.toString(), expectedFilters.next(), "filters don't match expected");
        }
    }

    @Test
    public void testMitochondrialRefConf() throws Exception {
        Utils.resetRandomGenerator();
        final File standardVcf = createTempFile("standard", ".vcf");
        final File unthresholded = createTempFile("unthresholded", ".vcf");


        final List<String> args = Arrays.asList("-I", NA12878_MITO_BAM.getAbsolutePath(),
                "-" + M2ArgumentCollection.TUMOR_SAMPLE_SHORT_NAME, "NA12878",
                "-R", MITO_REF.getAbsolutePath(),
                "-L", "chrM:1-1000",
                "--" + M2ArgumentCollection.MITOCHONDRIA_MODE_LONG_NAME,
                "-O", standardVcf.getAbsolutePath(),
                "-ERC", "GVCF",
                "-LODB", "-2.0",
                "-LODB", "0.0",
                "-min-AF", "0.01");
        runCommandLine(args);

        //check ref conf-specific headers are output
        final Pair<VCFHeader, List<VariantContext>> result = VariantContextTestUtils.readEntireVCFIntoMemory(standardVcf.getAbsolutePath());
        Assert.assertTrue(result.getLeft().hasFormatLine(GATKVCFConstants.TUMOR_LOG_10_ODDS_KEY));
        Assert.assertTrue(result.getLeft().getMetaDataLine(GATKVCFConstants.SYMBOLIC_ALLELE_DEFINITION_HEADER_TAG) != null);

        final List<VariantContext> variants = result.getRight();
        final Map<String, VariantContext> variantMap = variants.stream().collect(Collectors.toMap(vc -> keyForVariant(vc), Function.identity()));
        final List<String> variantKeys = new ArrayList<>(variantMap.keySet());

        final List<String> expectedKeys = Arrays.asList(
                "chrM:152-152 T*, [<NON_REF>, C]",
                "chrM:263-263 A*, [<NON_REF>, G]",
                "chrM:297-297 A*, [<NON_REF>, AC, C]",  //alt alleles get sorted when converted to keys
                //"chrM:301-301 A*, [<NON_REF>, AC, ACC]",
                //"chrM:302-302 A*, [<NON_REF>, AC, ACC, C]",  //one of these commented out variants has an allele that only appears in debug mode
                "chrM:310-310 T*, [<NON_REF>, C, TC]",
                "chrM:750-750 A*, [<NON_REF>, G]");
        Assert.assertTrue(expectedKeys.stream().allMatch(variantKeys::contains));
        //First entry should be a homRef block
        Assert.assertTrue(variantKeys.get(0).contains("*, [<NON_REF>]"));

        final CommandLineProgramTester validator = ValidateVariants.class::getSimpleName;
        final ArgumentsBuilder args2 = new ArgumentsBuilder();
        args2.addArgument("R", MITO_REF.getAbsolutePath());
        args2.addArgument("V", standardVcf.getAbsolutePath());
        args2.addArgument("L", IntervalUtils.locatableToString(new SimpleInterval("chrM:1-1000")));
        args2.add("-gvcf");
        validator.runCommandLine(args2);  //will throw a UserException if GVCF isn't contiguous

        final List<String> args3 = Arrays.asList("-I", NA12878_MITO_BAM.getAbsolutePath(),
                "-" + M2ArgumentCollection.TUMOR_SAMPLE_SHORT_NAME, "NA12878",
                "-R", MITO_REF.getAbsolutePath(),
                "-L", "chrM:1-1000",
                "--" + M2ArgumentCollection.MITOCHONDRIA_MODE_LONG_NAME,
                "-O", unthresholded.getAbsolutePath(),
                "-ERC", "GVCF",
                "-LODB", "-2.0",
                "-LODB", "0.0",
                "-min-AF", "0.00");
        runCommandLine(args3);
        final Pair<VCFHeader, List<VariantContext>> result_noThreshold = VariantContextTestUtils.readEntireVCFIntoMemory(unthresholded.getAbsolutePath());

        final Map<String, VariantContext> variantMap2 = result_noThreshold.getRight().stream().collect(Collectors.toMap(vc -> keyForVariant(vc), Function.identity()));

        //TLODs for variants should not change too much for variant allele, should change significantly for non-ref
        // however, there are edge cases where this need not be true (this might indicate the need to fix our
        // LOD calculation for the NON-REF allele), so we allow one anomalous site
        final long changedRegularAlleleLodCount = expectedKeys.stream()
                .filter(key -> !onlyNonRefTlodsChange(variantMap.get(key), variantMap2.get(key)))
                .count();

        Assert.assertTrue(changedRegularAlleleLodCount <= 1);

        final List<String> expectedRefKeys = Arrays.asList(
                //ref blocks will be dependent on TLOD band values
                "chrM:218-218 A*, [<NON_REF>]",
                "chrM:264-266 C*, [<NON_REF>]",
                "chrM:479-483 A*, [<NON_REF>]",
                "chrM:488-492 T*, [<NON_REF>]");

        //ref block boundaries aren't particularly stable, so try a few and make sure we check at least one
        boolean checkedAtLeastOneRefBlock = false;
        for (final String key : expectedRefKeys) {
            final VariantContext v1 = variantMap.get(key);
            final VariantContext v2 = variantMap2.get(key);
            if (v1 == null || v2 == null) {
                continue;
            }
            Assert.assertTrue(onlyNonRefTlodsChange(v1, v2));
            checkedAtLeastOneRefBlock = true;
        }
        Assert.assertTrue(checkedAtLeastOneRefBlock);
    }

    private boolean onlyNonRefTlodsChange(final VariantContext v1, final VariantContext v2) {
        if (v1 == null) {
            return false;
        }
        if (v2 == null) {
            return false;
        }
        if (!v1.getReference().equals(v2.getReference())) {
            return false;
        }
        if (!(v1.getAlternateAlleles().size() == v2.getAlternateAlleles().size())) {
            return false;
        }

        final double[] tlods1;
        final double[] tlods2;
        //ref blocks have TLOD in format field
        if (v1.getGenotype(0).isHomRef()) {
            tlods1 = new double[]{GATKProtectedVariantContextUtils.getAttributeAsDouble(v1.getGenotype(0), GATKVCFConstants.TUMOR_LOG_10_ODDS_KEY, 0)};
            tlods2 = new double[]{GATKProtectedVariantContextUtils.getAttributeAsDouble(v2.getGenotype(0), GATKVCFConstants.TUMOR_LOG_10_ODDS_KEY, 0)};
        } else {
            tlods1 = GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(v1, GATKVCFConstants.TUMOR_LOG_10_ODDS_KEY);
            tlods2 = GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(v2, GATKVCFConstants.TUMOR_LOG_10_ODDS_KEY);

        }
        for (int i = 0; i < v1.getAlternateAlleles().size(); i++) {
            if (!v1.getAlternateAllele(i).equals(v2.getAlternateAllele(i))) {
                return false;
            }
            //we expect the AF threshold to have a significant effect on the NON_REF TLOD, but only for that allele
            if (!v1.getAlternateAllele(i).equals(Allele.NON_REF_ALLELE)) {
                if (tlods1[i] > 0) {
                    if (Math.abs(tlods1[i] - tlods2[i]) / tlods1[i] > VARIANT_TLOD_MATCH_PCT) {
                        return false;
                    }
                } else if (Math.abs(tlods1[i] - tlods2[i]) > TLOD_MATCH_EPSILON) {
                    return false;
                }
            } else {
                if (Math.abs(tlods1[i] - tlods2[i]) < TLOD_MATCH_EPSILON) {
                    return false;
                }
            }
        }
        return true;
    }

    @Test
    @SuppressWarnings("deprecation")
    public void testAFAtHighDP() throws Exception {
        Utils.resetRandomGenerator();
        final File unfilteredVcf = createTempFile("unfiltered", ".vcf");

        final List<String> args = Arrays.asList("-I", DEEP_MITO_BAM.getAbsolutePath(),
                "-R", MITO_REF.getAbsolutePath(),
                "-L", "chrM:1-1018",
                "-ip", "300",
                "-O", unfilteredVcf.getAbsolutePath());
        runCommandLine(args);

        final List<VariantContext> variants = VariantContextTestUtils.streamVcf(unfilteredVcf).collect(Collectors.toList());

        for (final VariantContext vc : variants) {
            Assert.assertTrue(vc.isBiallelic()); //I do some lazy parsing below that won't hold for multiple alternate alleles
            Genotype g = vc.getGenotype(DEEP_MITO_SAMPLE_NAME);
            Assert.assertTrue(g.hasAD());
            final int[] ADs = g.getAD();
            Assert.assertTrue(g.hasExtendedAttribute(GATKVCFConstants.ALLELE_FRACTION_KEY));
            //Assert.assertEquals(Double.parseDouble(String.valueOf(vc.getGenotype(DEEP_MITO_SAMPLE_NAME).getExtendedAttribute(GATKVCFConstants.ALLELE_FRACTION_KEY,"0"))), (double)ADs[1]/(ADs[0]+ADs[1]), 1e-6);
            Assert.assertEquals(Double.parseDouble(String.valueOf(vc.getGenotype(DEEP_MITO_SAMPLE_NAME).getAttributeAsString(GATKVCFConstants.ALLELE_FRACTION_KEY, "0"))), (double) ADs[1] / (ADs[0] + ADs[1]), 2e-3);
        }
    }

    @DataProvider(name = "bamoutVariations")
    public Object[][] bamoutVariations() {
        return new Object[][]{
                // bamout, index, md5
                {true, true, true},
                {true, true, false},
                {true, false, true},
                {true, false, false},
        };
    }

    @Test(dataProvider = "bamoutVariations")
    public void testBamoutVariations(final boolean createBamout, final boolean createBamoutIndex, final boolean createBamoutMD5) {
        // hijack repeated reads test for bamout variations testing
        doMutect2Test(
                toolsTestDir + "mutect/repeated_reads.bam",
                "SM-612V3",
                "20:10018000-10020000",
                createBamout,
                createBamoutIndex,
                createBamoutMD5
        );
    }

    @Test()
    public void testBaseQualityFilter() throws IOException {
        // Create a test sam file
        final File samFile = File.createTempFile("synthetic", ".bam");
        final SAMFileHeader samHeader = M2TestingUtils.createSamHeader();
        final SAMFileGATKReadWriter writer = M2TestingUtils.getBareBonesSamWriter(samFile, samHeader);

        final byte poorQuality = 10;
        final byte goodQuality = 30;
        final int numReads = 20;
        final List<GATKRead> refReads = M2TestingUtils.createReads(numReads, M2TestingUtils.DEFAULT_REF_BASES, samHeader, poorQuality, "ref");
        final List<GATKRead> altReads = M2TestingUtils.createReads(numReads, M2TestingUtils.DEFAULT_ALT_BASES, samHeader, goodQuality, "alt");

        refReads.forEach(writer::addRead);
        altReads.forEach(writer::addRead);
        writer.close(); // closing the writer writes to the file
        // End creating sam file

        final File unfilteredVcf = File.createTempFile("unfiltered", ".vcf");
        final List<String> args = Arrays.asList(
                "-I", samFile.getAbsolutePath(),
                "-R", hg19_chr1_1M_Reference,
                "-O", unfilteredVcf.getAbsolutePath());
        runCommandLine(args);

        final File filteredVcf = File.createTempFile("filtered", ".vcf");
        final String[] filteringArgs = makeCommandLineArgs(Arrays.asList(
                "-R", hg19_chr1_1M_Reference,
                "-V", unfilteredVcf.getAbsolutePath(),
                "-O", filteredVcf.getAbsolutePath()),
                FilterMutectCalls.class.getSimpleName());
        new Main().instanceMain(filteringArgs);

        final Optional<VariantContext> vc = VariantContextTestUtils.streamVcf(filteredVcf).findAny();
        Assert.assertTrue(vc.isPresent());
        Assert.assertEquals(vc.get().getStart(), M2TestingUtils.DEFAULT_SNP_POSITION);
        Assert.assertFalse(vc.get().getFilters().contains(GATKVCFConstants.MEDIAN_BASE_QUALITY_FILTER_NAME));
    }

    public File createSamWithOverlappingReads(final int numAltPairs, final int refDepth) throws IOException {
        final byte altQuality = 50;
        final byte refQuality = 30;
        final File samFile = File.createTempFile("liquid-biopsy", ".bam");
        final SAMFileHeader samHeader = M2TestingUtils.createSamHeader();
        final SAMFileGATKReadWriter writer = M2TestingUtils.getBareBonesSamWriter(samFile, samHeader);

        final List<GATKRead> refReads = M2TestingUtils.createReads(refDepth, M2TestingUtils.DEFAULT_REF_BASES, samHeader, refQuality, "ref");
        refReads.forEach(writer::addRead);
        for (int i = 0; i < numAltPairs; i++) {
            // Create a read pair that completely overlap each other, which is not realistic but is easy to implement
            // and captures the essence of the issue
            final List<GATKRead> overlappingPair = ArtificialReadUtils.createPair(samHeader, "alt" + i, M2TestingUtils.DEFAULT_READ_LENGTH,
                    M2TestingUtils.DEFAULT_START_POSITION, M2TestingUtils.DEFAULT_START_POSITION, true, false);
            overlappingPair.forEach(read -> {
                read.setReadGroup(M2TestingUtils.DEFAULT_READ_GROUP_NAME);
                read.setMappingQuality(60);
                read.setBases(M2TestingUtils.DEFAULT_ALT_BASES);
                read.setBaseQualities(M2TestingUtils.getUniformBQArray(altQuality, M2TestingUtils.DEFAULT_READ_LENGTH));
                writer.addRead(read);
            });
        }

        writer.close(); // closing the writer writes to the file
        return samFile;
    }


    // Test that the strand bias annotations can count the number of reads, not fragments, when requested
    @Test
    public void testReadBasedAnnotations() throws IOException {
        // Case 1: with the read correction we lose the variant - blood biopsy-like case
        final int numAltPairs = 5;
        final int depth = 100;
        final int refDepth = depth - 2 * numAltPairs;
        final File samFileWithOverlappingReads = createSamWithOverlappingReads(numAltPairs, refDepth);

        final File unfilteredVcf = File.createTempFile("unfiltered", ".vcf");
        final File bamout = File.createTempFile("realigned", ".bam");
        final String[] args = makeCommandLineArgs(Arrays.asList(
                "-R", hg19_chr1_1M_Reference,
                "-I", samFileWithOverlappingReads.getAbsolutePath(),
                "-O", unfilteredVcf.getAbsolutePath(),
                "--bamout", bamout.getAbsolutePath(),
                "--annotation", StrandBiasBySample.class.getSimpleName(),
                "--" + AssemblyRegionWalker.MAX_STARTS_LONG_NAME, String.valueOf(depth)), Mutect2.class.getSimpleName());
        new Main().instanceMain(args);

        final Optional<VariantContext> vc = VariantContextTestUtils.streamVcf(unfilteredVcf).findAny();
        Assert.assertTrue(vc.isPresent());

        // Test case 2: we lose strand artifact. Make sure to reproduce the error and so on
        final Genotype g = vc.get().getGenotype(M2TestingUtils.DEFAULT_SAMPLE_NAME);
        final int[] contingencyTable = GATKProtectedVariantContextUtils.getAttributeAsIntArray(g, GATKVCFConstants.STRAND_BIAS_BY_SAMPLE_KEY, () -> null, -1);

        final int REF_FWD_INDEX = 0;
        final int REF_REV_INDEX = 1;
        final int ALT_FWD_INDEX = 2;
        final int ALT_REV_INDEX = 3;
        Assert.assertEquals(contingencyTable[REF_FWD_INDEX], refDepth / 2);
        Assert.assertEquals(contingencyTable[REF_REV_INDEX], refDepth / 2);
        Assert.assertEquals(contingencyTable[ALT_FWD_INDEX], numAltPairs);
        Assert.assertEquals(contingencyTable[ALT_REV_INDEX], numAltPairs);

        Assert.assertFalse(vc.get().getFilters().contains(GATKVCFConstants.STRAND_ARTIFACT_FILTER_NAME));
    }


    private File createSamWithNsandStrandBias(final int numAlts, final int numNs, final int numRefs) throws IOException {
        final byte altQuality = 50;
        final byte refQuality = 30;
        final File samFile = File.createTempFile("duplex", ".bam");
        final SAMFileHeader samHeader = M2TestingUtils.createSamHeader();
        final SAMFileGATKReadWriter writer = M2TestingUtils.getBareBonesSamWriter(samFile, samHeader);

        // create some alt reads with a strand bias
        final List<GATKRead> altReads = M2TestingUtils.createReads(numAlts, M2TestingUtils.DEFAULT_ALT_BASES, samHeader, altQuality, "alt");
        altReads.forEach(read -> {
            read.setReadGroup(M2TestingUtils.DEFAULT_READ_GROUP_NAME);
            read.setMappingQuality(60);
            read.setIsReverseStrand(false);
            read.setBases(M2TestingUtils.DEFAULT_ALT_BASES);
            read.setBaseQualities(M2TestingUtils.getUniformBQArray(altQuality, M2TestingUtils.DEFAULT_READ_LENGTH));
            writer.addRead(read);
        });

        // create some reads with Ns
        final byte[] DEFAULT_N_BASES = "CATCACACTNACTAAGCACACAGAGAATAAT".getBytes();

        final List<GATKRead> NReads = M2TestingUtils.createReads(numNs, DEFAULT_N_BASES, samHeader, altQuality, "N");
        NReads.forEach(writer::addRead);

        // create some ref reads
        final List<GATKRead> refReads = M2TestingUtils.createReads(numRefs, M2TestingUtils.DEFAULT_REF_BASES, samHeader, refQuality, "ref");
        refReads.forEach(writer::addRead);

        writer.close();
        return samFile;
    }

    private void doMutect2Test(
            final String inputBam,
            final String tumorSample,
            final String interval,
            final boolean createBamout,
            final boolean createBamoutIndex,
            final boolean createBamoutMD5) {
        final File tempDir = GATKBaseTest.createTempDir("mutect2");
        final File outputVcf = new File(tempDir, "output.vcf");
        File bamoutFile = null;

        final ArgumentsBuilder argBuilder = new ArgumentsBuilder();

        argBuilder.addInput(new File(inputBam));
        argBuilder.addReference(new File(b37_reference_20_21));
        argBuilder.addArgument("tumor", tumorSample);
        argBuilder.addOutput(new File(outputVcf.getAbsolutePath()));
        argBuilder.addArgument("L", interval);
        if (createBamout) {
            bamoutFile = new File(tempDir, "bamout.bam");
            argBuilder.addArgument(AssemblyBasedCallerArgumentCollection.BAM_OUTPUT_SHORT_NAME, bamoutFile.getAbsolutePath());
        }
        argBuilder.addBooleanArgument(StandardArgumentDefinitions.CREATE_OUTPUT_BAM_INDEX_LONG_NAME, createBamoutIndex);
        argBuilder.addBooleanArgument(StandardArgumentDefinitions.CREATE_OUTPUT_BAM_MD5_LONG_NAME, createBamoutMD5);

        runCommandLine(argBuilder);

        if (createBamout && createBamoutIndex) {
            Assert.assertNotNull(SamFiles.findIndex(bamoutFile));
        }

        if (createBamout && createBamoutIndex) {
            final File expectedMD5File = new File(bamoutFile.getAbsolutePath() + ".md5");
            Assert.assertEquals(expectedMD5File.exists(), createBamoutMD5);
        }
    }

    // tumor bam, normal bam, truth vcf, required sensitivity, tumor only
    @DataProvider(name = "dreamSyntheticData")
    public Object[][] dreamSyntheticData() {
        return new Object[][]{
                {DREAM_1_TUMOR, DREAM_1_NORMAL, DREAM_1_TRUTH, DREAM_1_MASK, 0.97, false},
                {DREAM_2_TUMOR, DREAM_2_NORMAL, DREAM_2_TRUTH, DREAM_2_MASK, 0.95, false},
                {DREAM_2_TUMOR, DREAM_2_NORMAL, DREAM_2_TRUTH, DREAM_2_MASK, 0.95, true},
                {DREAM_3_TUMOR, DREAM_3_NORMAL, DREAM_3_TRUTH, DREAM_3_MASK, 0.90, false},
                {DREAM_4_TUMOR, DREAM_4_NORMAL, DREAM_4_TRUTH, DREAM_4_MASK, 0.65, false},
                {DREAM_4_TUMOR, DREAM_4_NORMAL, DREAM_4_TRUTH, DREAM_4_MASK, 0.65, true},
        };
    }

    // tumor bam, normal bam
    @DataProvider(name = "dreamSyntheticDataSample1")
    public Object[][] dreamSyntheticDataSample1() {
        return new Object[][]{
                {DREAM_1_TUMOR, DREAM_1_NORMAL}
        };
    }

    // tumorBam1, tumorBam2, Optional<File> normalBam, truthVcf, mask, requiredSensitivity
    @DataProvider(name = "twoTumorData")
    public Object[][] twoTumorData() {
        return new Object[][]{
                {DREAM_1_TUMOR, DREAM_2_TUMOR, Optional.of(DREAM_1_NORMAL), DREAM_1_TRUTH, DREAM_1_MASK, 0.97},
                {DREAM_3_TUMOR, DREAM_4_TUMOR, Optional.of(DREAM_3_NORMAL), DREAM_3_TRUTH, DREAM_3_MASK, 0.90}
        };
    }

    @DataProvider(name = "vcfsForFiltering")
    public Object[][] vcfsForFiltering() {
        return new Object[][]{
                {NA12878_MITO_VCF.getPath(), Arrays.asList(
                        "[]",
                        "[numt_chimera]",
                        "[low_allele_frac, numt_novel, weak_evidence]",
                        "[]",
                        "[]",
                        "[]"),
                        new String[]{"--min-allele-fraction .5 --autosomal-coverage 30"}},
                {NA12878_MITO_GVCF.getPath(), Arrays.asList(
                        "[]",
                        "[base_qual, weak_evidence]",
                        "[numt_novel, weak_evidence]",
                        "[]",
                        "[base_qual, contamination, low_allele_frac, map_qual, numt_novel, position, weak_evidence]"),
                        new String[]{"-L MT:1 -L MT:37 -L MT:40 -L MT:152 -L MT:157 --min-allele-fraction .0009 --autosomal-coverage .5"}
                }
        };
    }

    //TODO: bring this to HaplotypeCallerIntegrationTest
    private Pair<Double, Double> calculateConcordance(final File outputVcf, final File truthVcf) {
        final Set<String> outputKeys = VariantContextTestUtils.streamVcf(outputVcf)
                .filter(vc -> vc.getFilters().isEmpty())
                .filter(vc -> !vc.isSymbolicOrSV())
                .map(vc -> keyForVariant(vc)).collect(Collectors.toSet());
        final Set<String> truthKeys = VariantContextTestUtils.streamVcf(truthVcf)
                .filter(vc -> vc.getFilters().isEmpty())
                .filter(vc -> !vc.isSymbolicOrSV())
                .map(vc -> keyForVariant(vc)).collect(Collectors.toSet());

        final long truePositives = outputKeys.stream().filter(truthKeys::contains).count();
        final long falsePositives = outputKeys.size() - truePositives;

        final double sensitivity = (double) truePositives / truthKeys.size();
        final double fdr = (double) falsePositives / outputKeys.size();
        return new ImmutablePair<>(sensitivity, fdr);
    }

    private static String keyForVariant(final VariantContext variant) {
        return String.format("%s:%d-%d %s, %s", variant.getContig(), variant.getStart(), variant.getEnd(), variant.getReference(),
                variant.getAlternateAlleles().stream().map(Allele::getDisplayString).sorted().collect(Collectors.toList()));
    }
}
