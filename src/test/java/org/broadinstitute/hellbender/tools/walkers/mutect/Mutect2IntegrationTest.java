package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamFiles;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.Main;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.tools.exome.orientationbiasvariantfilter.OrientationBiasUtils;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyBasedCallerArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.validation.ConcordanceSummaryRecord;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
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
    private static final String DREAM_VCFS_DIR = toolsTestDir + "mutect/dream/vcfs/";
    private static final String DREAM_MASKS_DIR = toolsTestDir + "mutect/dream/masks/";

    private static final File NO_CONTAMINATION_TABLE = new File(toolsTestDir, "mutect/no-contamination.table");
    private static final File FIVE_PCT_CONTAMINATION_TABLE = new File(toolsTestDir, "mutect/five-pct-contamination.table");
    private static final File TEN_PCT_CONTAMINATION_TABLE = new File(toolsTestDir, "mutect/ten-pct-contamination.table");

    private static final File NA12878_MITO_BAM = new File(toolsTestDir, "mutect/mito/NA12878.bam");
    private static final File MITO_REF = new File(toolsTestDir, "mutect/mito/Homo_sapiens_assembly38.mt_only.fasta");
    private static final File DEEP_MITO_BAM = new File(largeFileTestDir, "mutect/highDPMTsnippet.bam");
    private static final String DEEP_MITO_SAMPLE_NAME = "mixture";

    /**
     * Several DREAM challenge bams with synthetic truth data.  In order to keep file sizes manageable, bams are restricted
     * to chromosome 20, leaving ~100-200 variants, and then further restricted to 400-bp intervals centered around
     * these variants.
     *
     * Because this truth data is synthetic, ideal performance is perfect concordance with the truth vcfs.
     *
     * The characteristics of the data are as follows (note that we removed structural variants from the original
     * DREAM challenge vcfs):
     *
     * Sample 1: pure monoclonal sample, SNVs only
     * Sample 2: 80% pure monoclonal sample, SNVs only
     * Sample 3: pure triclonal sample, subclone minor allele frequencies are 1/2, 1/3, and 1/5, SNVs and indels
     * Sample 4: 80% biclonal sample, subclone minor allele fractions are 50% and 35%, SNVs and indels
     *
     * @throws Exception
     */
    @Test(dataProvider = "dreamSyntheticData")
    public void testDreamTumorNormal(final File tumorBam, final String tumorSample, final File normalBam, final String normalSample,
                                     final File truthVcf, final File mask, final double requiredSensitivity, final boolean tumorOnly) throws Exception {
        Utils.resetRandomGenerator();
        final File unfilteredVcf = createTempFile("unfiltered", ".vcf");
        final File filteredVcf = createTempFile("filtered", ".vcf");

        final File tumorNameFile = createTempFile("tumor_name", ".txt");
        final File normalNameFile = createTempFile("normal_name", ".txt");

        new Main().instanceMain(makeCommandLineArgs(Arrays.asList("-I", tumorBam.getAbsolutePath(), "-O", tumorNameFile.getAbsolutePath(), "-encode"), "GetSampleName"));
        new Main().instanceMain(makeCommandLineArgs(Arrays.asList("-I", normalBam.getAbsolutePath(), "-O", normalNameFile.getAbsolutePath(), "-encode"), "GetSampleName"));
        final String tumor = Files.readAllLines(tumorNameFile.toPath()).get(0);
        final String normal = Files.readAllLines(normalNameFile.toPath()).get(0);

        final List<String> args = Arrays.asList(
                "-I", tumorBam.getAbsolutePath(),
                "-" + M2ArgumentCollection.TUMOR_SAMPLE_SHORT_NAME, tumor,
                "-R", b37_reference_20_21,
                "-L", "20",
                "--" + M2ArgumentCollection.GERMLINE_RESOURCE_LONG_NAME, GNOMAD.getAbsolutePath(),
                "-XL", mask.getAbsolutePath(),
                "-O", unfilteredVcf.getAbsolutePath(),
                "--" + M2ArgumentCollection.DOWNSAMPLING_STRIDE_LONG_NAME, "20",
                "--max-reads-per-alignment-start", "4",
                "--" + M2ArgumentCollection.MAX_SUSPICIOUS_READS_PER_ALIGNMENT_START_LONG_NAME, "4").stream().collect(Collectors.toList());;

        // tumor-only calling with gnomAD
        if (!tumorOnly) {
            args.addAll(Arrays.asList("-I", normalBam.getAbsolutePath(), "-" + M2ArgumentCollection.NORMAL_SAMPLE_SHORT_NAME, normal));
        };

        runCommandLine(args);

        // run FilterMutectCalls
        new Main().instanceMain(makeCommandLineArgs(Arrays.asList("-V", unfilteredVcf.getAbsolutePath(), "-O", filteredVcf.getAbsolutePath()), FilterMutectCalls.class.getSimpleName()));

        // verify that alleles contained in likelihoods matrix but dropped from somatic calls do not show up in annotations
        // also check that alleles have been properly clipped after dropping any non-called alleles, i.e. if we had AAA AA A
        // and A got dropped, we need AAA AA -> AA A.  The condition we don't want is that all alleles share a common first base
        // and no allele has length 1.
        VariantContextTestUtils.streamVcf(unfilteredVcf)
                .forEach(vc -> {
                    final Genotype tumorGenotype = vc.getGenotype(tumorSample);
                    final int[] f1r2 = OrientationBiasUtils.getF1R2(tumorGenotype);
                    Assert.assertEquals(f1r2.length, vc.getNAlleles());
                    if (vc.getAlleles().stream().filter(a -> !a.isSymbolic()).map(a -> a.getBases()[0]).distinct().count() == 1) {
                        Assert.assertTrue(vc.getAlleles().stream().anyMatch(a -> a.getBases().length == 1));
                    }
                });

        // run Concordance
        final File concordanceSummary = createTempFile("concordance", ".txt");
        new Main().instanceMain(makeCommandLineArgs(Arrays.asList("-truth", truthVcf.getAbsolutePath(), "-eval", filteredVcf.getAbsolutePath(), "-L", "20", "-XL", mask.getAbsolutePath(), "-summary", concordanceSummary.getAbsolutePath()), "Concordance"));

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

    // make a pon with a tumor and then use this pon to call somatic variants on the same tumor
    // if the pon is doing its job all calls should be filtered by this pon
    @Test(dataProvider = "dreamSyntheticDataSample1")
    public void testPon(final File tumorBam, final String tumorSample, final File normalBam, final String normalSample) throws Exception {
        Utils.resetRandomGenerator();

        final File ponVcf = createTempFile("pon", ".vcf");
        final String[] createPonArgs = {
                "-I", tumorBam.getAbsolutePath(),
                "-" + M2ArgumentCollection.TUMOR_SAMPLE_SHORT_NAME, tumorSample,
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
                "-" + M2ArgumentCollection.TUMOR_SAMPLE_SHORT_NAME, tumorSample,
                "-I", normalBam.getAbsolutePath(),
                "-" + M2ArgumentCollection.NORMAL_SAMPLE_SHORT_NAME, normalSample,
                "-" + M2ArgumentCollection.PANEL_OF_NORMALS_SHORT_NAME, ponVcf.getAbsolutePath(),
                "-R", b37_reference_20_21,
                "-L", "20",
                "-O", unfilteredVcf.getAbsolutePath()
        };

        runCommandLine(callWithPonArgs);

        // run FilterMutectCalls
        new Main().instanceMain(makeCommandLineArgs(Arrays.asList("-V", unfilteredVcf.getAbsolutePath(), "-O", filteredVcf.getAbsolutePath()), FilterMutectCalls.class.getSimpleName()));

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

        final File tumorBam = new File(DREAM_BAMS_DIR, "tumor.bam");
        final String tumorName = "synthetic.challenge.set1.tumor";
        final File normalBam = new File(DREAM_BAMS_DIR, "normal.bam");
        final String normalName = "synthetic.challenge.set1.normal";

        final String[] args = {
                "-I", tumorBam.getAbsolutePath(),
                "-" + M2ArgumentCollection.TUMOR_SAMPLE_SHORT_NAME, tumorName,
                "-I", normalBam.getAbsolutePath(),
                "-" + M2ArgumentCollection.NORMAL_SAMPLE_SHORT_NAME, normalName,
                "-R", b37_reference_20_21,
                "-L", "20:10000000-10100000", // this is 1/3 of the chr 20 interval of our mini-dbSNP
                "-O", outputVcf.getAbsolutePath()
        };

        runCommandLine(args);
        VariantContextTestUtils.streamVcf(outputVcf).forEach(a -> Assert.assertTrue(a.getGenotype(tumorName).hasAD()));
        final long numVariants = VariantContextTestUtils.streamVcf(outputVcf).count();
        Assert.assertTrue(numVariants < 4);
    }

    // run tumor-only using our mini gnomAD on NA12878, which is not a tumor
    @Test
    public void testTumorOnly() throws Exception {
        Utils.resetRandomGenerator();
        final File unfilteredVcf = createTempFile("unfiltered", ".vcf");
        final File filteredVcf = createTempFile("filtered", ".vcf");

        final List<String> args = Arrays.asList("-I", NA12878_20_21_WGS_bam,
                "-" + M2ArgumentCollection.TUMOR_SAMPLE_SHORT_NAME, "NA12878",
                "-R", b37_reference_20_21,
                "-L", "20:10000000-10010000",
                "-O", unfilteredVcf.getAbsolutePath(),
                "--" + M2ArgumentCollection.GERMLINE_RESOURCE_LONG_NAME, GNOMAD.getAbsolutePath());
        runCommandLine(args);

        // run FilterMutectCalls
        new Main().instanceMain(makeCommandLineArgs(Arrays.asList("-V", unfilteredVcf.getAbsolutePath(), "-O", filteredVcf.getAbsolutePath()), "FilterMutectCalls"));

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

        for (final int maxMnpDistance : new int[] {0, 1, 2, 3, 5}) {
            final File outputVcf = createTempFile("unfiltered", ".vcf");

            final List<String> args = Arrays.asList("-I", bam.getAbsolutePath(),
                    "-" + M2ArgumentCollection.TUMOR_SAMPLE_SHORT_NAME, "NA12878",
                    "-R", b37_reference_20_21,
                    "-L", "20:10019000-10022000",
                    "-O", outputVcf.getAbsolutePath(),
                    "-" + M2ArgumentCollection.EMISSION_LOG_SHORT_NAME, "15",
                    "-" + M2ArgumentCollection.MAX_MNP_DISTANCE_SHORT_NAME, Integer.toString(maxMnpDistance));
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
    public void testGivenAllelesMode() throws Exception {
        Utils.resetRandomGenerator();
        final File unfilteredVcf = createTempFile("unfiltered", ".vcf");

        final File givenAllelesVcf = new File(toolsTestDir, "mutect/gga_mode.vcf");
        final List<String> args = Arrays.asList("-I", NA12878_20_21_WGS_bam,
                "-" + M2ArgumentCollection.TUMOR_SAMPLE_SHORT_NAME, "NA12878",
                "-R", b37_reference_20_21,
                "-L", "20:10000000-10010000",
                "-O", unfilteredVcf.getAbsolutePath(),
                "--genotyping-mode", "GENOTYPE_GIVEN_ALLELES",
                "--alleles", givenAllelesVcf.getAbsolutePath());
        runCommandLine(args);

        final Map<Integer, List<Allele>> altAllelesByPosition = VariantContextTestUtils.streamVcf(unfilteredVcf)
                .collect(Collectors.toMap(vc -> vc.getStart(), vc-> vc.getAlternateAlleles()));
        for (final VariantContext vc : new FeatureDataSource<VariantContext>(givenAllelesVcf)) {
            final List<Allele> altAllelesAtThisLocus = altAllelesByPosition.get(vc.getStart());
            vc.getAlternateAlleles().forEach(a -> Assert.assertTrue(altAllelesAtThisLocus.contains(a)));
        }
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
                "-" + M2ArgumentCollection.TUMOR_SAMPLE_SHORT_NAME, "NA12878",
                "-R", b37_reference_20_21,
                "-L", "20:10000000-20010000",
                "--" + M2ArgumentCollection.GERMLINE_RESOURCE_LONG_NAME, GNOMAD.getAbsolutePath(),
                "-O", unfilteredVcf.getAbsolutePath()
        };

        runCommandLine(args);

        // run FilterMutectCalls
        new Main().instanceMain(makeCommandLineArgs(Arrays.asList("-V", unfilteredVcf.getAbsolutePath(), "-O", filteredVcfNoContamination.getAbsolutePath(), "--contamination-table", NO_CONTAMINATION_TABLE.getAbsolutePath()), FilterMutectCalls.class.getSimpleName()));
        new Main().instanceMain(makeCommandLineArgs(Arrays.asList("-V", unfilteredVcf.getAbsolutePath(), "-O", filteredVcfFivePctContamination.getAbsolutePath(), "--contamination-table", FIVE_PCT_CONTAMINATION_TABLE.getAbsolutePath()), FilterMutectCalls.class.getSimpleName()));
        new Main().instanceMain(makeCommandLineArgs(Arrays.asList("-V", unfilteredVcf.getAbsolutePath(), "-O", filteredVcfTenPctContamination.getAbsolutePath(), "--contamination-table", TEN_PCT_CONTAMINATION_TABLE.getAbsolutePath()), FilterMutectCalls.class.getSimpleName()));


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
                    return AD[1] < 0.2 * AD[0];
                }).collect(Collectors.toList());

        Assert.assertTrue(missedObviousVariantsAtTenPercent.isEmpty());

        // If the filter is smart, it won't filter variants with allele fraction much higher than the contamination
        final List<VariantContext> highAlleleFractionFilteredVariantsAtFivePercent = VariantContextTestUtils.streamVcf(filteredVcfFivePctContamination)
                .filter(vc -> vc.getFilters().contains(GATKVCFConstants.CONTAMINATION_FILTER_NAME))
                .filter(VariantContext::isBiallelic)
                .filter(vc -> {
                    final int[] AD = vc.getGenotype(0).getAD();
                    return MathUtils.sum(AD) > 20 && AD[1] > AD[0];
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
                "-" + M2ArgumentCollection.TUMOR_SAMPLE_SHORT_NAME, "SM-612V3",
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
                "-" + M2ArgumentCollection.TUMOR_SAMPLE_SHORT_NAME, "SM-612V3",
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

        final File tumor = new File(DREAM_BAMS_DIR, "tumor_1.bam");
        final String tumorSample = "tumor sample";

        final File outputAtLowThreshold = createTempFile("output", ".vcf");
        final File outputAtHighThreshold = createTempFile("output", ".vcf");

        final String[] lowThresholdArgs = {
                "-I", tumor.getAbsolutePath(),
                "-" + M2ArgumentCollection.TUMOR_SAMPLE_SHORT_NAME, tumorSample,
                "-R", b37_reference_20_21,
                "-L", "20:10000000-13000000",
                "-O", outputAtLowThreshold.getAbsolutePath(),
                "--" + AssemblyBasedCallerArgumentCollection.MIN_BASE_QUALITY_SCORE_LONG_NAME, "20"
        };

        runCommandLine(lowThresholdArgs);

        final String[] highThresholdArgs = {
                "-I", tumor.getAbsolutePath(),
                "-" + M2ArgumentCollection.TUMOR_SAMPLE_SHORT_NAME, tumorSample,
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
                "-" + M2ArgumentCollection.TUMOR_SAMPLE_SHORT_NAME, "NA12878",
                "-R", MITO_REF.getAbsolutePath(),
                "-L", "chrM:1-1000",
                "-min-pruning", "5",
                "-O", unfilteredVcf.getAbsolutePath());
        runCommandLine(args);


        final List<VariantContext> variants = VariantContextTestUtils.streamVcf(unfilteredVcf).collect(Collectors.toList());
        final Set<String> variantKeys = variants.stream().map(vc -> keyForVariant(vc)).collect(Collectors.toSet());

        final List<String> expectedKeys = Arrays.asList(
                "chrM:152-152 [T*, C]",
                "chrM:263-263 [A*, G]",
                "chrM:301-301 [A*, AC]",
                "chrM:302-302 [A*, AC, C, ACC]",
                "chrM:310-310 [T*, TC]",
                "chrM:750-750 [A*, G]");
        Assert.assertTrue(expectedKeys.stream().allMatch(variantKeys::contains));
    }

   @Test
   @SuppressWarnings("deprecation")
   public void testAFfromADoverHighDP() throws Exception {
        Utils.resetRandomGenerator();
        final File unfilteredVcf = createTempFile("unfiltered", ".vcf");

        final List<String> args = Arrays.asList("-I", DEEP_MITO_BAM.getAbsolutePath(),
                "-" + M2ArgumentCollection.TUMOR_SAMPLE_SHORT_NAME, DEEP_MITO_SAMPLE_NAME,
                "-R", MITO_REF.getAbsolutePath(),
                "-L", "chrM:1-1018",
                "-ip", "300",
                "-min-pruning", "4",
                "--" + M2ArgumentCollection.GET_AF_FROM_AD_LONG_NAME,
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
            Assert.assertEquals(Double.parseDouble(String.valueOf(vc.getGenotype(DEEP_MITO_SAMPLE_NAME).getAttributeAsString(GATKVCFConstants.ALLELE_FRACTION_KEY,"0"))), (double)ADs[1]/(ADs[0]+ADs[1]), 1e-6);
        }
    }

    @DataProvider(name="bamoutVariations")
    public Object[][] bamoutVariations() {
        return new Object[][]{
                // bamout, index, md5
                { true, true, true },
                { true, true, false },
                { true, false, true },
                { true, false, false },
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
        final List<GATKRead> refReads = M2TestingUtils.createReads(numReads, M2TestingUtils.DEFAULT_REF_BASES, samHeader, poorQuality);
        final List<GATKRead> alt1Reads = M2TestingUtils.createReads(numReads, M2TestingUtils.DEFAULT_ALT_BASES, samHeader, goodQuality);

        refReads.forEach(writer::addRead);
        alt1Reads.forEach(writer::addRead);
        writer.close(); // closing the writer writes to the file
        // End creating sam file

        final File unfilteredVcf = File.createTempFile("unfiltered", ".vcf");
        final List<String> args = Arrays.asList(
                "-I", samFile.getAbsolutePath(),
                "-" + M2ArgumentCollection.TUMOR_SAMPLE_SHORT_NAME, M2TestingUtils.DEFAULT_SAMPLE_NAME,
                "-R", hg19_chr1_1M_Reference,
                "-O", unfilteredVcf.getAbsolutePath());
        runCommandLine(args);

        final File filteredVcf = File.createTempFile("filtered", ".vcf");
        final String[] filteringArgs = makeCommandLineArgs(Arrays.asList(
                "-V", unfilteredVcf.getAbsolutePath(),
                "-O", filteredVcf.getAbsolutePath()),
                FilterMutectCalls.class.getSimpleName());
        new Main().instanceMain(filteringArgs);

        final Optional<VariantContext> vc = VariantContextTestUtils.streamVcf(filteredVcf).findAny();
        Assert.assertTrue(vc.isPresent());
        Assert.assertEquals(vc.get().getStart(), M2TestingUtils.DEFAULT_SNP_POSITION);
        Assert.assertFalse(vc.get().getFilters().contains(GATKVCFConstants.MEDIAN_BASE_QUALITY_FILTER_NAME));
    }

    private void doMutect2Test(
            final String inputBam,
            final String tumorSample,
            final String interval,
            final boolean createBamout,
            final boolean createBamoutIndex,
            final boolean createBamoutMD5) {
        final File tempDir = GATKBaseTest.createTempDir("mutect2");
        final File outputVcf = new File(tempDir,"output.vcf");
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

    // tumor bam, tumor sample name, normal bam, normal sample name, truth vcf, required sensitivity, tumor only
    @DataProvider(name = "dreamSyntheticData")
    public Object[][] dreamSyntheticData() {
        return new Object[][]{
                {new File(DREAM_BAMS_DIR, "tumor_1.bam"), "tumor sample", new File(DREAM_BAMS_DIR, "normal_1.bam"), "synthetic.challenge.set1.normal", new File(DREAM_VCFS_DIR, "sample_1.vcf"), new File(DREAM_MASKS_DIR, "mask1.list"), 0.97, false},
                {new File(DREAM_BAMS_DIR, "tumor_2.bam"), "background.synth.challenge2.snvs.svs.tumorbackground", new File(DREAM_BAMS_DIR, "normal_2.bam"), "synthetic.challenge.set2.normal", new File(DREAM_VCFS_DIR, "sample_2.vcf"), new File(DREAM_MASKS_DIR, "mask2.list"), 0.95, false},
                {new File(DREAM_BAMS_DIR, "tumor_2.bam"), "background.synth.challenge2.snvs.svs.tumorbackground", new File(DREAM_BAMS_DIR, "normal_2.bam"), "synthetic.challenge.set2.normal", new File(DREAM_VCFS_DIR, "sample_2.vcf"), new File(DREAM_MASKS_DIR, "mask2.list"), 0.95, true},
                {new File(DREAM_BAMS_DIR, "tumor_3.bam"), "IS3.snv.indel.sv", new File(DREAM_BAMS_DIR, "normal_3.bam"), "G15512.prenormal.sorted", new File(DREAM_VCFS_DIR, "sample_3.vcf"), new File(DREAM_MASKS_DIR, "mask3.list"), 0.90, false},
                {new File(DREAM_BAMS_DIR, "tumor_4.bam"), "synthetic.challenge.set4.tumour", new File(DREAM_BAMS_DIR, "normal_4.bam"), "synthetic.challenge.set4.normal", new File(DREAM_VCFS_DIR, "sample_4.vcf"), new File(DREAM_MASKS_DIR, "mask4.list"), 0.65, false},
                {new File(DREAM_BAMS_DIR, "tumor_4.bam"), "synthetic.challenge.set4.tumour", new File(DREAM_BAMS_DIR, "normal_4.bam"), "synthetic.challenge.set4.normal", new File(DREAM_VCFS_DIR, "sample_4.vcf"), new File(DREAM_MASKS_DIR, "mask4.list"), 0.65, true}

        };
    }

    // tumor bam, tumor sample name, normal bam, normal sample name, truth vcf, required sensitivity
    @DataProvider(name = "dreamSyntheticDataSample1")
    public Object[][] dreamSyntheticDataSample1() {
        return new Object[][]{
                {new File(DREAM_BAMS_DIR, "tumor_1.bam"), "tumor sample", new File(DREAM_BAMS_DIR, "normal_1.bam"), "synthetic.challenge.set1.normal"}
        };
    }

    //TODO: bring this to HaplotypeCallerIntegrationTest
    private Pair<Double, Double> calculateConcordance(final File outputVcf, final File truthVcf ) {
        final Set<String> outputKeys = VariantContextTestUtils.streamVcf(outputVcf)
                .filter(vc -> vc.getFilters().isEmpty())
                .filter(vc -> ! vc.isSymbolicOrSV())
                .map(vc -> keyForVariant(vc)).collect(Collectors.toSet());
        final Set<String> truthKeys = VariantContextTestUtils.streamVcf(truthVcf)
                .filter(vc -> vc.getFilters().isEmpty())
                .filter(vc -> ! vc.isSymbolicOrSV())
                .map(vc -> keyForVariant(vc)).collect(Collectors.toSet());

        final long truePositives = outputKeys.stream().filter(truthKeys::contains).count();
        final long falsePositives = outputKeys.size() - truePositives;

        final double sensitivity = (double) truePositives / truthKeys.size();
        final double fdr = (double) falsePositives / outputKeys.size();
        return new ImmutablePair<>(sensitivity, fdr);
    }

    private static String keyForVariant( final VariantContext variant ) {
        return String.format("%s:%d-%d %s", variant.getContig(), variant.getStart(), variant.getEnd(), variant.getAlleles());
    }
}