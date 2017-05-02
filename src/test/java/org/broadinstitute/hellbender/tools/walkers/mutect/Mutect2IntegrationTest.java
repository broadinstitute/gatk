package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.variant.variantcontext.VariantContext;
import junit.framework.Assert;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.Main;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.tools.walkers.validation.ConcordanceSummaryRecord;
import org.broadinstitute.hellbender.utils.Utils;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

/**
 * Created by davidben on 9/1/16.
 */
public class Mutect2IntegrationTest extends CommandLineProgramTest {
    // positions 10,000,000 - 11,000,000 of chr 20 and with most annotations removed
    private static final File GNOMAD = new File(publicTestDir, "very-small-gnomad.vcf");
    private static final String DREAM_BAMS_DIR = largeFileTestDir + "mutect/dream_synthetic_bams/";
    private static final String DREAM_VCFS_DIR = publicTestDir + "org/broadinstitute/hellbender/tools/mutect/dream/vcfs/";
    private static final String DREAM_MASKS_DIR = publicTestDir + "org/broadinstitute/hellbender/tools/mutect/dream/masks/";

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
     * Sample2: 80% pure monoclonal sample, SNVs only
     * Sample3: pure triclonal sample, subclone minor allele frequencies are 1/2, 1/3, and 1/5, SNVs and indels
     * Sample 4: 80% biclonal sample, subclone minor allele fractions are 50% and 35%, SNVs and indels
     *
     * @throws Exception
     */
    @Test(dataProvider = "dreamSyntheticData")
    public void testDreamTumorNormal(final File tumorBam, final String tumorSample, final File normalBam, final String normalSample,
                                     final File truthVcf, final File mask, final double requiredSensitivity) throws Exception {
        Utils.resetRandomGenerator();
        final File unfilteredVcf = createTempFile("unfiltered", ".vcf");
        final File filteredVcf = createTempFile("filtered", ".vcf");

        final String[] args = {
                "-I", tumorBam.getAbsolutePath(),
                "-tumor", tumorSample,
                "-I", normalBam.getAbsolutePath(),
                "-normal", normalSample,
                "-R", b37_reference_20_21,
                "-L", "20",
                "-germline_resource", GNOMAD.getAbsolutePath(),
                "-XL", mask.getAbsolutePath(),
                "-O", unfilteredVcf.getAbsolutePath()
        };

        runCommandLine(args);

        // run FilterMutectCalls
        new Main().instanceMain(makeCommandLineArgs(Arrays.asList("-V", unfilteredVcf.getAbsolutePath(), "-O", filteredVcf.getAbsolutePath()), "FilterMutectCalls"));


        // run Concordance
        final File concordanceSummary = createTempFile("concordance", ".txt");
        new Main().instanceMain(makeCommandLineArgs(Arrays.asList("-truth", truthVcf.getAbsolutePath(), "-eval", unfilteredVcf.getAbsolutePath(), "-L", "20", "-XL", mask.getAbsolutePath(), "-summary", concordanceSummary.getAbsolutePath()), "Concordance"));

        final List<ConcordanceSummaryRecord> summaryRecords = new ConcordanceSummaryRecord.Reader(concordanceSummary).toList();
        summaryRecords.forEach(rec -> {
            if (rec.getTruePositives() + rec.getFalseNegatives() > 0) {
                Assert.assertTrue(rec.getSensitivity() > requiredSensitivity);
                Assert.assertTrue(rec.getPrecision() > 0.5);
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
                "-tumor", tumorSample,
                "-I", normalBam.getAbsolutePath(),
                "-normal", normalSample,
                "-R", b37_reference_20_21,
                "-L", "20",
                "-O", ponVcf.getAbsolutePath()
        };

        runCommandLine(createPonArgs);

        final File unfilteredVcf = createTempFile("unfiltered", ".vcf");
        final File filteredVcf = createTempFile("filtered", ".vcf");
        final String[] callWithPonArgs = {
                "-I", tumorBam.getAbsolutePath(),
                "-tumor", tumorSample,
                "-I", normalBam.getAbsolutePath(),
                "-normal", normalSample,
                "-normal_panel", ponVcf.getAbsolutePath(),
                "-R", b37_reference_20_21,
                "-L", "20",
                "-O", unfilteredVcf.getAbsolutePath()
        };

        runCommandLine(callWithPonArgs);

        // run FilterMutectCalls
        new Main().instanceMain(makeCommandLineArgs(Arrays.asList("-V", unfilteredVcf.getAbsolutePath(), "-O", filteredVcf.getAbsolutePath()), "FilterMutectCalls"));


        final long numVariants = StreamSupport.stream(new FeatureDataSource<VariantContext>(filteredVcf).spliterator(), false)
                .filter(vc -> vc.getFilters().isEmpty()).count();

        Assert.assertEquals(numVariants, 0);
    }

    // run tumor-only using the original DREAM synthetic sample 1 tumor and normal restricted to
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
                "-tumor", tumorName,
                "-I", normalBam.getAbsolutePath(),
                "-normal", normalName,
                "-R", b37_reference_20_21,
                "-L", "20:10000000-10100000", // this is 1/3 of the chr 20 interval of our mini-dbSNP
                "-O", outputVcf.getAbsolutePath()
        };

        runCommandLine(args);
        final long numVariants = StreamSupport.stream(new FeatureDataSource<VariantContext>(outputVcf).spliterator(), false).count();
        Assert.assertTrue(numVariants < 4);
    }

    // run tumor-only using our mini gnomAD on NA12878, which is not a tumor
    // we're just making sure nothing blows up
    @Test
    public void testTumorOnly() throws Exception {
        Utils.resetRandomGenerator();
        final File unfilteredVcf = createTempFile("unfiltered", ".vcf");
        final File filteredVcf = createTempFile("filtered", ".vcf");

        final String[] args = {
                "-I", NA12878_20_21_WGS_bam,
                "-tumor", "NA12878",
                "-R", b37_reference_20_21,
                "-L", "20:10000000-10010000",
                "-germline_resource", GNOMAD.getAbsolutePath(),
                "-O", unfilteredVcf.getAbsolutePath()
        };

        runCommandLine(args);

        // run FilterMutectCalls
        new Main().instanceMain(makeCommandLineArgs(Arrays.asList("-V", unfilteredVcf.getAbsolutePath(), "-O", filteredVcf.getAbsolutePath()), "FilterMutectCalls"));


        final long numVariantsBeforeFiltering = StreamSupport.stream(new FeatureDataSource<VariantContext>(filteredVcf).spliterator(), false).count();

        final long numVariantsPassingFilters = StreamSupport.stream(new FeatureDataSource<VariantContext>(filteredVcf).spliterator(), false)
                .filter(vc -> vc.getFilters().isEmpty()).count();

        // just a sanity check that this bam has some germline variants on this interval so that our test doesn't pass trivially!
        Assert.assertTrue(numVariantsBeforeFiltering > 50);

        // every variant on this interval in this sample is in gnomAD
        Assert.assertTrue(numVariantsPassingFilters < 2);
    }

    // test that ReadFilterLibrary.NON_ZERO_REFERENCE_LENGTH_ALIGNMENT removes reads that consume zero reference bases
    // e.g. read name HAVCYADXX150109:1:2102:20528:2129 with cigar 23S53I
    @Test
    public void testReadsThatConsumeZeroReferenceReads() throws Exception {
        final String CONSUMES_ZERO_REFERENCE_BASES = publicTestDir + "org/broadinstitute/hellbender/tools/mutect/na12878-chr20-consumes-zero-reference-bases.bam";
        final File outputVcf = createTempFile("output", ".vcf");
        final String[] args = {
                "-I", CONSUMES_ZERO_REFERENCE_BASES,
                "-tumor", "SM-612V3",
                "-R", b37_reference_20_21,
                "-O", outputVcf.getAbsolutePath()
        };
        runCommandLine(args);
    }

    // tumor bam, tumor sample name, normal bam, normal sample name, truth vcf, required sensitivity
    @DataProvider(name = "dreamSyntheticData")
    public Object[][] dreamSyntheticData() {
        return new Object[][]{
                {new File(DREAM_BAMS_DIR, "tumor_1.bam"), "synthetic.challenge.set1.tumor", new File(DREAM_BAMS_DIR, "normal_1.bam"), "synthetic.challenge.set1.normal", new File(DREAM_VCFS_DIR, "sample_1.vcf"), new File(DREAM_MASKS_DIR, "mask1.list"), 0.97},
                {new File(DREAM_BAMS_DIR, "tumor_2.bam"), "background.synth.challenge2.snvs.svs.tumorbackground", new File(DREAM_BAMS_DIR, "normal_2.bam"), "synthetic.challenge.set2.normal", new File(DREAM_VCFS_DIR, "sample_2.vcf"), new File(DREAM_MASKS_DIR, "mask2.list"), 0.95},
                {new File(DREAM_BAMS_DIR, "tumor_3.bam"), "IS3.snv.indel.sv", new File(DREAM_BAMS_DIR, "normal_3.bam"), "G15512.prenormal.sorted", new File(DREAM_VCFS_DIR, "sample_3.vcf"), new File(DREAM_MASKS_DIR, "mask3.list"), 0.90},
                {new File(DREAM_BAMS_DIR, "tumor_4.bam"), "synthetic.challenge.set4.tumour", new File(DREAM_BAMS_DIR, "normal_4.bam"), "synthetic.challenge.set4.normal", new File(DREAM_VCFS_DIR, "sample_4.vcf"), new File(DREAM_MASKS_DIR, "mask4.list"), 0.65}
        };
    }

    // tumor bam, tumor sample name, normal bam, normal sample name, truth vcf, required sensitivity
    @DataProvider(name = "dreamSyntheticDataSample1")
    public Object[][] dreamSyntheticDataSample1() {
        return new Object[][]{
                {new File(DREAM_BAMS_DIR, "tumor_1.bam"), "synthetic.challenge.set1.tumor", new File(DREAM_BAMS_DIR, "normal_1.bam"), "synthetic.challenge.set1.normal"}
        };
    }

    //TODO: bring this to HaplotypeCallerIntegrationTest
    private Pair<Double, Double> calculateConcordance(final File outputVcf, final File truthVcf ) {
        final Set<String> outputKeys = StreamSupport.stream(new FeatureDataSource<VariantContext>(outputVcf).spliterator(), false)
                .filter(vc -> vc.getFilters().isEmpty())
                .filter(vc -> ! vc.isSymbolicOrSV())
                .map(vc -> keyForVariant(vc)).collect(Collectors.toSet());
        final Set<String> truthKeys = StreamSupport.stream(new FeatureDataSource<VariantContext>(truthVcf).spliterator(), false)
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