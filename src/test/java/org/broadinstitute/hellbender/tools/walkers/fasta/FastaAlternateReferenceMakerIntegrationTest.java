package org.broadinstitute.hellbender.tools.walkers.fasta;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.FastaTestUtils;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;

public class FastaAlternateReferenceMakerIntegrationTest extends CommandLineProgramTest {
    private final File spanningDelTestFile = getTestFile("spanningDel.delOnly.starFirst.vcf");
    private final File iupacTestFile = getTestFile("NA12878.WGS.b37.chr20.firstMB.vcf.gz");
    private final File snp1 = getTestFile("snps.vcf");
    private final File snp2 = getTestFile("snps2.vcf");
    private final File snpsAndIndels = getTestFile("snpsAndIndels.vcf");
    private final File overlapsIndels = getTestFile("overlapsIndelsMask.vcf");

    @Test
    public void testAlternateReferenceContiguousSameContig() {
        // Show that FastaAlternateReferenceMaker behaves the same as FastaReferenceMaker across contiguous intervals on the same contig.
        // Note that there are variant locations in this interval.
       ArgumentsBuilder args = new ArgumentsBuilder();
        final File out = createTempFile("test", "fasta");
        args.addVCF(getTestFile("NA12878.chr1_10mb_11mb.slx.indels.vcf"))
                .addReference(new File(b37Reference))
                .add("L", "1:10,000,100-10,000,200")
                .add("L", "1:10,000,201-10,000,301")
                .addOutput(out);
        runCommandLine(args);
        FastaTestUtils.assertFastaFilesContainTheSameSequence(out.toPath(), getTestFile("expected_alternate_reference_contiguous_same_contig.fasta").toPath());
    }

    @Test
    public void testAlternateReferenceContiguousDiffContigs() {
        // Show that FastaAlternateReferenceMaker behaves the same as FastaReferenceMaker across contiguous intervals on different contigs.
        // Note that there are variant locations in this interval.
        ArgumentsBuilder args = new ArgumentsBuilder();
        final File out = createTempFile("test", "fasta");
        args.addVCF(getTestFile("NA12878.chr1_10mb_11mb.slx.indels.vcf"))
                .addReference(new File(b37Reference))
                .add("L", "1:10,000,100-10,000,200")
                .add("L", "2:10,000,201-10,000,301")
                .addOutput(out);
        runCommandLine(args);
        FastaTestUtils.assertFastaFilesContainTheSameSequence(out.toPath(), getTestFile("expected_alternate_reference_contiguous_diff_contig.fasta").toPath());
    }

    @DataProvider
    public Object[][] getSnpMaskVariants(){
        return new Object[][]{
                {snp1, snp2, getTestFile("expected_snp1_mask2.fasta"), false},
                {snp1, snp2, getTestFile("expected_snp1_mask2_priority.fasta"), true},
                {snp2, snp1, getTestFile("expected_snp2_mask1.fasta"), false},
                {snp2, snp1, getTestFile("expected_snp2_mask1_priority.fasta"), true},
                {snpsAndIndels, snp1, getTestFile("expected_snpsAndIndels_mask1.fasta"), false},
                {snpsAndIndels, snp1, getTestFile("expected_snpsAndIndels_mask1_priority.fasta"), true},
                {snpsAndIndels, overlapsIndels, getTestFile("expected_snpsAndIndels_maskOverlaps.fasta"), false},
                {snpsAndIndels, overlapsIndels, getTestFile("expected_snpsAndIndels_maskOverlaps_priority.fasta"), true},
        };
    }

    @Test(dataProvider = "getSnpMaskVariants")
    public void testSnpMask(File vcf, File mask, File expected, boolean priority) {
        ArgumentsBuilder args = new ArgumentsBuilder();
        final File out = createTempFile("test", "fasta");
        args.addVCF(vcf)
                .addReference(new File(b37Reference))
                .add(FastaAlternateReferenceMaker.SNP_MASK_LONG_NAME, mask)
                .add("L", "1:10,000,000-10,000,100")
                .add(FastaAlternateReferenceMaker.SNP_MASK_PRIORITY_LONG_NAME, priority)
                .addOutput(out);
        runCommandLine(args);
        FastaTestUtils.assertFastaFilesContainTheSameSequence(out.toPath(), expected.toPath());
    }

    @DataProvider
    public Object[][] noMask() {
        return new Object[][]{
                {snp1, getTestFile("expected_snp1.fasta")},
                {snp2, getTestFile("expected_snp2.fasta")},
        };
    }

    @Test(dataProvider = "noMask")
    public void testNoMask(File vcf, File expected) {
            ArgumentsBuilder args = new ArgumentsBuilder();
            final File out = createTempFile("test", "fasta");
            args.addVCF(vcf)
                    .addReference(new File(b37Reference))
                    .add("L", "1:10,000,000-10,000,100")
                    .addOutput(out);
            runCommandLine(args);
            FastaTestUtils.assertFastaFilesContainTheSameSequence(out.toPath(), expected.toPath());
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testBadIupacInput() {
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.addVCF(iupacTestFile)
                .addOutput( createTempFile("alternate", "fasta"))
                .addReference(new File(b37Reference))
                .add(FastaAlternateReferenceMaker.USE_IUPAC_SAMPLE_LONG_NAME, "SAMPLE_DOESNT_EXIST");
        runCommandLine(args);
    }

    @Test
    public void testIupac() {
        final File out =  createTempFile("alternate", "fasta");
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.addVCF(iupacTestFile)
                .addOutput(out)
                .addReference(new File(b37Reference))
                .add(FastaAlternateReferenceMaker.USE_IUPAC_SAMPLE_LONG_NAME, "NA12878")
                .add("L", "20:61050-66380");

        runCommandLine(args);
        FastaTestUtils.assertFastaFilesContainTheSameSequence(out.toPath(), getTestFile("expected_test_iupac.fasta").toPath());
    }

    @Test
    void testSpanDel() {
        final File out = createTempFile("alternate", "fasta");
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.addVCF(spanningDelTestFile)
                .addOutput(out)
                .addReference(new File(b37Reference))
                .add("L", "1:1273247");
        runCommandLine(args);
        FastaTestUtils.assertFastaFilesContainTheSameSequence(out.toPath(), getTestFile("iupac_sample_NA3.fasta").toPath());
    }

    @DataProvider(name = "iupacSample")
    public Object[][] getIupacSampleData() {
        return new Object[][]{
                {"NA1", getTestFile("iupac_sample_NA1.fasta")},
               // {"NA2", getTestFile("iupac_sample_NA2.fasta")}, this test disabled because gatk4 won't write a fasta with an empty reference sequence
                {"NA3", getTestFile("iupac_sample_NA3.fasta")}
        };
    }

    @Test(dataProvider = "iupacSample")
    void testSpanDelIUPAC(final String sample, final File expected) {
        final ArgumentsBuilder args = new ArgumentsBuilder();
        final File out = createTempFile("alternate", "fasta");
        args.addVCF(spanningDelTestFile)
                .addOutput(out)
                .addReference(new File(b37Reference))
                .add("L", "1:1273247")
                .add("use-iupac-sample", sample);
        runCommandLine(args);
        FastaTestUtils.assertFastaFilesContainTheSameSequence(out.toPath(), expected.toPath());
    }
}