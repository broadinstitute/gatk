package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.engine.ReferenceFileSource;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCaller;
import org.broadinstitute.hellbender.tools.walkers.mutect.InferOriginalRead;
import org.broadinstitute.hellbender.tools.walkers.mutect.M2TestingUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadCoordinateComparator;
import org.broadinstitute.hellbender.utils.read.SAMFileGATKReadWriter;
import org.broadinstitute.hellbender.utils.runtime.ProcessController;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.IntStream;

public class InferOriginalReadTest extends CommandLineProgramTest {
    private final String hg19 = "/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta";
    private final String tp53TestDir = "/dsde/working/tsato/consensus/tp53/test/";
    private final String tp53_AGA_TGA = tp53TestDir + "Jonna_Grimsby_A04_denovo_bloodbiopsy_1pct_rep1.tp53.AGG-TGA.grouped.bam";
    private final String tp53_full_grouped = tp53TestDir + "Jonna_Grimsby_A04_denovo_bloodbiopsy_1pct_rep1.tp53.umi_grouped.bam";

    private final String mutectTestDir = "/Users/tsato/workspace/gatk/src/test/java/org/broadinstitute/hellbender/tools/walkers/mutect/";
    private final String fgbioGroupByUMIScript = mutectTestDir + "fgbio_group_reads_by_umi.sh";
    private final String fgbioConsensusScript = mutectTestDir + "fgbio_consensus.sh";
    // /Users/tsato/workspace/gatk/src/main/java/org/broadinstitute/hellbender/tools/walkers/mutect/sort.sh
    private final String sortScript = mutectTestDir + "sort.sh";

    private final String homeDir = "/dsde/working/tsato/consensus/bams/synthetic-test/";
    private final String testDir = "/dsde/working/tsato/consensus/tp53/test/";

    @Test
    public void test(){
        final File out = new File("/dsde/working/tsato/consensus/tp53/test/out.bam");
        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addArgument("R", hg19)
                .addArgument("I", tp53_AGA_TGA)
                .addArgument("O", out.getAbsolutePath());
        runCommandLine(args, InferOriginalRead.class.getSimpleName());
        int d = 3;
    }

    // 1/22/20 --- Next goal is to thread this pipeline through
    @Test
    public void testTp53(){
        final File out = new File(testDir + "tp53_full_consensus.bam");
        final String testBam = "/dsde/working/tsato/consensus/tp53/Jonna_Grimsby_A04_denovo_bloodbiopsy_1pct_rep1.tp53.CTG-TTC.grouped.bam";
        final boolean test = false;
        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addArgument("R", hg19)
                .addArgument("I", test ? testBam : tp53_full_grouped)
                .addArgument("O", out.getAbsolutePath());
        runCommandLine(args, InferOriginalRead.class.getSimpleName());
        // runProcess(new ProcessController(), new String[]{ sortScript, out.getAbsolutePath(), testDir, "Jonna_Grimsby_A04_denovo_bloodbiopsy_1pct_rep1.tp53.CTG-TTC." });
    }

    @Test
    public void testCleanDeletion(){
        final File out = new File(homeDir + "clean_deletion.consensus.bam");
        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addArgument("R", hg19) // need to take care of the reference and whatnot....
                .addArgument("I", homeDir + "clean_deletion.grouped.bam")
                .addArgument("O", out.getAbsolutePath());
        runCommandLine(args, InferOriginalRead.class.getSimpleName());
        int d = 3;
    }

    @Test
    public void testHalotypeCaller(){
        final File out = new File(testDir + "hc.vcf");
        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addArgument("R", hg19) // need to take care of the reference and whatnot....
                .addArgument("I", testDir + "Jonna_Grimsby_A04_denovo_bloodbiopsy_1pct_rep1.tp53.AGG-TGA.bam")
                .addArgument("O", out.getAbsolutePath());
        runCommandLine(args, HaplotypeCaller.class.getSimpleName());
        int d = 3;
    }

    @DataProvider(name = "whatever")
    public Object[][] whateverData() {
        return new Object[][]{
                {Arrays.asList(2), Arrays.asList(8), Arrays.asList(0), Arrays.asList(0), "ref2_insertion8" },
        };
    }

    @Test(dataProvider = "whatever")
    public void testBetter(final List<Integer> refCounts,
                           final List<Integer> insertionCounts,
                           final List<Integer> twoBPDeletionCounts,
                           final List<Integer> threeBPDeletionCounts,
                           final String testName) throws IOException {
        final File originalSam = createSyntheticSam(2, 8, 0, 0, "sample",
                homeDir + "ref2_insertion8.bam", true);

        // Run fgbio GroupByUMI on this subset bam
        final String bam = homeDir + testName + ".bam";
        final String groupedBam = homeDir + testName + ".grouped.bam";
        runProcess(new ProcessController(), new String[]{ fgbioGroupByUMIScript, bam, groupedBam, homeDir });

        // For comparison, maybe unnecessary
        runProcess(new ProcessController(), new String[]{ fgbioConsensusScript, groupedBam, testName, homeDir });

        final File out = new File(homeDir + testName + ".my.consensus.bam");
        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addArgument("R", hg19)
                .addArgument("I", groupedBam)
                .addArgument("O", out.getAbsolutePath());
        runCommandLine(args, InferOriginalRead.class.getSimpleName());

        // Don't want to have to do this but whatever--can't the
        runProcess(new ProcessController(), new String[]{ sortScript, groupedBam, testName });
        int d = 3;
    }

    @Test
    public void makeTestSamFiles() throws IOException {
        // 2 ref, 8 insertion reads
        final File ref2Insertion8 = createSyntheticSam(2, 8, 0, 0, "sample",
                homeDir + "ref2_insertion8.bam", true);

        // 5 ref, 5 deletion reads
        final File ref5Deletion5 = createSyntheticSam(5, 0, 5, 0, "sample",
                homeDir + "ref5_deletion5.bam", true);

        // just 5 deletion reads
        final File cleanDeletion = createSyntheticSam(0, 0, 0, 5, "sample",
                homeDir + "clean_deletion.bam", true);

        final File pcrError = createSyntheticSam(0, 0, 3, 0, "sample",
                homeDir + "pcr_error.bam", false);

        // Three 2-bp deletions and Seven 3-bp deletions
        final File twoTypesDeletions = createSyntheticSam(0, 0, 3, 7, "sample",
                homeDir + "two_types_deletion.bam", true);

        // Run fgbio GroupByUMI on this subset bam
        final String bam = homeDir + "two_types_deletion.bam";
        final String groupedBam = homeDir + "two_types_deletion.grouped.bam";
        runProcess(new ProcessController(), new String[]{ fgbioGroupByUMIScript, bam, groupedBam, homeDir });
        runProcess(new ProcessController(), new String[]{ fgbioConsensusScript, groupedBam, "two_types_deletion", homeDir });

        final String bam2 = homeDir + "ref5_deletion5.bam";
        final String groupedBam2 = homeDir + "ref5_deletion5.grouped.bam";
        runProcess(new ProcessController(), new String[]{ fgbioGroupByUMIScript, bam2, groupedBam2, homeDir });
        runProcess(new ProcessController(), new String[]{ fgbioConsensusScript, groupedBam2, "ref5_deletion5", homeDir });

        final String bam3 = homeDir + "clean_deletion.bam";
        final String groupedBam3 = homeDir + "clean_deletion.grouped.bam";
        runProcess(new ProcessController(), new String[]{ fgbioGroupByUMIScript, bam3, groupedBam3, homeDir });
        runProcess(new ProcessController(), new String[]{ fgbioConsensusScript, groupedBam3, "ref5_deletion5", homeDir });

        int d = 3;
    }

    private File createSyntheticSam(final int refCount, final int insertionCount, final int twoBPDeletionCount,
                                    final int threeBPDeletionCount, final String sampleName, final String bamName,
                                    final boolean includeF2R1) throws IOException {
        final ReferenceDataSource ref = new ReferenceFileSource(Paths.get(hg19));

        final File samFile = new File(bamName);
        final SAMFileHeader samHeader = M2TestingUtils.createSamHeader(sampleName);
        samHeader.setSequenceDictionary(ref.getSequenceDictionary());
        final SAMFileGATKReadWriter writer = M2TestingUtils.getBareBonesSamWriter(samFile, samHeader); // ts: is reference=null ok?

        final List<GATKRead> reads = new ArrayList<>(20);
        reads.addAll(M2TestingUtils.createReads(refCount, M2TestingUtils.DEFAULT_REF_BASES,
                samHeader, (byte)30, "ref", M2TestingUtils.CIGAR_REF, M2TestingUtils.DEFAULT_UMI, true, true));

        if (insertionCount > 0){
            // i = 0 creates F1R2 reads, i = 1 creates F1R2 reads. Not reader friendly.
            IntStream.range(0, 2).forEach(i -> {
                final List<GATKRead> insertionReads = M2TestingUtils.createReads(insertionCount, M2TestingUtils.TTT_INSERTION,
                        samHeader, (byte)30, "del", M2TestingUtils.CIGAR_TTT_INSERTION, M2TestingUtils.DEFAULT_UMI, true, i % 2 == 0);
                reads.addAll(insertionReads);
            });

        }

        if (twoBPDeletionCount > 0){
            IntStream.range(0, 2).forEach(i -> {
                final List<GATKRead> twoBPDeletionReads = M2TestingUtils.createReads(twoBPDeletionCount, M2TestingUtils.AC_DELETION,
                        samHeader, (byte) 30, "del_two", M2TestingUtils.CIGAR_AC_DELETION, M2TestingUtils.DEFAULT_UMI, true, i % 2 == 0);
                reads.addAll(twoBPDeletionReads);
            });
        }

        if (threeBPDeletionCount > 0){
            IntStream.range(0, 2).forEach(i -> {
                final List<GATKRead> threeBPDeletionReads = M2TestingUtils.createReads(threeBPDeletionCount, M2TestingUtils.ACA_DELETION,
                        samHeader, (byte) 30, "del_three", M2TestingUtils.CIGAR_ACA_DELETION, M2TestingUtils.DEFAULT_UMI, true, i % 2 == 0);
                reads.addAll(threeBPDeletionReads);
            });
        }

        reads.sort(new ReadCoordinateComparator(samHeader));
        reads.forEach(writer::addRead);

        writer.close(); // closing the writer writes to the file
        return samFile;
    }
}