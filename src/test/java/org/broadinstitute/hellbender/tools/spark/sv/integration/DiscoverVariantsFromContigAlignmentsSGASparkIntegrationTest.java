package org.broadinstitute.hellbender.tools.spark.sv.integration;


import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.LocatedFileStatus;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.fs.RemoteIterator;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.test.MiniClusterUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class DiscoverVariantsFromContigAlignmentsSGASparkIntegrationTest extends CommandLineProgramTest {

    private final String contigs =
            "1\t>contig-0 250 0|CATTCTATAACAGCTATGAAACCATTTGTTATAGCTAATTCGTAGAAGCCAGGCTAAATGAGCCCTAATCAACTAGAACAAGCATTCAGAACTCAAAATCCATTTACATGTTAATTAAGTAAATAGTATAAGCTACACGTTGAGCCTACAACTCTTGCTGGAGTCTCCAGCTGTAAGAGCCTTAACTCCTTCCTTTCATAAAGGCCAAGAAAAGAGTGAGTTAGTCCCTATTCGGGGAGTGTTGGTTGGT|>contig-1 151 0|AGGCACAATCTCCCCAGTCCAGAATCTAAATGGCATTATATTAACATTTTCCTATACCACTAATTAAACATAATGCTACATTCATTTTATTTTCTTCTGTCTCATTCTATAACAGCTATGAAACCATTTGTTATAGCTAATTCGTAGAAGC|>contig-2 152 0|TCACATTTTGGCATTAAAAAAAAGTATCATCTTAGAGTACCAAATATGTGTCTTCTTATAGATGTGTGCCACGGGCCTCGCAGTGAACATCACAGCTCTGCAGATGCACAAACCTTGTCTTTATGAGTCTATATATTAAATAAAAGTGGTAG|>contig-3 195 0|AAGTCAAGTAATGGGATAAAAGAATGATGATCGAAGGAGATGATTCATGATAACTATGTGCCCTTATGTGAAACCTAGCAACAATTTCTTACTTATGATAAAGGCACAATCTCCCCAGTCCAGAATCTAAATGGCATTATATTAACATTTTCCTATACCACTAATTAAACATAATGCTACATTCATTTTATTTTC|>contig-4 151 0|CCTTATGTGAAACCTAGCAACAATTTCTTACTTATGATAAAGGCACAATCTCCCCAGTCCAGAATCTAAATGGCATTATATTAACATTTTCCTATACCACTAATTAAACATAATGCTACATTCATTTTATTTTCTTCTGCCTCATTCTATA|>contig-5 151 0|GAGTCTATATATTAAATAAAAGTGGTAGAATTTTGGATGGGTGCAAGTCAAGTAATGGGATAAAAGAATGATGATCGAAGGAGATGATTCATGATAACTATGTGCGCTTATGTGAAACCTAGCAACAATTTCTTACTTATGATAAAGGCAC|>contig-6 241 0|CTTATGTGAAACCTAGCAACAATTTCTTACTTATGATAAAGGCACAATCTCCCCAGTCCAGAATCTAAATGGCATTATATTAACATTTTCCTATACCACTAATTAAACATAATGCTACATTCATTTTATTTTCTTCTGCCTCATTCTATAACAGCTATGAAACCATTTGTTATAGCTAATTCGTAGAAGCCAGGCTAAATGAGCCCTAATCAACTAGAACAAGCATTCAGAACTCAAAATC|>contig-7 151 0|TATTAAATAAAAGTGGTAGAATTTTGGATGGGTGCAAGTCAAGTAATGGGATAAAAGAATGATGATCGAAGGAGATGACTCATGATAACTAAGTGCCCTTATGTGAAACCTAGCAACAATTTCTTACTTATGATAAAGGCACAATCTCCCC|>contig-8 194 0|CCACTAATTAAACATAATGCTACATTCATTTTATTTTCTTCTGCCTCATTCTATAACAGCTATGAAACCATTTGTTATAGCTAATTCGTAGAAGCCAGGCTAAATGAGCCCTAATCAACTAGAACAAGCATTCAGAACTCAAAATCCATTTACATGTTAATTAAGTAAATAGTATAAGCTACACGTTGAGCCTA|>contig-9 2009 0|CATGAAGTACATCTATAACAACGTATTATTAATTCTGCTATTTAATTGATATTTATATATTAAAAATCTGTATTATTTGTATATATTTACATACACATTAAGTACAACTATAACAATGTATTATTTAATTATGCTACTTAATTGATATTTATGTATTAAAAATCTGTATTATTTGTATATTTTGTATATTTACCTGACCCTTGGGCAGGAATGTTCCAAACAATAAAACGCAGTCAATCAACTCCAACTGATTAATACACACCACGCACTGTTATACATTAAAACACATCTTACCTGTCAACATGTCTTGCCACCCTCCCCAAGATCAGTTTCACTCATCCCTGATTGCCTCAAAGTCAATCAGTGACCCATCCAATTGCATGTATAAAGTCAATCCTGAGTAGAATTATAACTCTACCTCCCATGTATATCCCATTTTTATTGGAAGTGAAACAACAGCCCAATTCTCTCTCTGGGTTTGTGCCCTGCCCTTGCTAGACTGACGCCTTAGCACACACAGATGTCTTCTTAGTTTTACTGACAGAGGAAAAAACACCTTTCTAAAAGATGAAAGCCACGCATCACCCATTACTTCCAAGTTAGGAATATTCAACAGATTCTTTTTTTCTTTTTATGAGACAGGGTCTTGCACTGCCGCCCAGTGGCACAATCTTGGCTCACTGCAACCTCTACCTCCTGAGTTCAAGTGATTCTTGCGCCTCAGCCTCCCAAGTAGCTGGGATTACCCTTGTACCAGCACGCCTGGCTAATTTTTTGGTATTTTTAGTAGAGATGGGGTTTCACCAGGTAGGCCAGGCTGGTCTTGAACTCCTGGCCTCAAGTGATACACCTGCCTCAGCCTCTCTAAGTGCTGGGATTACAGGCGTGAGCCACTGCCACTGGTCCTTCAACAGATTGATTCTAATTAGCCAATCAAAGACAAGGATCCACAATGTCAATACAAATGGGGACTAACTATTGTTTTACTTCCCTATACACGTTGCTCTTTCAGTAGATGCAATAAAGTACTGGTAAAACCAGAGGTGGCTACCATCACGATGATGTCAACAGGAGGGACAGTCAGCACTAAGCCCAGAAGGTGTCAAACACTCCGCAGGAGAAATGCGCCATGCAACGGACATGAAGATGATCTGACACTCTTCACGTGGTTTTCAGATGGAAACGTGGCTACGAAAGCATCAACCTCATTATCATCCATCATTAAGGCCATCTCACTCAGTACTGCTGCTTTCAAAGTCCACGCTCCCAAAGCAAATTGGATTTCTGTACACAATACTCTTACAGGATGAAACCCAACCAACTTGTTGACGTAAGTATCCCATCTACTTCGCAATTTTATTTATCTTCCAAAATTAAAGGACTGGCACCCTGATTTATTAAAAGTGAATTGGTTCTAGGGACCATATCCCCTCTGAGTTACTGACAGAGCAGCTTCTGGCCTGTGAAGCTCAAAGCCATGCCTAGATGTGAGGATCCCATCAACTATAGCCAGCAGTGGTCTCTGCTCCCAGAGCTTACCTTCAGTTGTAACAAACATCTACAGCTTCACTTTGCTTTCCTTCTCTGTATCTCCTTCACCAATGTGTACATATTCATGTCACATGCTCTTTAACGTTTCTAAGACACAACAATACCCATTATTAATCTCTTACTCAGGCAGTGTTAATTCCAATAAGACTGCAAGGCAGCAGTGCCTCCAGAGCCTGCACTGTGCCGGGAATTGGCACTGGGATTGCAAGCACCGTGGTCATTGCTACCAAGGGAGGCACAGAATCCCTTCACCCCATAGTCACGGGAGCATTGGCAACAAAATTTACCATGCCTTGGCAATTTATGCTGACCCTCACATTTTGGCATTAAAAAAAAGTATCATCTTAGAGTACCAAATATGTGTCTTCTTATAGATGTGTGCCACGGGCCTCGCAGTGAACATCACAGCTCTGCAGATGCACAAACCTTGTCTTTATGAGTCTATATATTAAATAAAA|>contig-10 151 0|TGTGTCTTCTTATAGATGTGTGCCACGGGCCTCGCAGTGAACATCACAGCTCTGCAGATGCACAAACCTTGTCTTTATGAGTCTATATATTAAATAAAAGAGGTAGAATTTTGGATGGGTGCAAGTCAAGTAATGGGATAAAAGAATGATG|>contig-11 192 0|GTGAACATCACAGCTCTGCAGATGCACAAACCTTGTCTTTATGAGTCTATATATTAAATAAAAGTGGTAGAATTTTGGATGGGTGCAAGTCAAGTAATGGGATAAAAGAATGATGATCGAAGGAGATGATTCATGATAACTATGTGCCCTTATGTGAAACCTAGCAACAATTTCTTACTTATGATAAAGGCA|>contig-12 195 0|AAAAAAAGTATCATCTTAGAGTACCAAATATGTGTCTTCTTATAGATGTGTGCCACGGGCCTCGCAGTGAACATCACAGCTCTGCAGATGCACAAACCTTGTCTTTATGAGTCTATATATTAAATAAAAGTGGTAGAATTTTGGATGGGTGCAAGTCAAGTAATGGGATAAAAGAATGATGATCGAAGGAGATGA|>contig-13 165 0|GTGTGCCACGGGCCTCGCAGTGAACATCACAGCTCTGCAGATGCACAAACCTTGTCTTTATGAGTCTATATATTAAATAAAAGTGGTAGAATTTTGGATGGGTGCAAGTCAAGTAATGGGATAAAAGAATGATGATCGAAGGAGATGATTCATGATAACTATGTG|>contig-14 151 0|TAGAATTTTGGATGGGTGCAAGTCAAGTAATGGGATAAAAGAATGATGATCGAAGGAGATGATTCATGATAACTATGTGCCCTTATGTGAAACCTAGCAACAATTTCTTACTTATGATAAAGGCACAATCTCCCCAGTCCAGAATCTAAAT|>contig-15 151 0|ACAAACCTTGTCTTTATGAGTCTATATATTAAATAAAAGTGGTAGAATTTTGGATGCGTGCAAGTCAAGTAATGGGATAAAAGAATGATGATCGAAGGAGATGATTCATGATAACTATGTGCCCTTATGTGAAACCTAGCAACAATTTCTT|>contig-16 169 0|ACATGTTAATTAAGTAAATAGTATAAGCTACACGTTGAGCCTACAACTCTTGCTGGAGTCTCCAGCTGTAAGAGCCTTAACTCCTTCCTTTCATAAAGGCCAAGAAAAGAGTGAGTTAGTCCCTATTCGGGGAGTGTTGGTTGGTTCTCCCAATCCTGTTAGGGGCAGT|>contig-17 151 0|TAGAAGCCAGGCTAAATGAGCCCTAATCAACTAGAACAAGCATTCAGAACTCAAAATCCCTTTACATGTTAATTAAGTAAATAGTATAAGCTACACGTTGAGCCTACAACTCTTGCTGGAGTCTCCAGCTGTAAGAGCCTTAACTCCTTCC";
    private final String[] alignments = {
            "asm000001:tig00000\t1-250%CTG=21START=27375445END=27375694%250M%+%60%0",
            "asm000001:tig00001\t1-151%CTG=21START=27375343END=27375493%151M%+%60%1",
            "asm000001:tig00002\t1-152%CTG=21START=27375074END=27375225%152M%+%60%0",
            "asm000001:tig00003\t1-195%CTG=21START=27375242END=27375436%195M%+%60%0",
            "asm000001:tig00004\t1-151%CTG=21START=27375303END=27375453%151M%+%60%0",
            "asm000001:tig00005\t1-151%CTG=21START=27375198END=27375348%151M%+%60%1",
            "asm000001:tig00006\t1-241%CTG=21START=27375304END=27375544%241M%+%60%0",
            "asm000001:tig00007\t1-151%CTG=21START=27375207END=27375357%151M%+%60%2",
            "asm000001:tig00008\t1-194%CTG=21START=27375399END=27375594%194M%+%60%0",
            "asm000001:tig00009\t1-950%CTG=21START=27373209END=27374158%950M1059S%+%60%0\t944-1491%CTG=21START=27374159END=27374706%943H548M518H%-%60%1\t1492-2009%CTG=21START=27374701END=27375218%1491H518M%+%60%0",
            "asm000001:tig00010\t1-151%CTG=21START=27375120END=27375370%151M%+%60%1",
            "asm000001:tig00011\t1-192%CTG=21START=27375156END=27375347%192M%+%60%0",
            "asm000001:tig00012\t1-195%CTG=21START=27375090END=27375284%195M%+%60%0",
            "asm000001:tig00013\t1-165%CTG=21START=27375137END=27375301%165M%+%60%0",
            "asm000001:tig00014\t1-151%CTG=21START=27375223END=27375373%151M%+%60%0",
            "asm000001:tig00015\t1-151%CTG=21START=27375181END=27375331%151M%+%60%1",
            "asm000001:tig00016\t1-169%CTG=21START=27375550END=27375718%169M%+%60%0",
            "asm000001:tig00017\t1-151%CTG=21START=27375487END=27375637%151M%+%60%1"
    };

    private static final class DiscoverVariantsFromContigAlignmentsSGASparkIntegrationTestArgs {
        final String inputAssemblies;
        final String inputAlignments;
        final String outputDir;

        DiscoverVariantsFromContigAlignmentsSGASparkIntegrationTestArgs(final String inputAssemblies, final String inputAlignments, final String outputDir){
            this.inputAssemblies = inputAssemblies;
            this.inputAlignments = inputAlignments;
            this.outputDir = outputDir;
        }

        String getCommandLineNoApiKey() {
            return  " --inputAssemblies " + inputAssemblies +
                    " --inputAlignments " + inputAlignments +
                    " --fastaReference " + SVIntegrationTestDataProvider.reference +
                    " -R " + SVIntegrationTestDataProvider.reference_2bit +
                    " -O " + outputDir + "/variants.vcf" +
                    " --logContigAlignmentSimpleStats";
        }

    }

    @DataProvider(name = "discoverVariantsFromContigAlignmentsSGASparkIntegrationTest")
    public Object[][] createTestData() throws IOException {
        List<Object[]> tests = new ArrayList<>();

        final File tempWorkingDir = BaseTest.createTempDir("discoverVariantsFromContigAlignmentsSGASparkIntegrationTest");
        tempWorkingDir.deleteOnExit();

        final File assemblyWithPackedFastaWithLength = Files.createDirectory(Paths.get(tempWorkingDir.getAbsolutePath()+"/"+"assemblyWithPackedFastaWithLength")).toFile();
        try(  final PrintWriter out = new PrintWriter( new File(assemblyWithPackedFastaWithLength, "contents.txt") )  ){
            out.println( contigs );
        }
        final File alignmentsDir = Files.createDirectory(Paths.get(tempWorkingDir.getAbsolutePath()+"/"+"alignments")).toFile();
        try(  final PrintWriter out = new PrintWriter( new File(alignmentsDir, "alignments.txt") )  ){
            for(final String aln : alignments) out.println( aln );
        }

        tests.add(new Object[]{new DiscoverVariantsFromContigAlignmentsSGASparkIntegrationTest.DiscoverVariantsFromContigAlignmentsSGASparkIntegrationTestArgs(
                assemblyWithPackedFastaWithLength.getAbsolutePath(),
                alignmentsDir.getAbsolutePath(),
                tempWorkingDir.getAbsolutePath())});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "discoverVariantsFromContigAlignmentsSGASparkIntegrationTest", groups = "sv")
    public void testDiscoverVariantsFromContigAlignmentsSGASparkRunnableLocal(final DiscoverVariantsFromContigAlignmentsSGASparkIntegrationTest.DiscoverVariantsFromContigAlignmentsSGASparkIntegrationTestArgs params) throws Exception {

        final List<String> args = Arrays.asList( new ArgumentsBuilder().add(params.getCommandLineNoApiKey()).getArgsArray() );
        runCommandLine(args);
        StructuralVariationDiscoveryPipelineSparkIntegrationTest.svDiscoveryVCFEquivalenceTest(args.get(args.indexOf("-O")+1), SVIntegrationTestDataProvider.EXPECTED_SIMPLE_INV_VCF, Arrays.asList("ALIGN_LENGTHS", "CTG_NAMES"), false);
        Assert.assertTrue(Files.exists(Paths.get(args.get(args.indexOf("--inputAlignments")+1)+"_withMoreThanTwoAlignments")));
        Assert.assertTrue(Files.exists(Paths.get(args.get(args.indexOf("--inputAlignments")+1)+"_withTwoAlignments")));
    }

    @Test(dataProvider = "discoverVariantsFromContigAlignmentsSGASparkIntegrationTest", groups = "sv")
    public void testDiscoverVariantsRunnableMiniCluster(final DiscoverVariantsFromContigAlignmentsSGASparkIntegrationTest.DiscoverVariantsFromContigAlignmentsSGASparkIntegrationTestArgs params) throws Exception {

        MiniClusterUtils.runOnIsolatedMiniCluster(cluster -> {

            final List<String> argsToBeModified = Arrays.asList( new ArgumentsBuilder().add(params.getCommandLineNoApiKey()).getArgsArray() );
            final Path workingDirectory = MiniClusterUtils.getWorkingDir(cluster);

            int idx = 0;

            idx = argsToBeModified.indexOf("--inputAssemblies");
            Path path = new Path(workingDirectory, "assemblies_0");
            File file = new File(argsToBeModified.get(idx+1));
            cluster.getFileSystem().copyFromLocalFile(new Path(file.toURI()), path);
            argsToBeModified.set(idx+1, path.toUri().toString());

            idx = argsToBeModified.indexOf("--inputAlignments");
            path = new Path(workingDirectory, "alignments");
            file = new File(argsToBeModified.get(idx+1));
            cluster.getFileSystem().copyFromLocalFile(new Path(file.toURI()), path);
            argsToBeModified.set(idx+1, path.toUri().toString());

            idx = argsToBeModified.indexOf("-R");
            path = new Path(workingDirectory, "reference.2bit");
            file = new File(argsToBeModified.get(idx+1));
            cluster.getFileSystem().copyFromLocalFile(new Path(file.toURI()), path);
            argsToBeModified.set(idx+1, path.toUri().toString());

            idx = argsToBeModified.indexOf("--fastaReference");
            path = new Path(workingDirectory, "reference.fasta");
            file = new File(argsToBeModified.get(idx+1));
            cluster.getFileSystem().copyFromLocalFile(new Path(file.toURI()), path);
            argsToBeModified.set(idx+1, path.toUri().toString());

            path = new Path(workingDirectory, "reference.fasta.fai");
            cluster.getFileSystem().copyFromLocalFile(new Path(SVIntegrationTestDataProvider.reference_fai.toURI()), path);
            path = new Path(workingDirectory, "reference.dict");
            cluster.getFileSystem().copyFromLocalFile(new Path(SVIntegrationTestDataProvider.reference_dict.toURI()), path);

            // outputs, prefix with hdfs address
            idx = argsToBeModified.indexOf("-O");
            path = new Path(workingDirectory, "variants.vcf");
            final String vcfOnHDFS = path.toUri().toString();
            argsToBeModified.set(idx+1, vcfOnHDFS);

            runCommandLine(argsToBeModified);
            StructuralVariationDiscoveryPipelineSparkIntegrationTest.svDiscoveryVCFEquivalenceTest(vcfOnHDFS, SVIntegrationTestDataProvider.EXPECTED_SIMPLE_INV_VCF, Arrays.asList("ALIGN_LENGTHS", "CTG_NAMES"), true);

            final RemoteIterator<LocatedFileStatus> it = FileSystem.get(workingDirectory.toUri(), new Configuration()).listFiles(workingDirectory, true);
            while(it.hasNext()) System.err.println(it.next());

            final FileSystem fs = FileSystem.get(workingDirectory.toUri(), new Configuration());
            Assert.assertTrue(fs.exists( new Path(argsToBeModified.get(argsToBeModified.indexOf("--inputAlignments")+1)+"_withTwoAlignments") ));
            Assert.assertTrue(fs.exists( new Path(argsToBeModified.get(argsToBeModified.indexOf("--inputAlignments")+1)+"_withMoreThanTwoAlignments") ));
        });
    }


}
