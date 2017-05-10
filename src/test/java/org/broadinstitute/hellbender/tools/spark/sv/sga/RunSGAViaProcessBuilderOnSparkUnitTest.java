package org.broadinstitute.hellbender.tools.spark.sv.sga;

import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.SkipException;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;
import scala.Tuple2;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;

/**
 * To make this test runnable, "-Dpath.to.sga" should be set properly.
 */
// TODO: investigate why SGA, run with exactly the same parameters on the same machine, yields slightly different results: either SNP between the contigs or slightly different contig lengths
// TODO: test if error messages are handled correctly
public class RunSGAViaProcessBuilderOnSparkUnitTest extends CommandLineProgramTest {

    private static Path sgaPath;

    private static final SGAModule indexer = new SGAModule("index");

    // data block
    private static File TEST_DATA_DIR;
    private static List<Tuple2<Long, File>> rawFASTQFiles;
    private static Map<Long, File> expectedAssembledFASTAFiles;

    // parameter block
    private static final int lengthDiffTolerance = 1;   // tolerance on how much of a difference can each pair of assembled contigs (actual vs expected) length could be
    private static final int editDistanceTolerance = 3; // tolerance of the summed edit distance between expected contigs and corresponding actual assembled contigs

    // set up path and data
    @BeforeClass
    private static void setup() throws IOException{
        String sgaPathString = System.getProperty("path.to.sga");
        if(null==sgaPathString){
            sgaPathString = System.getenv("path.to.sga");
        }
        sgaPath =  (null==sgaPathString) ? null : Paths.get(sgaPathString);

        if(null==sgaPath){ // early return if SGA isn't available for these tests
            return;
        }

        rawFASTQFiles = new ArrayList<>();
        expectedAssembledFASTAFiles = new HashMap<>();

        TEST_DATA_DIR = new File(getTestDataDir(), "spark/sv/RunSGAViaProcessBuilderOnSpark/");
        final File tempWorkingDir = BaseTest.createTempDir("whatever");

        Long i = 0L;

        final String[] extension = new String[1]; // look for FASTQ files in test data dir
        extension[0] = "fastq";

        final Iterator<File> it = FileUtils.iterateFiles(TEST_DATA_DIR, extension, false);

        while(it.hasNext()){

            final File fastqFile = it.next();
            FileUtils.copyFileToDirectory(fastqFile, tempWorkingDir);
            final File tempFASTQFile = new File(tempWorkingDir, fastqFile.getName());
            rawFASTQFiles.add(new Tuple2<>(i, tempFASTQFile));

            final File contigFile = new File(TEST_DATA_DIR, fastqFile.getName().replace("fastq", "pp.ec.filter.pass.merged.rmdup-contigs.fa"));
            expectedAssembledFASTAFiles.put(i++, contigFile);
        }
    }

    @BeforeMethod
    private static void checkPathToSGAAvailability(){
        if(null==sgaPath){
            throw new SkipException("Skipping test because \"-Dpath.to.sga\" is set neither in system property nor as an environment variable.");
        }
    }

    @Test(groups = "sv")
    public void writeToLocalTest() throws IOException{
        final File originalFASTQFile = rawFASTQFiles.get(0)._2();
        final String pathToFile = originalFASTQFile.toURI().toString();
        final String fileContents = StringUtils.join(FileUtils.readLines(originalFASTQFile, StandardCharsets.UTF_8), "\n");
        final Tuple2<Long, File> copied = RunSGAViaProcessBuilderOnSpark.writeToLocal(new Tuple2<>(pathToFile, fileContents), ".raw");
        Assert.assertTrue( FileUtils.contentEqualsIgnoreEOL(copied._2(), originalFASTQFile, StandardCharsets.UTF_8.name()) );
    }

    // Test if disabling stdout stderr works
    @Test(groups = "sv")
    public void stdioCaptureTest() throws IOException{

        final File dir = BaseTest.createTempDir("xyz");

        FileUtils.copyFileToDirectory(rawFASTQFiles.get(0)._2(), dir);

        final String fileName = rawFASTQFiles.get(0)._2().getName();

        File rawFASTQFile = new File(dir, fileName);

        final List<SGAModule.RuntimeInfo> runtimeInfo = new ArrayList<>();

        final ArrayList<String> indexerArgs = new ArrayList<>(4);
        indexerArgs.add("--algorithm"); indexerArgs.add("ropebwt");
        indexerArgs.add("--check");
        indexerArgs.add(fileName);

        final ArrayList<String> moduleArgs = new ArrayList<>(6);
        moduleArgs.add("--pe-mode");    moduleArgs.add("2");
        moduleArgs.add("--out");        moduleArgs.add(FilenameUtils.getBaseName(rawFASTQFile.getName()) + ".pp.fa");
        moduleArgs.add("--permute-ambiguous");

        RunSGAViaProcessBuilderOnSpark.runSGAModule(sgaPath, "preprocess", ".pp.fa", rawFASTQFile, dir, indexer, indexerArgs, moduleArgs, runtimeInfo, false);

        for(final SGAModule.RuntimeInfo info : runtimeInfo){
            Assert.assertTrue(info.stdoutMsg.isEmpty() && info.stderrMsg.isEmpty());
        }

        runtimeInfo.clear();

        // deliberately let the process fail and see if stdout stderr is captured
        // sga preprocess produces runtime metrics to stderr, and
        // sga index is provided with wrong arguments this time, so it should also fill stderr with message
        indexerArgs.set(1, "wrongAlgo");
        RunSGAViaProcessBuilderOnSpark.runSGAModule(sgaPath, "preprocess", ".pp.fa", rawFASTQFile, dir, indexer, indexerArgs, moduleArgs, runtimeInfo, true);
        for(final SGAModule.RuntimeInfo info : runtimeInfo){
            Assert.assertFalse(info.stderrMsg.isEmpty());
        }

        FileUtils.cleanDirectory(dir);
    }

    @Test(groups = "sv")
    public void mapToListTest(){
        Map<String, String> indexerArgsMap = new LinkedHashMap<>();

        List<String> convertedList = RunSGAViaProcessBuilderOnSpark.turnCmdLineArgsKeyValuePairIntoList(indexerArgsMap);
        Assert.assertTrue(convertedList.isEmpty());

        indexerArgsMap.put("--algorithm", "ropebwt");
        indexerArgsMap.put("--check", "");

        convertedList = RunSGAViaProcessBuilderOnSpark.turnCmdLineArgsKeyValuePairIntoList(indexerArgsMap);
        Assert.assertEquals(Arrays.asList("--algorithm", "ropebwt", "--check"), convertedList);
    }

    @Test(groups = "sv")
    public void assemblyOneStepTest() throws IOException{

        for(final Tuple2<Long, File> breakpoint : rawFASTQFiles){

            final RunSGAViaProcessBuilderOnSpark.SGAAssemblyResult result = RunSGAViaProcessBuilderOnSpark.runSGAModulesInSerial(sgaPath.toAbsolutePath().toString(), breakpoint._2(), true, false);

            Assert.assertTrue(result.assembledContigs!=null);

            final List<String> actualAssembledContigs = result.assembledContigs.toListOfStrings();
            final List<String> expectedAssembledContigs = FileUtils.readLines(expectedAssembledFASTAFiles.get(breakpoint._1()), StandardCharsets.UTF_8);

            compareContigs(actualAssembledContigs, expectedAssembledContigs);
        }
    }

    @Test(groups = "sv")
    public void assemblyStepByStepTest() throws IOException{

        for(final Tuple2<Long, File> breakpoint : rawFASTQFiles){

            final File expectedAssembledContigsFile = expectedAssembledFASTAFiles.get(breakpoint._1());

            stepByStepTestWorker(breakpoint._2(), expectedAssembledContigsFile);
        }
    }

    private static void stepByStepTestWorker(final File rawFASTQFile, final File expectedAssembledContigsFile) throws IOException{

        final File workingDir = rawFASTQFile.getParentFile();
        final List<SGAModule.RuntimeInfo> runtimeInfo = new ArrayList<>();

        ArrayList<Integer> editDistancesBetweenSeq = new ArrayList<>();
        final Integer zero = 0;

        final ArrayList<String> indexerArgs = new ArrayList<>();
        indexerArgs.add("--algorithm"); indexerArgs.add("ropebwt");
        indexerArgs.add("--check");
        indexerArgs.add("");

        final List<String> moduleArgs = new ArrayList<>(20);

        final String filenamePrefix = FilenameUtils.getBaseName(rawFASTQFile.getName());
        moduleArgs.clear();
        moduleArgs.add("--pe-mode");    moduleArgs.add("2");
        moduleArgs.add("--out");        moduleArgs.add(FilenameUtils.getBaseName(rawFASTQFile.getName()) + ".pp.fa");

        final File actualPreppedFile = new File(workingDir, RunSGAViaProcessBuilderOnSpark.runSGAModule(sgaPath, "preprocess", ".pp.fa", rawFASTQFile, workingDir, indexer, indexerArgs, moduleArgs, runtimeInfo, false));
        final File expectedPreppedFile = new File(TEST_DATA_DIR, filenamePrefix + ".pp.fq");
        final String preppedFileName = compareNamesAndComputeSeqEditDist(actualPreppedFile, expectedPreppedFile, true, editDistancesBetweenSeq);
        for(final Integer d : editDistancesBetweenSeq){ Assert.assertEquals(d, zero); }

        editDistancesBetweenSeq = new ArrayList<>();
        moduleArgs.clear();
        final File preppedFile = new File(workingDir, preppedFileName);
        final File actualCorrectedFile = new File(workingDir, RunSGAViaProcessBuilderOnSpark.runSGAModule(sgaPath, "correct", ".ec.fa", preppedFile, workingDir, indexer, indexerArgs, moduleArgs, runtimeInfo, false));
        final File expectedCorrectedFile = new File(TEST_DATA_DIR, filenamePrefix + ".pp.ec.fa");
        final String correctedFileName = compareNamesAndComputeSeqEditDist(actualCorrectedFile, expectedCorrectedFile, true, editDistancesBetweenSeq);
        for(final Integer d : editDistancesBetweenSeq){ Assert.assertEquals(d, zero); }

        editDistancesBetweenSeq = new ArrayList<>();
        moduleArgs.clear();
        moduleArgs.add("--kmer-threshold");         moduleArgs.add(String.valueOf(RunSGAViaProcessBuilderOnSpark.FILTER_STEP_KMER_FREQUENCY_THREASHOLD));
        moduleArgs.add("--homopolymer-check");
        moduleArgs.add("--low-complexity-check");
        final File correctedFile = new File(workingDir, correctedFileName);
        final File actualFilterPassingFile = new File(workingDir, RunSGAViaProcessBuilderOnSpark.runSGAModule(sgaPath, "filter", ".filter.pass.fa", correctedFile, workingDir, null, null, moduleArgs, runtimeInfo, false));
        final File expectedFilterPassingFile = new File(TEST_DATA_DIR, filenamePrefix + ".pp.ec.filter.pass.fa");
        final String filterPassingFileName = compareNamesAndComputeSeqEditDist(actualFilterPassingFile, expectedFilterPassingFile, true, editDistancesBetweenSeq);
        for(final Integer d : editDistancesBetweenSeq){ Assert.assertEquals(d, zero); }

        editDistancesBetweenSeq = new ArrayList<>();
        moduleArgs.clear();
        moduleArgs.add("--min-overlap");     moduleArgs.add(String.valueOf(RunSGAViaProcessBuilderOnSpark.MIN_OVERLAP_IN_FILTER_OVERLAP_ASSEMBLE));
        final File filterPassingFile = new File(workingDir, filterPassingFileName);
        final File actualMergedFile = new File(workingDir, RunSGAViaProcessBuilderOnSpark.runSGAModule(sgaPath, "fm-merge", ".merged.fa", filterPassingFile, workingDir, indexer, indexerArgs, moduleArgs, runtimeInfo, false));
        final File expectedMergedFile = new File(TEST_DATA_DIR, filenamePrefix + ".pp.ec.filter.pass.merged.fa");
        final String mergedFileName = compareNamesAndComputeSeqEditDist(actualMergedFile, expectedMergedFile, false, editDistancesBetweenSeq);
        for(final Integer d : editDistancesBetweenSeq){ Assert.assertEquals(d, zero); }

        editDistancesBetweenSeq = new ArrayList<>();
        moduleArgs.clear();
        final File mergedFile = new File(workingDir, mergedFileName);
        final File actualDeduplicatedFile = new File(workingDir, RunSGAViaProcessBuilderOnSpark.runSGAModule(sgaPath, "rmdup", ".rmdup.fa", mergedFile, workingDir, indexer, indexerArgs, moduleArgs, runtimeInfo, false));
        final File expectedDeduplicatedFilie = new File(TEST_DATA_DIR, filenamePrefix + ".pp.ec.filter.pass.merged.rmdup.fa");
        final String readyForOverlappingFileName = compareNamesAndComputeSeqEditDist(actualDeduplicatedFile, expectedDeduplicatedFilie, false, editDistancesBetweenSeq);
        for(final Integer d : editDistancesBetweenSeq){ Assert.assertEquals(d, zero); }

        // do overlapping but don't test output
        moduleArgs.clear();
        moduleArgs.add("--min-overlap");   moduleArgs.add(String.valueOf(RunSGAViaProcessBuilderOnSpark.MIN_OVERLAP_IN_FILTER_OVERLAP_ASSEMBLE));
        moduleArgs.add("--error-rate");    moduleArgs.add(String.valueOf(RunSGAViaProcessBuilderOnSpark.OVERLAP_STEP_ERROR_RATE));
        final File readyForOverlappingFile = new File(workingDir, readyForOverlappingFileName);
        final File readyForAssemblyFile = new File(workingDir, RunSGAViaProcessBuilderOnSpark.runSGAModule(sgaPath, "overlap", ".asqg.gz", readyForOverlappingFile, workingDir, null, null, moduleArgs, runtimeInfo, false));

        moduleArgs.clear();
        moduleArgs.add("--cut-terminal");         moduleArgs.add(String.valueOf(RunSGAViaProcessBuilderOnSpark.CUT_OFF_TERMINAL_BRANCHES_IN_N_ROUNDS));
        moduleArgs.add("--max-gap-divergence");   moduleArgs.add(String.valueOf(RunSGAViaProcessBuilderOnSpark.ASSEMBLE_STEP_MAX_GAP_DIVERGENCE));
        moduleArgs.add("--min-overlap");          moduleArgs.add(String.valueOf(RunSGAViaProcessBuilderOnSpark.MIN_OVERLAP_IN_FILTER_OVERLAP_ASSEMBLE));
        moduleArgs.add("--min-branch-length");    moduleArgs.add(String.valueOf(RunSGAViaProcessBuilderOnSpark.ASSEMBLE_STEP_MIN_BRANCH_TAIL_LENGTH));
        moduleArgs.add("--out-prefix");           moduleArgs.add(FilenameUtils.getBaseName(readyForAssemblyFile.getName()).replace(".asqg", ""));
        if(RunSGAViaProcessBuilderOnSpark.TURN_ON_TRANSITIVE_REDUCTION){
            moduleArgs.add("--transitive-reduction");
        }

        // final assembled contig test
        final String assembledFileNameTobeTreated = RunSGAViaProcessBuilderOnSpark.runSGAModule(sgaPath, "assemble", "-contigs.fa", readyForAssemblyFile, workingDir, null, null, moduleArgs, runtimeInfo, false);
        final File actualAssembledContigsFile = new File(workingDir, assembledFileNameTobeTreated.replace(".asqg", ""));
        final List<String> actualAssembledContigs   = (null==actualAssembledContigsFile  ) ? null : Files.readAllLines(Paths.get(actualAssembledContigsFile.getAbsolutePath()));
        final List<String> expectedAssembledContigs = (null==expectedAssembledContigsFile) ? null : Files.readAllLines(Paths.get(expectedAssembledContigsFile.getAbsolutePath()));

        compareContigs(actualAssembledContigs, expectedAssembledContigs);
    }

    private static void compareContigs(final List<String> actualFASTAContents,
                                       final List<String> expectedFASTAContents){

        // first create mapping from read length to sequence
        // (bad, but contig names are not guaranteed to be reproducible in different runs)
        final SortedMap<Integer, String> actualMap = collectContigsByLength(actualFASTAContents);

        final SortedMap<Integer, String> expectedMap = collectContigsByLength(expectedFASTAContents);

        // essentially, make sure the number of assembled contigs with unique length values are the same
        final SortedSet<Integer> actualLengthVals = new TreeSet<>(actualMap.keySet());
        final SortedSet<Integer> expectedLengthVals = new TreeSet<>(expectedMap.keySet());
        Assert.assertEquals(actualLengthVals.size(), expectedLengthVals.size());

        // first compare contig length values, with some tolerance
        final Iterator<Integer> itActual = actualLengthVals.iterator();
        final Iterator<Integer> itExpected = expectedLengthVals.iterator();
        // if the two sets are not exactly the same but difference is within tolerance, create new map
        final List<Tuple2<String, String>> actualVSexpectedPairs = new ArrayList<>();
        while(itActual.hasNext()){
            final Integer actualLength = itActual.next();
            final Integer expectedLength = itExpected.next(); // safe operation as the previous assertion guarantees two sets are of same size
            Assert.assertTrue( Math.abs(actualLength - expectedLength)<=lengthDiffTolerance );
            actualVSexpectedPairs.add(new Tuple2<>(actualMap.get(actualLength), expectedMap.get(expectedLength)));
        }

        // then compare sequences, aware of RC
        for(final Tuple2<String, String> pairOfSimilarLength: actualVSexpectedPairs){

            final String actualString = pairOfSimilarLength._1();
            final String expectedString = pairOfSimilarLength._2();

            final String rcOfActualString = new String(BaseUtils.simpleReverseComplement(actualString.getBytes(StandardCharsets.UTF_8)));

            final int dist = StringUtils.getLevenshteinDistance(actualString, expectedString);
            final int rcDist = StringUtils.getLevenshteinDistance(rcOfActualString, expectedString);

            boolean sequencesAreTheSame = false;
            boolean sequencesAreCloseEnough = false;

            final int minDist = (dist < rcDist ? dist : rcDist);

            if(actualString.equals(expectedString) || rcOfActualString.equals(expectedString)){
                sequencesAreCloseEnough = sequencesAreTheSame = true;
            }else{
                // two things: the unit test tests if sequences are similar enough to each other
                // but if the minimum edit distance is different, we want to see which one it is.
                sequencesAreCloseEnough = (minDist <= editDistanceTolerance);
                if(0!=minDist){
                    System.err.println("Contig that has nonzero edit distance is of length " + Integer.toString(actualString.length()) +
                                        "\nEdit distance:" + String.valueOf(minDist) +
                                        "\nactual sequence: " + actualString +
                                        "\nreverse complement of actual sequence: " + rcOfActualString +
                                        "\nexpected sequence: " + expectedString);
                }
            }

            Assert.assertTrue(sequencesAreTheSame || sequencesAreCloseEnough,
                              "Culprit contigs (expected, actual) of length: " + String.valueOf(expectedString.length()) + ", " + String.valueOf(actualString.length()) +
                              "With edit distance " + String.valueOf(minDist));
        }
    }

    private static String compareNamesAndComputeSeqEditDist(final File actualFile, final File expectedFile, final boolean fastqFilesWellFormed, final ArrayList<Integer> distances) throws IOException{

        final List<String> actualNames = new ArrayList<>();
        final List<String> actualSeq = new ArrayList<>();
        extractNamesAndSeqFromFASTQ(actualFile, fastqFilesWellFormed, actualNames, actualSeq);

        final List<String> expectedNames = new ArrayList<>();
        final List<String> expectedSeq = new ArrayList<>();
        extractNamesAndSeqFromFASTQ(expectedFile, fastqFilesWellFormed, expectedNames, expectedSeq);

        Assert.assertEquals(expectedNames, actualNames);

        // in stead of more stringent tests, e.g. "Assert.assertEquals(expectedSeq, actualSeq);",
        // compute distance and let caller decide how large an edit distance is allowed
        for(int i=0; i<expectedSeq.size(); ++i){
            distances.add(StringUtils.getLevenshteinDistance(expectedSeq.get(i), actualSeq.get(i)));
        }

        //return RunSGAViaProcessBuilderOnSpark.extractBaseNameWithoutExtension(actualFile) + ".fa";
        return FilenameUtils.getBaseName(actualFile.getName()) + ".fa";
    }

    // utility function: for extracting read names and sequences from a fastq/fasta file
    private static void extractNamesAndSeqFromFASTQ(final File FASTAFile, final boolean fastqFilesWellFormed, List<String> readNames, List<String> sequences) throws IOException{

        if(fastqFilesWellFormed){
            try(final FastqReader reader = new FastqReader(FASTAFile)){
                while(reader.hasNext()){
                    final FastqRecord record = reader.next();
                    readNames.add(record.getReadName());
                    sequences.add(record.getReadString());
                }
            }
        }else{
            try(final BufferedReader reader = new BufferedReader(new FileReader(FASTAFile))){
                String line = "";
                int i=0;
                while ((line = reader.readLine()) != null){
                    final int l = i%4;
                    if(0==l){
                        readNames.add(line);
                        ++i;
                    }else{
                        sequences.add(line);
                        reader.readLine(); reader.readLine();
                        i+=3;
                    }
                }
            }
        }
    }

    // utility function: generate mapping from contig length to its DNA sequence
    // this is possible because test result contigs have unique length values
    private static SortedMap<Integer, String> collectContigsByLength(final List<String> contigFileContents){

        final Iterator<String> it = contigFileContents.iterator();

        final SortedMap<Integer, String> result = new TreeMap<>();
        int i=0;
        while(it.hasNext()){
            final String line = it.next();
            if((++i)%2==1) {
                continue;
            } else{
                result.put(line.length(), line);
            }
        }
        return result;
    }
}