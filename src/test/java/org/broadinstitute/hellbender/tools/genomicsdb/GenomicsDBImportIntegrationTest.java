package org.broadinstitute.hellbender.tools.genomicsdb;

import com.fasterxml.jackson.core.type.TypeReference;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.intel.genomicsdb.GenomicsDBUtils;
import com.intel.genomicsdb.model.GenomicsDBExportConfiguration;
import com.intel.genomicsdb.model.GenomicsDBVidMapProto;
import com.intel.genomicsdb.reader.GenomicsDBFeatureReader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.PositionalBufferedStream;
import htsjdk.variant.bcf2.BCF2Codec;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.Main;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.BaseTest;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.*;
import java.util.*;
import java.util.concurrent.CompletionException;
import java.util.stream.Collectors;
import java.util.stream.Stream;

@Test(groups = {"variantcalling"})
public final class GenomicsDBImportIntegrationTest extends CommandLineProgramTest {
    private static final String HG_00096 = largeFileTestDir + "gvcfs/HG00096.g.vcf.gz";
    private static final String HG_00268 = largeFileTestDir + "gvcfs/HG00268.g.vcf.gz";
    private static final String NA_19625 = largeFileTestDir + "gvcfs/NA19625.g.vcf.gz";
    //The following 3 files were obtained by running CombineGVCFs on the above 3 files (separately). This introduces spanning
    //deletions in the files. Hence, these files can be used to test for spanning deletions in the input VCF.
    private static final String HG_00096_after_combine_gvcfs = largeFileTestDir + "gvcfs/HG00096_after_combine_gvcfs.g.vcf.gz";
    private static final String HG_00268_after_combine_gvcfs = largeFileTestDir + "gvcfs/HG00268_after_combine_gvcfs.g.vcf.gz";
    private static final String NA_19625_after_combine_gvcfs = largeFileTestDir + "gvcfs/NA19625_after_combine_gvcfs.g.vcf.gz";
    private static final String NA_24385 = largeFileTestDir + "NA24385.vcf.gz";
    private static final String NA_12878_PHASED = largeFileTestDir + "NA12878.phasedData.Chr20.vcf"; //NOTE: this is not phased according to the vcf spec but it reflects phasing currently produced by haplotype caller
    private static final String MULTIPLOID_DATA_HG37 = largeFileTestDir + "gvcfs/HapMap5plex.ploidy10.b37.g.vcf";
    private static final String NA12878_HG37 = toolsTestDir + "haplotypecaller/expected.testGVCFMode.gatk4.g.vcf";
    private static final String MULTIPLOID_EXPECTED_RESULT = toolsTestDir + "GenomicsDBImport/expected.testGenomicsDBImportWithNonDiploidData.vcf";
    private static final String MNP_GVCF = toolsTestDir + "GenomicsDBImport/mnp.input.g.vcf";
    private static final String ARTIFICIAL_PHASED = getTestDataDir() + "/ArtificalPhasedData.1.g.vcf";
    private static final String HG_00268_WITH_SPACES = largeFileTestDir + "gvcfs/HG00268.spaceInSampleName.g.vcf";
    private static final List<String> LOCAL_GVCFS = Arrays.asList(HG_00096, HG_00268, NA_19625);
    private static final List<String> LOCAL_GVCFS_AFTER_COMBINE_GVCFS = Arrays.asList(HG_00096_after_combine_gvcfs,
            HG_00268_after_combine_gvcfs,
            NA_19625_after_combine_gvcfs);
    private static final String GENOMICSDB_TEST_DIR = toolsTestDir + "GenomicsDBImport/";
    private static final String COMBINEGVCFS_TEST_DIR = toolsTestDir + "walkers/CombineGVCFs/";
    private static final String COMBINED = largeFileTestDir + "gvcfs/combined.gatk3.7.g.vcf.gz";
    private static final String COMBINED_WITH_GENOTYPES = largeFileTestDir + "gvcfs/combined_with_genotypes.g.vcf.gz";
    //This file was obtained from combined.gatk3.7.g.vcf.gz by dropping all the samples
    private static final String COMBINED_SITES_ONLY = largeFileTestDir + "gvcfs/combined.gatk3.7_sites_only.g.vcf.gz";
    //Consider a gVCF with a REF block chr20:50-150. Importing this data into GenomicsDB using multiple intervals
    //-L chr20:1-100 and -L chr20:101-200 will cause the REF block to be imported into both the arrays
    //Now, when reading data from the workspace (assume full scan) - the data is split into 2 REF block intervals chr20:50-100
    //and chr20:101-150 one from each array
    //The following COMBINED_MULTI_INTERVAL gvcf is identical to the gVCF in the previous line except at the partition break
    //position
    //The previous file has the following line:
    //chr20   17970000        .       G       <NON_REF>       .       .       END=17970001
    //
    //while this file has:
    //chr20   17970000        .       G       <NON_REF>       .       .       .
    //chr20   17970001        .       G       <NON_REF>       .       .       .
    //
    private static final String COMBINED_MULTI_INTERVAL = largeFileTestDir + "gvcfs/combined_multi_interval.gatk3.7.g.vcf.gz";
    private static final String COMBINED_WITHSPACES = largeFileTestDir + "gvcfs/combined.gatk3.7.smaller_interval.g.vcf";
    private static final ArrayList<SimpleInterval> INTERVAL =
            new ArrayList<SimpleInterval>(Arrays.asList(new SimpleInterval("chr20", 17960187, 17981445)));
    private static final ArrayList<SimpleInterval> MULTIPLE_INTERVALS = new ArrayList<SimpleInterval>(Arrays.asList(
        new SimpleInterval("chr20", 17960187, 17970000),
        new SimpleInterval("chr20", 17970001, 17980000),
        new SimpleInterval("chr20", 17980001, 17981445)
    ));
    private static final ArrayList<SimpleInterval> MULTIPLE_INTERVALS_THAT_WORK_WITH_COMBINE_GVCFS =
        new ArrayList<SimpleInterval>(Arrays.asList(
            new SimpleInterval("chr20", 17960187, 17969999),
            new SimpleInterval("chr20", 17970000, 17980000),
            new SimpleInterval("chr20", 17980001, 17981445)
    ));
    private static final ArrayList<SimpleInterval> MULTIPLE_NON_ADJACENT_INTERVALS_THAT_WORK_WITH_COMBINE_GVCFS =
        new ArrayList<SimpleInterval>(Arrays.asList(
            new SimpleInterval("chr20", 17960187, 17969999),
            new SimpleInterval("chr20", 17980001, 17981445),
            new SimpleInterval("chr21", 29477554, 29486255)
    ));
    private static final ArrayList<SimpleInterval> INTERVAL_3736 =
            new ArrayList<SimpleInterval>(Arrays.asList(new SimpleInterval("chr6",130365070,146544250)));
    private static final ArrayList<SimpleInterval> INTERVAL_NONDIPLOID =
            new ArrayList<SimpleInterval>(Arrays.asList(new SimpleInterval("20", 10000000, 10100000)));
    private static final ArrayList<SimpleInterval> SMALLER_INTERVAL =
            new ArrayList<SimpleInterval>(Arrays.asList(new SimpleInterval("chr20", 17960187, 17961973)));
    private static final VCFHeader VCF_HEADER = VariantContextTestUtils.getCompleteHeader();
    private static final String SAMPLE_NAME_KEY = "SN";
    private static final String ANOTHER_ATTRIBUTE_KEY = "AA";

    private static final List<String> GVCFS_WITH_NEW_MQ = Arrays.asList(toolsTestDir + "/haplotypecaller/expected.testGVCFMode.gatk4.g.vcf", getTestDataDir() + "/walkers/CombineGVCFs/YRIoffspring.chr20snippet.g.vcf");
    private static final String COMBINED_WITH_NEW_MQ = toolsTestDir + "/walkers/GenomicsDBImport/newMQcalc.combined.g.vcf";
    private static final List<SimpleInterval> INTERVAL2 = Arrays.asList(new SimpleInterval("20", 1, 11_000_000));
    private static final List<String> ATTRIBUTES_TO_IGNORE = Arrays.asList("RAW_MQ","RAW_MQandDP");  //CombineGVCFs doesn't support the old RAW_MQ anymore

    @Override
    public String getTestedClassName() {
        return GenomicsDBImport.class.getSimpleName();
    }

    @DataProvider(name="batchSizes")
    public Object[][] batchSizes() {
        return new Object[][] {
                new Object[]{1},
                new Object[]{2},
                new Object[]{3},
                new Object[]{4},
                new Object[]{100},
        };
    }

    @Test
    public void testGenomicsDBImportFileInputs() throws IOException {
        testGenomicsDBImporter(LOCAL_GVCFS, INTERVAL, COMBINED, b38_reference_20_21, true, 1);
    }

    @Test
    public void testGenomicsDBImportFileInputs_newMQ() throws IOException {
        testGenomicsDBImporter_newMQ(GVCFS_WITH_NEW_MQ, INTERVAL2, COMBINED_WITH_NEW_MQ, b37_reference_20_21, true, Collections.emptyList());
    }

    @Test
    public void testGenomicsDBImportFileInputsWithMultipleIntervals() throws IOException {
        testGenomicsDBImporter(LOCAL_GVCFS, MULTIPLE_INTERVALS, COMBINED_MULTI_INTERVAL, b38_reference_20_21, true, 1);
    }

    private void testGenomicsDBImportWith1000Intervals() throws IOException {
        final String workspace = createTempDir("genomicsdb-tests-").getAbsolutePath() + "/workspace";
        LinkedList<SimpleInterval> intervals = new LinkedList<SimpleInterval>();
        //[ 17960187, 17981445 ]
        int base = 17960187;
        for (int i = 0; i < 1000; ++i)
            intervals.add(new SimpleInterval("chr20", base + 20 * i, base + 20 * i + 10)); //intervals of size 10 separated by 10
        writeToGenomicsDB(new ArrayList<String>(Arrays.asList(LOCAL_GVCFS.get(0))), intervals, workspace, 0, false, 0, 1);
    }

    @Test
    public void testGenomicsDBImportFileInputsAgainstCombineGVCF() throws IOException {
        testGenomicsDBAgainstCombineGVCFs(LOCAL_GVCFS, INTERVAL, b38_reference_20_21, new String[0]);
    }

    @Test
    public void testGenomicsDBImportFileInputsAgainstCombineGVCFWithMultipleIntervals() throws IOException {
        testGenomicsDBAgainstCombineGVCFs(LOCAL_GVCFS, MULTIPLE_INTERVALS_THAT_WORK_WITH_COMBINE_GVCFS, b38_reference_20_21, new String[0]);
    }

    @Test
    public void testGenomicsDBImportFileInputsAgainstCombineGVCFWithMultipleNonAdjacentIntervals() throws IOException {
        testGenomicsDBAgainstCombineGVCFs(LOCAL_GVCFS, MULTIPLE_NON_ADJACENT_INTERVALS_THAT_WORK_WITH_COMBINE_GVCFS,
            b38_reference_20_21, new String[0]);
    }

    @Test
    public void testGenomicsDBImportFileInputsAgainstCombineGVCFWithMultipleNonAdjacentIntervalsForFilesProducedAfterCombineGVCFs()
        throws IOException {
        //this test covers the scenario where the input vcfs have spanning deletions
        testGenomicsDBAgainstCombineGVCFs(LOCAL_GVCFS_AFTER_COMBINE_GVCFS, MULTIPLE_NON_ADJACENT_INTERVALS_THAT_WORK_WITH_COMBINE_GVCFS,
            b38_reference_20_21, new String[0]);
    }

    @Test
    public void testGenomicsDBImportFileInputsAgainstCombineGVCFWithNonDiploidData() throws IOException {
        testGenomicsDBImporterWithGenotypes(Arrays.asList(NA12878_HG37, MULTIPLOID_DATA_HG37), INTERVAL_NONDIPLOID,
                MULTIPLOID_EXPECTED_RESULT, b37_reference_20_21,
                true,
                false,
                false);
    }

    @Test
    public void testGenomicsDBImportPhasedData() throws IOException {
        testGenomicsDBImporterWithGenotypes(Arrays.asList(NA_12878_PHASED), INTERVAL, NA_12878_PHASED, b37_reference_20_21);
    }

    @Test
    public void testGenomicsDBImportPhasedDataWithMultipleIntervals() throws IOException {
        testGenomicsDBImporterWithGenotypes(Arrays.asList(NA_12878_PHASED), MULTIPLE_INTERVALS, NA_12878_PHASED, b37_reference_20_21);
    }

    @Test
    public void testGenomicsDBImportArtificialPhasedData() throws IOException {
        ArrayList<SimpleInterval> intervals = new ArrayList<SimpleInterval>(Arrays.asList(new SimpleInterval("1", 10109, 10297)));
        testGenomicsDBImporterWithGenotypes(Arrays.asList(ARTIFICIAL_PHASED), intervals, ARTIFICIAL_PHASED, b37_reference_20_21);
    }

    @Test
    public void testGenomicsDBThreeLargeSamplesWithGenotypes() throws IOException {
        ArrayList<SimpleInterval> intervals = new ArrayList<SimpleInterval>(Arrays.asList(new SimpleInterval("chr20", 1, 64444167)));
        testGenomicsDBImporterWithGenotypes(LOCAL_GVCFS, intervals, COMBINED_WITH_GENOTYPES, b38_reference_20_21, true, true, false);
    }

    @Test
    public void testGenomicsDBThreeLargeSamplesSitesOnlyQuery() throws IOException {
        ArrayList<SimpleInterval> intervals = new ArrayList<SimpleInterval>(Arrays.asList(
                    new SimpleInterval("chr20", 1, 64444167),
                    new SimpleInterval("chr21", 1, 46709983)));
        testGenomicsDBImporterWithGenotypes(LOCAL_GVCFS, intervals, COMBINED_SITES_ONLY, b38_reference_20_21, true, true, true);
    }

    @Test(expectedExceptions={UserException.BadInput.class}, expectedExceptionsMessageRegExp=".*GenomicsDBImport does not support GVCFs.*")
    public void testGenomicsDbImportThrowsOnMnp() throws IOException {
        for (int threads = 1; threads <= 2; ++threads) {
            testGenomicsDBImporter(
                    Collections.singletonList(MNP_GVCF),
                    Collections.singletonList(new SimpleInterval("20", 69700, 69900)),
                    null, // Should never produce a VCF
                    b38_reference_20_21,
                    true,
                    threads
            );
        }
    }

    private void testGenomicsDBImporterWithGenotypes(final List<String> vcfInputs, final List<SimpleInterval> intervals,
                                                     final String expectedCombinedVCF,
                                                      final String referenceFile) throws IOException {
        testGenomicsDBImporterWithGenotypes(vcfInputs, intervals,
                expectedCombinedVCF, referenceFile,
                false,
                true,
                false);
    }

    private void testGenomicsDBImporterWithGenotypes(final List<String> vcfInputs, final List<SimpleInterval> intervals,
                                                      final String expectedCombinedVCF, final String referenceFile,
                                                     final boolean testAll) throws IOException {
        testGenomicsDBImporterWithGenotypes(vcfInputs, intervals,
                expectedCombinedVCF, referenceFile,
                testAll,
                false,
                false);
    }

    private void testGenomicsDBImporterWithGenotypes(final List<String> vcfInputs, final List<SimpleInterval> intervals,
                                                      final String expectedCombinedVCF, final String referenceFile,
                                                     final boolean testAll,
                                                     final boolean produceGTField,
                                                     final boolean sitesOnlyQuery) throws IOException {
        final String workspace = createTempDir("genomicsdb-tests-").getAbsolutePath() + "/workspace";

        writeToGenomicsDB(vcfInputs, intervals, workspace, 0, false, 0, 1);
        checkJSONFilesAreWritten(workspace);
        checkGenomicsDBAgainstExpected(workspace, intervals, expectedCombinedVCF, referenceFile, testAll, ATTRIBUTES_TO_IGNORE, produceGTField, sitesOnlyQuery);
    }

    private File runCombineGVCFs(final List<String> inputs, final List<SimpleInterval> intervals, final String reference, final String[] extraArgs) {
        final File output = createTempFile("genotypegvcf", ".vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addReference(new File(reference))
                .addOutput(output);
        for (String input: inputs) {
            args.addArgument("V", input);
        }
        intervals.forEach(interval -> args.addArgument("L", interval.toString()));
        Arrays.stream(extraArgs).forEach(args::add);

        Utils.resetRandomGenerator();
        new Main().instanceMain(makeCommandLineArgs(args.getArgsList(), "CombineGVCFs"));
        return output;
    }

    private void testGenomicsDBAgainstCombineGVCFs(final List<String> vcfInputs, final List<SimpleInterval> intervals,
                                                   final String referenceFile, final String[] CombineGVCFArgs) throws IOException {
        final String workspace = createTempDir("genomicsdb-tests-").getAbsolutePath() + "/workspace";

        writeToGenomicsDB(vcfInputs, intervals, workspace, 0, false, 0, 1);
        checkJSONFilesAreWritten(workspace);
        for(SimpleInterval currInterval : intervals) {
            List<SimpleInterval> tmpList = new ArrayList<SimpleInterval>(Arrays.asList(currInterval));
            File expectedCombinedVCF = runCombineGVCFs(vcfInputs, tmpList, referenceFile, CombineGVCFArgs);
            checkGenomicsDBAgainstExpected(workspace, tmpList, expectedCombinedVCF.getAbsolutePath(), referenceFile, true, ATTRIBUTES_TO_IGNORE);
        }
    }

    @Test(groups = {"bucket"})
    public void testGenomicsDBImportGCSInputs() throws IOException {
        testGenomicsDBImporter(resolveLargeFilesAsCloudURIs(LOCAL_GVCFS), INTERVAL, COMBINED, b38_reference_20_21, true, 1);
    }

    @Test
    public void testGenomicsDBAbsolutePathDepndency() throws IOException {
        final File workspace = createTempDir("genomicsdb-tests-");
        final File workspace2 = createTempDir("genomicsdb-secondary-tests-");

        writeToGenomicsDB(LOCAL_GVCFS, INTERVAL, workspace.getAbsolutePath() + "/workspace", 0, false, 0, 1);
        checkJSONFilesAreWritten(workspace.getAbsolutePath() + "/workspace");
        Files.move(workspace.toPath(), workspace2.toPath(), StandardCopyOption.REPLACE_EXISTING);
        checkGenomicsDBAgainstExpected(workspace2.getAbsolutePath() + "/workspace", INTERVAL, COMBINED, b38_reference_20_21, true, ATTRIBUTES_TO_IGNORE);
    }

    @Test (enabled = true)
    public void testGenomicsDBAlleleSpecificAnnotations() throws IOException {
        testGenomicsDBAgainstCombineGVCFs(Arrays.asList(COMBINEGVCFS_TEST_DIR+"NA12878.AS.chr20snippet.g.vcf", COMBINEGVCFS_TEST_DIR+"NA12892.AS.chr20snippet.g.vcf"),
                new ArrayList<SimpleInterval>(Arrays.asList(new SimpleInterval("20", 10433000, 10700000))),
                b37_reference_20_21,
                new String[]{"-G", "StandardAnnotation", "-G", "AS_StandardAnnotation"});
    }

    /**
     * Converts a list of large file paths into equivalent cloud paths
     * This must be done non-statically because any failure during static initialization results in hard to understand
     * TestNG errors and it is possible for {@link BaseTest#getGCPTestInputPath()} to fail if the environment isn't
     * fully set up.
     *
     * The cloud bucket must be organized the same way as the local test files in order to resolve correctly.
     */
    private static List<String> resolveLargeFilesAsCloudURIs(final List<String> filenames){
        return filenames.stream()
                .map( filename -> filename.replace(publicTestDir, getGCPTestInputPath()))
                .peek( filename -> Assert.assertTrue(BucketUtils.isCloudStorageUrl(filename)))
                .collect(Collectors.toList());
    }

    @Test(dataProvider = "batchSizes")
    public void testGenomicsDBImportFileInputsInBatches(final int batchSize) throws IOException {
        testGenomicsDBImporterWithBatchSize(LOCAL_GVCFS, INTERVAL, COMBINED, batchSize);
    }

    @Test(dataProvider = "batchSizes")
    public void testGenomicsDBImportFileInputsInBatchesWithMultipleIntervals(final int batchSize) throws IOException {
        testGenomicsDBImporterWithBatchSize(LOCAL_GVCFS, MULTIPLE_INTERVALS, COMBINED_MULTI_INTERVAL, batchSize);
    }

    @Test(groups = {"bucket"}, dataProvider = "batchSizes")
    public void testGenomicsDBImportGCSInputsInBatches(final int batchSize) throws IOException {
        testGenomicsDBImporterWithBatchSize(resolveLargeFilesAsCloudURIs(LOCAL_GVCFS), INTERVAL, COMBINED, batchSize);
    }

    @DataProvider
    public Object[][] getThreads(){
        return new Object[][] {
                {1}, {2}, {5}
        };
    }

    @Test(groups = {"bucket"}, dataProvider = "getThreads")
    public void testDifferentThreadValuesFromABucket(final int threads) throws IOException {
        final List<String> vcfInputs = resolveLargeFilesAsCloudURIs(LOCAL_GVCFS);
        final String workspace = createTempDir("genomicsdb-tests-").getAbsolutePath() + "/workspace";

        writeToGenomicsDB(vcfInputs, INTERVAL, workspace, 0, false, 0, threads);
        checkJSONFilesAreWritten(workspace);
        checkGenomicsDBAgainstExpected(workspace, INTERVAL, COMBINED, b38_reference_20_21, true, ATTRIBUTES_TO_IGNORE);
    }

    @Test(dataProvider = "getThreads")
    public void testDifferentThreadValuesLocally(final int threads) throws IOException {
        final List<String> vcfInputs = LOCAL_GVCFS;
        final String workspace = createTempDir("genomicsdb-tests-").getAbsolutePath() + "/workspace";

        writeToGenomicsDB(vcfInputs, INTERVAL, workspace, 0, false, 0, threads);
        checkJSONFilesAreWritten(workspace);
        checkGenomicsDBAgainstExpected(workspace, INTERVAL, COMBINED, b38_reference_20_21, true, ATTRIBUTES_TO_IGNORE);
    }
    /**
     *
     * @throws CommandLineException.OutOfRangeArgumentValue  Value must be >= 1024 bytes
     */
    @Test(expectedExceptions = CommandLineException.OutOfRangeArgumentValue.class)
    public void testZeroVCFBufferSize() throws IOException {
        testGenomicsDBImportWithZeroBufferSize(LOCAL_GVCFS, INTERVAL, COMBINED);
    }


    private void testGenomicsDBImporter(final List<String> vcfInputs, final List<SimpleInterval> intervals,
                                        final String expectedCombinedVCF, final String referenceFile,
                                        final boolean testAll, final int threads) throws IOException {
        final String workspace = createTempDir("genomicsdb-tests-").getAbsolutePath() + "/workspace";
        writeToGenomicsDB(vcfInputs, intervals, workspace, 0, false, 0, 1);

        checkGenomicsDBAgainstExpected(workspace, intervals, expectedCombinedVCF, referenceFile, testAll, ATTRIBUTES_TO_IGNORE);
    }

    private void testGenomicsDBImporter_newMQ(final List<String> vcfInputs, final List<SimpleInterval> intervals,
                                        final String expectedCombinedVCF, final String referenceFile,
                                        final boolean testAll, final List<String> attributesToIgnore) throws IOException {
        final String workspace = createTempDir("genomicsdb-tests-").getAbsolutePath() + "/workspace";

        writeToGenomicsDB(vcfInputs, intervals, workspace, 0, false, 0, 1);
        checkJSONFilesAreWritten(workspace);

        checkGenomicsDBAgainstExpected(workspace, intervals, expectedCombinedVCF, referenceFile, testAll, attributesToIgnore);
    }

    private void testGenomicsDBImporterWithBatchSize(final List<String> vcfInputs, final List<SimpleInterval> intervals,
                                                     final String expectedCombinedVCF, final int batchSize) throws IOException {
        final String workspace = createTempDir("genomicsdb-batchsize-tests-").getAbsolutePath() + "/workspace-" + batchSize;

        writeToGenomicsDB(vcfInputs, intervals, workspace, batchSize, false, 0, 1);
        checkJSONFilesAreWritten(workspace);
        checkGenomicsDBAgainstExpected(workspace, intervals, expectedCombinedVCF, b38_reference_20_21, true, ATTRIBUTES_TO_IGNORE);
    }

    private void testGenomicsDBImportWithZeroBufferSize(final List<String> vcfInputs, final List<SimpleInterval> intervals,
                                                        final String expectedCombinedVCF) throws IOException {
        final String workspace = createTempDir("genomicsdb-buffersize-tests-").getAbsolutePath() + "/workspace";

        writeToGenomicsDB(vcfInputs, intervals, workspace, 0, true, 0, 1);
        checkJSONFilesAreWritten(workspace);
        checkGenomicsDBAgainstExpected(workspace, intervals, expectedCombinedVCF, b38_reference_20_21, true, ATTRIBUTES_TO_IGNORE);

    }

    private void writeToGenomicsDB(final List<String> vcfInputs, final List<SimpleInterval> intervals, final String workspace,
                                   final int batchSize, final Boolean useBufferSize, final int bufferSizePerSample, int threads) {
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addArgument(GenomicsDBImport.WORKSPACE_ARG_LONG_NAME, workspace);
        intervals.forEach(interval -> args.addArgument("L", IntervalUtils.locatableToString(interval)));
        vcfInputs.forEach(vcf -> args.addArgument("V", vcf));
        args.addArgument("batch-size", String.valueOf(batchSize));
        args.addArgument(GenomicsDBImport.VCF_INITIALIZER_THREADS_LONG_NAME, String.valueOf(threads));
        if (useBufferSize)
            args.addArgument("genomicsdb-vcf-buffer-size", String.valueOf(bufferSizePerSample));

        runCommandLine(args);
    }

    private static void checkJSONFilesAreWritten(final String workspace) {
        Assert.assertTrue(BucketUtils.fileExists(IOUtils.appendPathToDir(workspace, GenomicsDBConstants.DEFAULT_VIDMAP_FILE_NAME)));
        Assert.assertTrue(BucketUtils.fileExists(IOUtils.appendPathToDir(workspace, GenomicsDBConstants.DEFAULT_CALLSETMAP_FILE_NAME)));
        Assert.assertTrue(BucketUtils.fileExists(IOUtils.appendPathToDir(workspace, GenomicsDBConstants.DEFAULT_VCFHEADER_FILE_NAME)));
    }

    private static void checkGenomicsDBAgainstExpected(final String workspace, final List<SimpleInterval> intervals,
                                                       final String expectedCombinedVCF, final String referenceFile,
                                                       final boolean testAll, final List<String> attributesToIgnore) throws IOException {
        checkGenomicsDBAgainstExpected(workspace, intervals,
                expectedCombinedVCF, referenceFile,
                testAll,
                attributesToIgnore,
                false,
                false);
    }

    private static void checkGenomicsDBAgainstExpected(final String workspace, final List<SimpleInterval> intervals,
                                                       final String expectedCombinedVCF, final String referenceFile,
                                                       final boolean testAll,
                                                       final List<String> attributesToIgnore,
                                                       final boolean produceGTField,
                                                       final boolean sitesOnlyQuery) throws IOException {
        final GenomicsDBFeatureReader<VariantContext, PositionalBufferedStream> genomicsDBFeatureReader =
                getGenomicsDBFeatureReader(workspace, referenceFile, produceGTField, sitesOnlyQuery);

        final AbstractFeatureReader<VariantContext, LineIterator> combinedVCFReader =
                AbstractFeatureReader.getFeatureReader(expectedCombinedVCF, new VCFCodec(), true);


        intervals.forEach(interval -> {
            try (CloseableTribbleIterator<VariantContext> actualVcs =
                         genomicsDBFeatureReader.query(interval.getContig(), interval.getStart(), interval.getEnd());

                 CloseableTribbleIterator<VariantContext> expectedVcs =
                         combinedVCFReader.query(interval.getContig(), interval.getStart(), interval.getEnd())) {

                BaseTest.assertCondition(actualVcs, expectedVcs, (a, e) -> {
                        // Test that the VCs match
                    if (testAll) {
                        // To correct a discrepancy between genotypeGVCFs which outputs empty genotypes as "./." and GenomicsDB
                        // which returns them as "." we simply remap the empty ones to be consistent for comparison
                        List<Genotype> genotypes = a.getGenotypes().stream()
                                .map(g -> g.getGenotypeString().equals(".")?new GenotypeBuilder(g).alleles(GATKVariantContextUtils.noCallAlleles(2)).make():g)
                                .collect(Collectors.toList());
                        a = new VariantContextBuilder(a).genotypes(genotypes).make();
                        VariantContextTestUtils.assertVariantContextsAreEqualAlleleOrderIndependent(a, e, attributesToIgnore, VCF_HEADER);

                        // Test only that the genotypes match
                    } else {
                        List<Genotype> genotypes = e.getGenotypes().stream()
                                .map(g -> g.getGenotypeString().equals(".")?new GenotypeBuilder(g).alleles(Collections.emptyList()).make():g)
                                .collect(Collectors.toList());
                        e = new VariantContextBuilder(e).genotypes(genotypes).make();
                        VariantContextTestUtils.assertVariantContextsHaveSameGenotypes(a, e);
                    }
                });
            } catch (IOException e) {
                Assert.fail(e.getMessage(), e);
            }
        });
    }

    @DataProvider
    public Iterator<Object[]> getOrderingTests(){
        final File outOfOrderSampleMap = getSampleMapFile(
                        "HG00268\t" + HG_00268 + "\n" +
                        "NA19625\t" + NA_19625 + "\n" +
                        "HG00096\t" + HG_00096);

        final List<Integer> batchSizes = Arrays.asList(0, 1, 2, 3, 4);
        final List<Object[]> results = new ArrayList<>();
        for( final Integer batchSize: batchSizes){
            // -V in order
            results.add(new Object[] {new ArgumentsBuilder()
                    .addArgument(GenomicsDBImport.BATCHSIZE_ARG_LONG_NAME, String.valueOf(batchSize))
                    .addVCF(new File(HG_00096))
                    .addVCF(new File(HG_00268))
                    .addVCF(new File(NA_19625))});

            // -V out of order
            results.add(new Object[] {new ArgumentsBuilder()
                    .addArgument(GenomicsDBImport.BATCHSIZE_ARG_LONG_NAME, String.valueOf(batchSize))
                    .addVCF(new File(HG_00268))
                    .addVCF(new File(NA_19625))
                    .addVCF(new File(HG_00096))});

            //in order sample map
            results.add(new Object[] {new ArgumentsBuilder()
                    .addArgument(GenomicsDBImport.BATCHSIZE_ARG_LONG_NAME, String.valueOf(batchSize))
                    .addFileArgument(GenomicsDBImport.SAMPLE_NAME_MAP_LONG_NAME, createInOrderSampleMap())});

            //out of order sample map
            results.add(new Object[] {new ArgumentsBuilder()
                    .addArgument(GenomicsDBImport.BATCHSIZE_ARG_LONG_NAME, String.valueOf(batchSize))
                    .addFileArgument(GenomicsDBImport.SAMPLE_NAME_MAP_LONG_NAME, outOfOrderSampleMap)});

            //out of order sample map with multiple threads
            results.add(new Object[] {new ArgumentsBuilder()
                    .addArgument(GenomicsDBImport.BATCHSIZE_ARG_LONG_NAME, String.valueOf(batchSize))
                    .addFileArgument(GenomicsDBImport.SAMPLE_NAME_MAP_LONG_NAME, outOfOrderSampleMap)
                    .addArgument(GenomicsDBImport.VCF_INITIALIZER_THREADS_LONG_NAME, "2")});
        }
        return results.iterator();
    }

    @Test
    public void testSampleNameWithSpaces() throws IOException {
        final File outOfOrderSampleMap = getSampleMapFile(
                "HG00268 withSpaces\t" + HG_00268_WITH_SPACES + "\n" +
                        "NA19625\t" + NA_19625 + "\n" +
                        "HG00096\t" + HG_00096 );

        final String workspace = createTempDir("gendbtest").getAbsolutePath() + "/workspace";

        ArgumentsBuilder args = new ArgumentsBuilder()
                .addArgument(GenomicsDBImport.BATCHSIZE_ARG_LONG_NAME, String.valueOf(2))
                .addFileArgument(GenomicsDBImport.SAMPLE_NAME_MAP_LONG_NAME, outOfOrderSampleMap)
                .addArgument("L", IntervalUtils.locatableToString(SMALLER_INTERVAL.get(0)))
                .addArgument(GenomicsDBImport.WORKSPACE_ARG_LONG_NAME, workspace);

        runCommandLine(args);
        checkJSONFilesAreWritten(workspace);
        checkGenomicsDBAgainstExpected(workspace, SMALLER_INTERVAL, COMBINED_WITHSPACES, b38_reference_20_21, true, ATTRIBUTES_TO_IGNORE);
    }

    @Test(dataProvider = "getOrderingTests")
    public void testSampleNameOrdering(final ArgumentsBuilder args) throws IOException {
        final String workspace = createTempDir("gendbtest").getAbsolutePath() + "/workspace";

        args.addArgument("L", IntervalUtils.locatableToString(INTERVAL.get(0)))
            .addArgument(GenomicsDBImport.WORKSPACE_ARG_LONG_NAME, workspace);

        runCommandLine(args);
        checkJSONFilesAreWritten(workspace);
        checkGenomicsDBAgainstExpected(workspace, INTERVAL, COMBINED, b38_reference_20_21, true, ATTRIBUTES_TO_IGNORE);
    }

    private static File createInOrderSampleMap() {
        final String sampleFileContents =
                "HG00096\t" +HG_00096 +"\n" +
                "HG00268\t"+ HG_00268 + "\n" +
                "NA19625\t"+ NA_19625;

        return getSampleMapFile(sampleFileContents);
    }

    private static File getSampleMapFile(final String sampleFileContents) {
        final File sampleNameMap = IOUtils.writeTempFile(sampleFileContents, "sampleNameMap", ".txt");
        sampleNameMap.deleteOnExit();
        return sampleNameMap;
    }

    private static File getSampleMapFile(final Map<String, String> mapping){
        return getSampleMapFile(mapping.entrySet()
                .stream()
                .map( pair -> pair.getKey() + "\t" + pair.getValue())
                .collect(Collectors.joining("\n")));
    }

    @DataProvider
    public static Iterator<Object[]> getRenameCombinations() {
        final Map<String,String> noRemapping = new LinkedHashMap<>();
        noRemapping.put("s1", "s1");
        noRemapping.put("s2", "s2");
        noRemapping.put("s3", "s3");

        final Map<String,String> sameInput = new LinkedHashMap<>();
        sameInput.put("s1", "s1");
        sameInput.put("s2", "s1");
        sameInput.put("s3", "s1");


        final Map<String,String> sameInputWeirdOrder = new LinkedHashMap<>();
        sameInputWeirdOrder.put("s3", "s1");
        sameInputWeirdOrder.put("s1", "s1");
        sameInputWeirdOrder.put("s2", "s1");

        final Map<String,String> swizzled = new LinkedHashMap<>();
        swizzled.put("s2","s1");
        swizzled.put("s3","s2");
        swizzled.put("s1","s3");

        final Map<String,String> multipleOutOfOrderRenamingsAcrossBatches = new LinkedHashMap<>();
        multipleOutOfOrderRenamingsAcrossBatches.put("s1", "s1");
        multipleOutOfOrderRenamingsAcrossBatches.put("s2", "s2");
        multipleOutOfOrderRenamingsAcrossBatches.put("s1_Renamed", "s1");
        multipleOutOfOrderRenamingsAcrossBatches.put("Renamed_s2", "s2");
        multipleOutOfOrderRenamingsAcrossBatches.put("s4", "s3");
        multipleOutOfOrderRenamingsAcrossBatches.put("s3", "s3");
        multipleOutOfOrderRenamingsAcrossBatches.put("someOtherSample", "s4");


        final List<Integer> batchSizes = Arrays.asList(0, 1, 4);
        final List<Integer> threads = Arrays.asList(1, 2);
        final List<Map<String, String>> mappings = Arrays.asList(noRemapping, sameInput, sameInputWeirdOrder, swizzled, multipleOutOfOrderRenamingsAcrossBatches);
        final List<Object[]> out = new ArrayList<>();
        for(final Map<String,String> mapping : mappings){
            for(final int batchSize :batchSizes){
                for(final int threading : threads){
                    out.add( new Object[]{mapping, threading, batchSize});
                }
            }
        }
        return out.iterator();
    }


    @Test(dataProvider = "getRenameCombinations")
    public void testRenamingSamples(final Map<String, String> renamingMap, final int threads, final int batchSize) throws IOException {
        final LinkedHashMap<String, String> sampleMap = new LinkedHashMap<>(renamingMap);
        sampleMap.replaceAll( (newSampleName, originalSampleName)-> createInputVCF(originalSampleName).getAbsolutePath());

        final File sampleMapFile = getSampleMapFile(sampleMap);

        final String workspace = createTempDir("workspace").getAbsolutePath();
        Files.delete(Paths.get(workspace));
        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addArgument(GenomicsDBImport.SAMPLE_NAME_MAP_LONG_NAME, sampleMapFile.getAbsolutePath())
                .addArgument(GenomicsDBImport.WORKSPACE_ARG_LONG_NAME, new File(workspace).getAbsolutePath())
                .addArgument(GenomicsDBImport.VCF_INITIALIZER_THREADS_LONG_NAME, String.valueOf(threads))
                .addArgument(GenomicsDBImport.BATCHSIZE_ARG_LONG_NAME, String.valueOf(batchSize))
                .addArgument("L", IntervalUtils.locatableToString(INTERVAL.get(0)));

        runCommandLine(args);
        final Set<String> expectedSampleNames = sampleMap.keySet();
        try(final GenomicsDBFeatureReader<VariantContext, PositionalBufferedStream> reader = getGenomicsDBFeatureReader(workspace, b37_reference_20_21)) {
            final CloseableTribbleIterator<VariantContext> iterator = reader.iterator();
            Assert.assertTrue(iterator.hasNext(), "expected to see a variant");
            Assert.assertTrue(expectedSampleNames.size() > 0);
            Assert.assertEquals(expectedSampleNames.size(), renamingMap.size());
            iterator.forEachRemaining(vc -> {
                   Assert.assertEquals(vc.getSampleNames(), expectedSampleNames);
                expectedSampleNames.forEach( sample -> {
                   Assert.assertEquals(vc.getGenotype(sample).getAnyAttribute(SAMPLE_NAME_KEY), renamingMap.get(sample));
                   Assert.assertEquals(vc.getGenotype(sample).getAnyAttribute(ANOTHER_ATTRIBUTE_KEY), 10); //check another attribute just to make sure we're not mangling things
                });
            });
        }

    }

    private static File createInputVCF(final String sampleName) {
        final String contig = "chr20";
        final SAMSequenceDictionary dict = new SAMSequenceDictionary(
                Collections.singletonList(new SAMSequenceRecord(contig, 64444167)));

        final VCFFormatHeaderLine formatField = new VCFFormatHeaderLine(SAMPLE_NAME_KEY, 1, VCFHeaderLineType.String,
                                                                        "the name of the sample this genotype came from");
        final Set<VCFHeaderLine> headerLines = new HashSet<>();
        headerLines.add(formatField);
        headerLines.add(new VCFFormatHeaderLine(ANOTHER_ATTRIBUTE_KEY, 1, VCFHeaderLineType.Integer, "Another value"));
        headerLines.add(VCFStandardHeaderLines.getFormatLine("GT"));

        final File out = createTempFile(sampleName +"_", ".vcf");
        try (final VariantContextWriter writer = GATKVariantContextUtils.createVCFWriter(out, dict, false,
                                                                                         Options.INDEX_ON_THE_FLY)) {
            final VCFHeader vcfHeader = new VCFHeader(headerLines, Collections.singleton(sampleName));
            vcfHeader.setSequenceDictionary(dict);
            writer.writeHeader(vcfHeader);
            final Allele Aref = Allele.create("A", true);
            final Allele C = Allele.create("C");
            final List<Allele> alleles = Arrays.asList(Aref, C);
            final VariantContext variant = new VariantContextBuilder("invented", contig, INTERVAL.get(0).getStart(), INTERVAL.get(0).getStart(), alleles)
                    .genotypes(new GenotypeBuilder(sampleName, alleles).attribute(SAMPLE_NAME_KEY, sampleName)
                                       .attribute(ANOTHER_ATTRIBUTE_KEY, 10).make())
                    .make();
            writer.add(variant);
            return out;
        }
    }

    @Test(expectedExceptions = CommandLineException.class)
    public void testCantSpecifyVCFAndSampleNameFile(){
        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addArgument(GenomicsDBImport.SAMPLE_NAME_MAP_LONG_NAME, createInOrderSampleMap().getAbsolutePath())
                .addArgument(StandardArgumentDefinitions.VARIANT_LONG_NAME, HG_00096)
                .addArgument(GenomicsDBImport.WORKSPACE_ARG_LONG_NAME, createTempDir("workspace").getAbsolutePath())
                .addArgument("L",  IntervalUtils.locatableToString(INTERVAL.get(0)));
        runCommandLine(args);
    }

    @Test(expectedExceptions = CommandLineException.MissingArgument.class)
    public void testRequireOneOfVCFOrSampleNameFile(){
        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addArgument(GenomicsDBImport.WORKSPACE_ARG_LONG_NAME, createTempDir("workspace").getAbsolutePath())
                .addArgument("L", "1:1-10");

        runCommandLine(args);
    }

    @Test
    public void testGenomicsDBImportWithoutDBField() throws IOException {
        //Test for https://github.com/broadinstitute/gatk/issues/3736
        final List<String> vcfInputs = Arrays.asList(NA_24385);
        final String workspace = createTempDir("genomicsdb-tests").getAbsolutePath() + "/workspace";
	writeToGenomicsDB(vcfInputs, INTERVAL_3736, workspace, 0, false, 0, 1);
    }

    @Test
    public void testLongWorkspacePath() throws IOException {
        //Test for https://github.com/broadinstitute/gatk/issues/4160
        final List<String> vcfInputs = LOCAL_GVCFS;
        final String workspace = createTempDir("long_aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa_genomicsdb").getAbsolutePath() + "/should_not_fail_aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa";
        writeToGenomicsDB(vcfInputs, INTERVAL, workspace, 0, false, 0, 1);
    }

    @Test
    public void testCommandIncludedInOutputHeader() throws IOException {
        final List<String> vcfInputs = LOCAL_GVCFS;
        final String workspace = createTempDir("genomicsdb-tests").getAbsolutePath() + "/workspace";

        writeToGenomicsDB(vcfInputs, INTERVAL, workspace, 0, false, 0, 1);
        try(final GenomicsDBFeatureReader<VariantContext, PositionalBufferedStream> genomicsDBFeatureReader =
                    getGenomicsDBFeatureReader(workspace, b38_reference_20_21))
        {
            final VCFHeader header = (VCFHeader) genomicsDBFeatureReader.getHeader();
            final Optional<VCFHeaderLine> commandLineHeaderLine = header.getMetaDataInSortedOrder().stream()
                    .filter(line -> line.getValue().contains(GenomicsDBImport.class.getSimpleName()))
                    .findAny();

            Assert.assertTrue(commandLineHeaderLine.isPresent(), "no headerline was present containing information about the GenomicsDBImport command");
        }


    }

    @Test
    public void testPreserveContigOrderingInHeader() throws IOException {
        final String workspace = createTempDir("testPreserveContigOrderingInHeader-").getAbsolutePath() + "/workspace";
        ArrayList<SimpleInterval> intervals = new ArrayList<SimpleInterval>(Arrays.asList(new SimpleInterval("chr20", 17959479, 17959479)));
        writeToGenomicsDB(Arrays.asList(GENOMICSDB_TEST_DIR + "testHeaderContigLineSorting1.g.vcf",
                GENOMICSDB_TEST_DIR + "testHeaderContigLineSorting2.g.vcf"), intervals, workspace, 0, false, 0, 1);

        try ( final GenomicsDBFeatureReader<VariantContext, PositionalBufferedStream> genomicsDBFeatureReader =
                      getGenomicsDBFeatureReader(workspace, b38_reference_20_21);

             final AbstractFeatureReader<VariantContext, LineIterator> inputGVCFReader =
                      AbstractFeatureReader.getFeatureReader(GENOMICSDB_TEST_DIR + "testHeaderContigLineSorting1.g.vcf", new VCFCodec(), true);
        ) {
            final SAMSequenceDictionary dictionaryFromGenomicsDB = ((VCFHeader)genomicsDBFeatureReader.getHeader()).getSequenceDictionary();
            final SAMSequenceDictionary dictionaryFromInputGVCF =  ((VCFHeader)inputGVCFReader.getHeader()).getSequenceDictionary();

            Assert.assertEquals(dictionaryFromGenomicsDB, dictionaryFromInputGVCF, "Sequence dictionary from GenomicsDB does not match original sequence dictionary from input GVCF");
        }

    }
    private static GenomicsDBFeatureReader<VariantContext, PositionalBufferedStream> getGenomicsDBFeatureReader(
            final String workspace, final String reference,
            final boolean produceGTField) throws IOException {
        return getGenomicsDBFeatureReader(workspace, reference,
                produceGTField, false);
    }

    private static GenomicsDBFeatureReader<VariantContext, PositionalBufferedStream> getGenomicsDBFeatureReader(
            final String workspace, final String reference,
            final boolean produceGTField,
            final boolean sitesOnlyQuery) throws IOException {
       String workspaceAbsPath = BucketUtils.makeFilePathAbsolute(workspace);GenomicsDBExportConfiguration.ExportConfiguration .Builder exportConfigurationBuilder = GenomicsDBExportConfiguration.ExportConfiguration.newBuilder()
                .setWorkspace(workspace)
                .setReferenceGenome(reference)
                .setVidMappingFile(IOUtils.appendPathToDir(workspaceAbsPath, GenomicsDBConstants.DEFAULT_VIDMAP_FILE_NAME))
                .setCallsetMappingFile(IOUtils.appendPathToDir(workspaceAbsPath, GenomicsDBConstants.DEFAULT_CALLSETMAP_FILE_NAME))
                .setVcfHeaderFilename(IOUtils.appendPathToDir(workspaceAbsPath, GenomicsDBConstants.DEFAULT_VCFHEADER_FILE_NAME))
                .setProduceGTField(produceGTField)
                .setSitesOnlyQuery(sitesOnlyQuery)
               .setGenerateArrayNameFromPartitionBounds(true);
        GenomicsDBVidMapProto.VidMappingPB vidMapPB = null;
        try {
            vidMapPB = org.broadinstitute.hellbender.tools.genomicsdb.GenomicsDBUtils.getProtobufVidMappingFromJsonFile(IOUtils.appendPathToDir(workspace, GenomicsDBConstants.DEFAULT_VIDMAP_FILE_NAME));
        }
        catch (final IOException e) {
            throw new UserException("Could not open vid json file "+GenomicsDBConstants.DEFAULT_VIDMAP_FILE_NAME, e);
        }
        HashMap<String, Integer> fieldNameToIndexInVidFieldsList =
                org.broadinstitute.hellbender.tools.genomicsdb.GenomicsDBUtils.getFieldNameToListIndexInProtobufVidMappingObject(vidMapPB);

        vidMapPB = org.broadinstitute.hellbender.tools.genomicsdb.GenomicsDBUtils.updateINFOFieldCombineOperation(vidMapPB, fieldNameToIndexInVidFieldsList,
                GATKVCFConstants.RAW_MAPPING_QUALITY_WITH_DEPTH_KEY, "element_wise_sum");

        if(vidMapPB != null) {
            exportConfigurationBuilder.setVidMapping(vidMapPB);
        }

        return new GenomicsDBFeatureReader<>(exportConfigurationBuilder.build(), new BCF2Codec(), Optional.empty());
    }

    private static GenomicsDBFeatureReader<VariantContext, PositionalBufferedStream> getGenomicsDBFeatureReader(
            final String workspace, final String reference) throws IOException {
        return getGenomicsDBFeatureReader(workspace, reference, false);
    }

    @Test(expectedExceptions = GenomicsDBImport.UnableToCreateGenomicsDBWorkspace.class)
    public void testYouCantWriteIntoAnExistingDirectory(){
        // this actually creates the directory on disk, not just the file name.
        final String workspace = createTempDir("workspace").getAbsolutePath();
        writeToGenomicsDB(LOCAL_GVCFS, INTERVAL, workspace, 0, false, 0, 1);
    }

    /*@Test(groups = {"bucket"})
    public void testWriteToAndQueryFromGCS() throws IOException {
        final String workspace = BucketUtils.randomRemotePath(getGCPTestStaging(), "", "") + "/";
        BucketUtils.deleteOnExit(workspace);
        writeToGenomicsDB(LOCAL_GVCFS, INTERVAL, workspace, 0, false, 0, 1);
        checkJSONFilesAreWritten(workspace);
        checkGenomicsDBAgainstExpected(workspace, INTERVAL, COMBINED, b38_reference_20_21, true);
    }

    @Test(groups = {"bucket"}, expectedExceptions = GenomicsDBImport.UnableToCreateGenomicsDBWorkspace.class)
    public void testWriteToExistingGCSDirectory() throws IOException {
        final String workspace = BucketUtils.randomRemotePath(getGCPTestStaging(), "", "") + "/";
        BucketUtils.deleteOnExit(workspace);
        int rc = GenomicsDBUtils.createTileDBWorkspace(workspace, false);
        Assert.assertEquals(rc, 0);
        writeToGenomicsDB(LOCAL_GVCFS, INTERVAL, workspace, 0, false, 0, 1);
    }*/
}
