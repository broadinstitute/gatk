package org.broadinstitute.hellbender.tools.genomicsdb;

import com.intel.genomicsdb.GenomicsDBFeatureReader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.PositionalBufferedStream;
import htsjdk.variant.bcf2.BCF2Codec;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.test.VariantContextTestUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

public final class GenomicsDBImportIntegrationTest extends CommandLineProgramTest {

    private static final String HG_00096 = largeFileTestDir + "gvcfs/HG00096.g.vcf.gz";
    private static final String HG_00268 = largeFileTestDir + "gvcfs/HG00268.g.vcf.gz";
    private static final String NA_19625 = largeFileTestDir + "gvcfs/NA19625.g.vcf.gz";
    private static final List<String> LOCAL_GVCFS = Arrays.asList(HG_00096, HG_00268, NA_19625);
    private static final String GENOMICSDB_TEST_DIR = publicTestDir + "org/broadinstitute/hellbender/tools/GenomicsDBImport/";

    private static final String COMBINED = largeFileTestDir + "gvcfs/combined.gatk3.7.g.vcf.gz";
    private static final SimpleInterval INTERVAL = new SimpleInterval("chr20", 17960187, 17981445);
    private static final VCFHeader VCF_HEADER = VariantContextTestUtils.getCompleteHeader();


    private static final String SAMPLE_NAME_KEY = "SN";
    private static final String ANOTHER_ATTRIBUTE_KEY = "AA";


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
        testGenomicsDBImporter(LOCAL_GVCFS, INTERVAL, COMBINED);
    }

    @Test(groups = {"bucket"})
    public void testGenomicsDBImportGCSInputs() throws IOException {
        testGenomicsDBImporter(resolveLargeFilesAsCloudURIs(LOCAL_GVCFS), INTERVAL, COMBINED);
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
        checkGenomicsDBAgainstExpected(workspace, INTERVAL, COMBINED);
    }

    @Test(dataProvider = "getThreads")
    public void testDifferentThreadValuesLocally(final int threads) throws IOException {
        final List<String> vcfInputs = LOCAL_GVCFS;
        final String workspace = createTempDir("genomicsdb-tests-").getAbsolutePath() + "/workspace";

        writeToGenomicsDB(vcfInputs, INTERVAL, workspace, 0, false, 0, threads);
        checkJSONFilesAreWritten(workspace);
        checkGenomicsDBAgainstExpected(workspace, INTERVAL, COMBINED);
    }
    /**
     *
     * @throws CommandLineException.OutOfRangeArgumentValue  Value must be >= 1024 bytes
     */
    @Test(expectedExceptions = CommandLineException.OutOfRangeArgumentValue.class)
    public void testZeroVCFBufferSize() throws IOException {
        testGenomicsDBImportWithZeroBufferSize(LOCAL_GVCFS, INTERVAL, COMBINED);
    }


    private void testGenomicsDBImporter(final List<String> vcfInputs, final SimpleInterval interval, final String expectedCombinedVCF) throws IOException {
        final String workspace = createTempDir("genomicsdb-tests-").getAbsolutePath() + "/workspace";

        writeToGenomicsDB(vcfInputs, interval, workspace, 0, false, 0, 1);
        checkJSONFilesAreWritten(workspace);
        checkGenomicsDBAgainstExpected(workspace, interval, expectedCombinedVCF);
    }

    private void testGenomicsDBImporterWithBatchSize(final List<String> vcfInputs, final SimpleInterval interval, final String expectedCombinedVCF, final int batchSize) throws IOException {
        final String workspace = createTempDir("genomicsdb-batchsize-tests-").getAbsolutePath() + "/workspace-" + batchSize;

        writeToGenomicsDB(vcfInputs, interval, workspace, batchSize, false, 0, 1);
        checkJSONFilesAreWritten(workspace);
        checkGenomicsDBAgainstExpected(workspace, interval, expectedCombinedVCF);
    }

    private void testGenomicsDBImportWithZeroBufferSize(final List<String> vcfInputs, final SimpleInterval interval, final String expectedCombinedVCF) throws IOException {
        final String workspace = createTempDir("genomicsdb-buffersize-tests-").getAbsolutePath() + "/workspace";

        writeToGenomicsDB(vcfInputs, interval, workspace, 0, true, 0, 1);
        checkJSONFilesAreWritten(workspace);
        checkGenomicsDBAgainstExpected(workspace, interval, expectedCombinedVCF);

    }

    private void writeToGenomicsDB(final List<String> vcfInputs, final SimpleInterval interval, final String workspace,
                                   final int batchSize, final Boolean useBufferSize, final int bufferSizePerSample, int threads) {
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addArgument("genomicsDBWorkspace", workspace);
        args.addArgument("L", IntervalUtils.locatableToString(interval));
        vcfInputs.forEach(vcf -> args.addArgument("V", vcf));
        args.addArgument("batchSize", String.valueOf(batchSize));
        args.addArgument(GenomicsDBImport.VCF_INITIALIZER_THREADS_LONG_NAME, String.valueOf(threads));
        if (useBufferSize)
            args.addArgument("genomicsDBVCFBufferSize", String.valueOf(bufferSizePerSample));

        runCommandLine(args);
    }

    private static void checkJSONFilesAreWritten(final String workspace) {
        Assert.assertTrue(new File(workspace, GenomicsDBConstants.DEFAULT_VIDMAP_FILE_NAME).exists());
        Assert.assertTrue(new File(workspace, GenomicsDBConstants.DEFAULT_CALLSETMAP_FILE_NAME).exists());
	Assert.assertTrue(new File(workspace, GenomicsDBConstants.DEFAULT_VCFHEADER_FILE_NAME).exists());
    }

    private static void checkGenomicsDBAgainstExpected(final String workspace, final SimpleInterval interval, final String expectedCombinedVCF) throws IOException {
        final GenomicsDBFeatureReader<VariantContext, PositionalBufferedStream> genomicsDBFeatureReader =
                new GenomicsDBFeatureReader<>(
                        new File(workspace, GenomicsDBConstants.DEFAULT_VIDMAP_FILE_NAME).getAbsolutePath(),
                        new File(workspace, GenomicsDBConstants.DEFAULT_CALLSETMAP_FILE_NAME).getAbsolutePath(),
                        workspace,
                        GenomicsDBConstants.DEFAULT_ARRAY_NAME,
                        b38_reference_20_21,
			new File(workspace, GenomicsDBConstants.DEFAULT_VCFHEADER_FILE_NAME).getAbsolutePath(),
			new BCF2Codec());

        final AbstractFeatureReader<VariantContext, LineIterator> combinedVCFReader =
                AbstractFeatureReader.getFeatureReader(expectedCombinedVCF, new VCFCodec(), true);

        try (CloseableTribbleIterator<VariantContext> actualVcs =
                     genomicsDBFeatureReader.query(interval.getContig(), interval.getStart(), interval.getEnd());

             CloseableTribbleIterator<VariantContext> expectedVcs =
                     combinedVCFReader.query(interval.getContig(), interval.getStart(), interval.getEnd())) {

            BaseTest.assertCondition(actualVcs, expectedVcs, (a, e) -> {
                VariantContextTestUtils.assertVariantContextsAreEqualAlleleOrderIndependent(a, e, Collections.emptyList(), VCF_HEADER);
            });
        }
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
                    .addArgument(GenomicsDBImport.BATCHSIZE_ARG_NAME, String.valueOf(batchSize))
                    .addVCF(new File(HG_00096))
                    .addVCF(new File(HG_00268))
                    .addVCF(new File(NA_19625))});

            // -V out of order
            results.add(new Object[] {new ArgumentsBuilder()
                    .addArgument(GenomicsDBImport.BATCHSIZE_ARG_NAME, String.valueOf(batchSize))
                    .addVCF(new File(HG_00268))
                    .addVCF(new File(NA_19625))
                    .addVCF(new File(HG_00096))});

            //in order sample map
            results.add(new Object[] {new ArgumentsBuilder()
                    .addArgument(GenomicsDBImport.BATCHSIZE_ARG_NAME, String.valueOf(batchSize))
                    .addFileArgument(GenomicsDBImport.SAMPLE_NAME_MAP_LONG_NAME, createInOrderSampleMap())});

            //out of order sample map
            results.add(new Object[] {new ArgumentsBuilder()
                    .addArgument(GenomicsDBImport.BATCHSIZE_ARG_NAME, String.valueOf(batchSize))
                    .addFileArgument(GenomicsDBImport.SAMPLE_NAME_MAP_LONG_NAME, outOfOrderSampleMap)});

            //out of order sample map with multiple threads
            results.add(new Object[] {new ArgumentsBuilder()
                    .addArgument(GenomicsDBImport.BATCHSIZE_ARG_NAME, String.valueOf(batchSize))
                    .addFileArgument(GenomicsDBImport.SAMPLE_NAME_MAP_LONG_NAME, outOfOrderSampleMap)
                    .addArgument(GenomicsDBImport.VCF_INITIALIZER_THREADS_LONG_NAME, "2")});
        }
        return results.iterator();
    }

    @Test(dataProvider = "getOrderingTests")
    public void testSampleNameOrdering(final ArgumentsBuilder args) throws IOException {
        final String workspace = createTempDir("gendbtest").getAbsolutePath() + "/workspace";

        args.addArgument("L", IntervalUtils.locatableToString(INTERVAL))
            .addArgument(GenomicsDBImport.WORKSPACE_ARG_NAME, workspace);

        runCommandLine(args);
        checkJSONFilesAreWritten(workspace);
        checkGenomicsDBAgainstExpected(workspace, INTERVAL, COMBINED);
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
                .addArgument(GenomicsDBImport.WORKSPACE_ARG_NAME, new File(workspace).getAbsolutePath())
                .addArgument(GenomicsDBImport.VCF_INITIALIZER_THREADS_LONG_NAME, String.valueOf(threads))
                .addArgument(GenomicsDBImport.BATCHSIZE_ARG_NAME, String.valueOf(batchSize))
                .addArgument("L", IntervalUtils.locatableToString(INTERVAL));

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
            final VariantContext variant = new VariantContextBuilder("invented", contig, INTERVAL.getStart(), INTERVAL.getStart(), alleles)
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
                .addArgument(GenomicsDBImport.WORKSPACE_ARG_NAME, createTempDir("workspace").getAbsolutePath())
                .addArgument("L",  IntervalUtils.locatableToString(INTERVAL));
        runCommandLine(args);
    }

    @Test(expectedExceptions = CommandLineException.MissingArgument.class)
    public void testRequireOneOfVCFOrSampleNameFile(){
        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addArgument(GenomicsDBImport.WORKSPACE_ARG_NAME, createTempDir("workspace").getAbsolutePath())
                .addArgument("L", "1:1-10");

        runCommandLine(args);
    }

    @Test
    public void testCommandIncludedInOutputHeader() throws IOException {
        final List<String> vcfInputs = LOCAL_GVCFS;
        final String workspace = createTempDir("genomicsdb-tests").getAbsolutePath() + "/workspace";

        writeToGenomicsDB(vcfInputs, INTERVAL, workspace, 0, false, 0, 1);
        try(final GenomicsDBFeatureReader<VariantContext, PositionalBufferedStream> genomicsDBFeatureReader =
                new GenomicsDBFeatureReader<>(
                        new File(workspace, GenomicsDBConstants.DEFAULT_VIDMAP_FILE_NAME).getAbsolutePath(),
                        new File(workspace, GenomicsDBConstants.DEFAULT_CALLSETMAP_FILE_NAME).getAbsolutePath(),
                        workspace,
                        GenomicsDBConstants.DEFAULT_ARRAY_NAME,
                        b38_reference_20_21,
			new File(workspace, GenomicsDBConstants.DEFAULT_VCFHEADER_FILE_NAME).getAbsolutePath(), 
			new BCF2Codec()))
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

        writeToGenomicsDB(Arrays.asList(GENOMICSDB_TEST_DIR + "testHeaderContigLineSorting1.g.vcf", GENOMICSDB_TEST_DIR + "testHeaderContigLineSorting2.g.vcf"),
                          new SimpleInterval("chr20", 17959479, 17959479), workspace, 0, false, 0, 1);

        try ( final GenomicsDBFeatureReader<VariantContext, PositionalBufferedStream> genomicsDBFeatureReader =
                new GenomicsDBFeatureReader<>(
                        new File(workspace, GenomicsDBConstants.DEFAULT_VIDMAP_FILE_NAME).getAbsolutePath(),
                        new File(workspace, GenomicsDBConstants.DEFAULT_CALLSETMAP_FILE_NAME).getAbsolutePath(),
                        workspace,
                        GenomicsDBConstants.DEFAULT_ARRAY_NAME,
                        b38_reference_20_21, 
			new File(workspace, GenomicsDBConstants.DEFAULT_VCFHEADER_FILE_NAME).getAbsolutePath(),
			new BCF2Codec());

             final AbstractFeatureReader<VariantContext, LineIterator> inputGVCFReader =
                      AbstractFeatureReader.getFeatureReader(GENOMICSDB_TEST_DIR + "testHeaderContigLineSorting1.g.vcf", new VCFCodec(), true);
        ) {
            final SAMSequenceDictionary dictionaryFromGenomicsDB = ((VCFHeader)genomicsDBFeatureReader.getHeader()).getSequenceDictionary();
            final SAMSequenceDictionary dictionaryFromInputGVCF =  ((VCFHeader)inputGVCFReader.getHeader()).getSequenceDictionary();

            Assert.assertEquals(dictionaryFromGenomicsDB, dictionaryFromInputGVCF, "Sequence dictionary from GenomicsDB does not match original sequence dictionary from input GVCF");
        }

    }

    private static GenomicsDBFeatureReader<VariantContext, PositionalBufferedStream> getGenomicsDBFeatureReader(final String workspace, final String reference) throws IOException {
        return new GenomicsDBFeatureReader<>(
                new File(workspace, GenomicsDBConstants.DEFAULT_VIDMAP_FILE_NAME).getAbsolutePath(),
                new File(workspace, GenomicsDBConstants.DEFAULT_CALLSETMAP_FILE_NAME).getAbsolutePath(),
                workspace,
                GenomicsDBConstants.DEFAULT_ARRAY_NAME,
                reference,
                null,
                new BCF2Codec());
    }

    @Test(expectedExceptions = GenomicsDBImport.UnableToCreateGenomicsDBWorkspace.class)
    public void testYouCantWriteIntoAnExistingDirectory(){
        // this actually creates the directory on disk, not just the file name.
        final String workspace = createTempDir("workspace").getAbsolutePath();
        writeToGenomicsDB(LOCAL_GVCFS, INTERVAL, workspace, 0, false, 0, 1);
    }
}
