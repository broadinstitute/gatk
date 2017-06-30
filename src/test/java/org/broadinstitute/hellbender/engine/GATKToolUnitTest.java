package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.SAMFormatException;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.ValidationStringency;
import htsjdk.tribble.Feature;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFIDHeaderLine;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineArgumentParser;
import org.broadinstitute.barclay.argparser.CommandLineParser;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.cmdline.TestProgramGroup;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.reference.ReferenceUtils;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.test.TestResources;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.Iterator;
import java.io.IOException;
import java.util.Set;

public final class GATKToolUnitTest extends BaseTest{

    public static final String bqsrTestDir = TestResources.toolsTestDir + "BQSR/";

    public static final String BQSR_WGS_B37_CH20_21_10M_100_CRAM = bqsrTestDir +
            "CEUTrio.HiSeq.WGS.b37.NA12878.20.21.10m-10m100.cram";

    public static final String hg19MicroDictFileName = ReferenceUtils.getFastaDictionaryFileName(TestResources.hg19MicroReference);
    public static final String hg19MiniDictFileName  = ReferenceUtils.getFastaDictionaryFileName(TestResources.hg19MiniReference);
    public static final String v37_chr17DictFileName = ReferenceUtils.getFastaDictionaryFileName(TestResources.v37_chr17_1Mb_Reference);
    public static final String b37_20_21DictFile     = ReferenceUtils.getFastaDictionaryFileName(TestResources.b37_reference_20_21);

    @CommandLineProgramProperties(
            summary = "TestGATKToolWithSequenceDictionary",
            oneLineSummary = "TestGATKToolWithSequenceDictionary",
            programGroup = TestProgramGroup.class
    )
    private static final class TestGATKToolWithSequenceDictionary extends GATKTool {

        @Override
        public void traverse() {
            //no op
        }
    }

    @CommandLineProgramProperties(
            summary = "TestGATKToolWithReads",
            oneLineSummary = "TestGATKToolWithReads",
            programGroup = TestProgramGroup.class
    )
    private static final class TestGATKToolWithReads extends GATKTool{

        @Override
        public boolean requiresReads() {
            return true;
        }

        @Override
        public void traverse() {
            //no op
        }
    }

    @CommandLineProgramProperties(
            summary = "TestGATKToolWithFeatures",
            oneLineSummary = "TestGATKToolWithFeatures",
            programGroup = TestProgramGroup.class
    )
    private static final class TestGATKToolWithFeatures extends GATKTool{

        @Argument(fullName="mask", shortName="mask", doc="Input mask", optional=true)
        public FeatureInput<Feature> mask;

        @Override
        public boolean requiresFeatures() {
            return true;
        }

        @Override
        public void traverse() {
            //no op
        }
    }

    @CommandLineProgramProperties(
            summary = "TestGATKToolValidationStringency",
            oneLineSummary = "TestGATKToolValidationStringency",
            programGroup = TestProgramGroup.class
    )
    private static final class TestGATKToolValidationStringency extends GATKTool {
        private int count = 0;

        @Override
        public boolean requiresReads() {
            return true;
        }

        @Override
        public void traverse() {
            Iterator<GATKRead> iterator = reads.iterator();
            while (iterator.hasNext()) {
                GATKRead read = iterator.next();
                count++;
            }
        }

        public int getCount() { return count; }
    }

    @CommandLineProgramProperties(
            summary = "TestGATKToolWithNothing",
            oneLineSummary = "TestGATKToolWithNothing",
            programGroup = TestProgramGroup.class
    )
    private static final class TestGATKToolWithNothing extends GATKTool{

        @Override
        public void traverse() {
            //no op
        }
    }

    @CommandLineProgramProperties(
            summary = "TestGATKToolWithVariants",
            oneLineSummary = "TestGATKToolWithVariants",
            programGroup = TestProgramGroup.class
    )
    private static final class TestGATKToolWithVariants extends GATKTool{

        @Argument(fullName="output", shortName="out", doc="Input variants", optional=true)
        public File out;

        @Override
        public SAMSequenceDictionary getBestAvailableSequenceDictionary() {
            return new SAMSequenceDictionary();
        }

        @Override
        public boolean requiresFeatures() {
            return true;
        }

        @Override
        public void traverse() {
            //no op
        }
    }

    @DataProvider
    public Object[][] sequenceDictionaryTestValuesCompatible() {

        return new Object[][] {
                { v37_chr17DictFileName, "--input", TestResources.NA12878_chr17_1k_BAM, null, null },
                { b37_20_21DictFile, "--input", BQSR_WGS_B37_CH20_21_10M_100_CRAM, "--reference", TestResources.b37_reference_20_21 },
                { hg19MiniDictFileName, "--reference", TestResources.hg19MicroReference, null, null },
                { hg19MicroDictFileName, "--reference", TestResources.hg19MiniReference, null, null },
                { v37_chr17DictFileName, "--reference", TestResources.v37_chr17_1Mb_Reference, null, null },
        };
    }

    @DataProvider
    public Object[][] sequenceDictionaryTestValuesIncompatible() {

        return new Object[][] {
                { v37_chr17DictFileName, "--input", TestResources.NA12878_20_21_WGS_bam, null, null },
                { b37_20_21DictFile, "--input", BQSR_WGS_B37_CH20_21_10M_100_CRAM, "--reference", TestResources.v37_chr17_1Mb_Reference },
                { b37_20_21DictFile, "--reference", TestResources.hg19MiniReference, null, null },
                { v37_chr17DictFileName, "--reference", TestResources.hg19_chr1_1M_Reference, null, null },
        };
    }

    private void testGATKToolWithSequenceDictionaryHelper(String masterSequenceDictionaryFile,
                                                         String otherSeqArg, String otherSequenceFile,
                                                         String cramRefArg, String cramRefFile) throws Exception {

        final GATKTool tool = new TestGATKToolWithSequenceDictionary();
        final CommandLineParser clp = new CommandLineArgumentParser(tool);

        final String[] args = (cramRefArg == null)
                ? new String[] {"--sequenceDictionary", masterSequenceDictionaryFile, otherSeqArg, otherSequenceFile }
                : new String[] {"--sequenceDictionary", masterSequenceDictionaryFile, otherSeqArg, otherSequenceFile, cramRefArg, cramRefFile };

        clp.parseArguments(System.out, args);

        // This method would throw if sequence dictionary validation failed.
        // Here we are testing that it throws when the input arguments are incompatible,
        // and that it functions normally when the arguments are correct.
        tool.onStartup();
    }

    @Test(dataProvider= "sequenceDictionaryTestValuesCompatible")
    public void TestGATKToolWithSequenceDictionaryOk(String masterSequenceDictionaryFile,
                                                     String otherSeqArg, String otherSequenceFile,
                                                     String cramRefArg, String cramRefFile) throws Exception {

        testGATKToolWithSequenceDictionaryHelper(masterSequenceDictionaryFile,
                                                 otherSeqArg,
                                                 otherSequenceFile,
                                                 cramRefArg,
                                                 cramRefFile);
    }

    @Test(dataProvider="sequenceDictionaryTestValuesIncompatible",
            expectedExceptions = UserException.IncompatibleSequenceDictionaries.class)
    public void TestGATKToolWithSequenceDictionaryException(String masterSequenceDictionaryFile,
                                                            String otherSeqArg, String otherSequenceFile,
                                                            String cramRefArg, String cramRefFile) throws Exception {

        testGATKToolWithSequenceDictionaryHelper(masterSequenceDictionaryFile,
                                                 otherSeqArg,
                                                 otherSequenceFile,
                                                 cramRefArg,
                                                 cramRefFile);
    }

    @DataProvider
    public Object[][] provideForGetMasterSequenceTest() {

        final String dictFileName = ReferenceUtils.getFastaDictionaryFileName(TestResources.v37_chr17_1Mb_Reference);
        final SAMSequenceDictionary dict = ReferenceUtils.loadFastaDictionary(new File(dictFileName));

        return new Object[][] {
                { null, TestResources.NA12878_chr17_1k_BAM, null },
                { dictFileName, TestResources.NA12878_chr17_1k_BAM, dict },
        };
    }

    @Test(dataProvider = "provideForGetMasterSequenceTest")
    public void testGetMasterSequenceDictionary(String masterSequenceFileName, String inputFileName, SAMSequenceDictionary expectedDict) {
        final GATKTool tool = new TestGATKToolWithSequenceDictionary();
        final CommandLineParser clp = new CommandLineArgumentParser(tool);

        final String[] args = (masterSequenceFileName == null)
                ? new String[] { "--input", inputFileName, }
                : new String[] { "--sequenceDictionary", masterSequenceFileName, "--input", inputFileName, };

        clp.parseArguments(System.out, args);

        // This initializes our tool
        tool.onStartup();

        // Make sure that the master sequence dictionary dictionary is the expected dictionary
        Assert.assertEquals( tool.getMasterSequenceDictionary(), expectedDict );
    }

    @DataProvider
    public Object[][] provideMultipleSequenceDictionaries() {

        final SAMSequenceDictionary v37_chr17Dict = ReferenceUtils.loadFastaDictionary(new File(v37_chr17DictFileName));
        final SAMSequenceDictionary b37_20_21Dict = ReferenceUtils.loadFastaDictionary(new File(b37_20_21DictFile));
        final SAMSequenceDictionary hg19MiniDict = ReferenceUtils.loadFastaDictionary(new File(hg19MiniDictFileName));
        final SAMSequenceDictionary hg19MicroDict = ReferenceUtils.loadFastaDictionary(new File(hg19MicroDictFileName));

        return new Object[][] {
                { v37_chr17DictFileName, "--input", TestResources.NA12878_chr17_1k_BAM, null, null, v37_chr17Dict },
                { b37_20_21DictFile, "--input", BQSR_WGS_B37_CH20_21_10M_100_CRAM, "--reference", TestResources.b37_reference_20_21, b37_20_21Dict },
                { hg19MiniDictFileName, "--reference", TestResources.hg19MicroReference, null, null, hg19MiniDict },
                { hg19MicroDictFileName, "--reference", TestResources.hg19MiniReference, null, null, hg19MicroDict },
        };
    }

    @Test(dataProvider = "provideMultipleSequenceDictionaries")
    public void testGetBestAvailableSequenceDictionaryWithMasterDictionary(String sequenceFileName,
                                                               String otherSeqArg, String otherSequenceFile,
                                                               String cramRefArg, String cramRefFile,
                                                               SAMSequenceDictionary expectedDict) {

        final GATKTool tool = new TestGATKToolWithSequenceDictionary();
        final CommandLineParser clp = new CommandLineArgumentParser(tool);

        final String[] args = (cramRefArg == null)
            ? new String[] { "--sequenceDictionary", sequenceFileName, otherSeqArg, otherSequenceFile, }
            : new String[] { "--sequenceDictionary", sequenceFileName, otherSeqArg, otherSequenceFile, cramRefArg, cramRefFile };

        clp.parseArguments(System.out, args);

        // This initializes our tool
        tool.onStartup();

        // Make sure that the best available dictionary is the expected dictionary
        Assert.assertEquals( tool.getBestAvailableSequenceDictionary(), expectedDict );
    }

    @Test
    public void testReadsHeader() throws Exception {
        final GATKTool tool = new TestGATKToolWithReads();
        final CommandLineParser clp = new CommandLineArgumentParser(tool);
        final File bamFile = new File(TestResources.publicTestDir + "org/broadinstitute/hellbender/engine/reads_data_source_test1.bam");
        final String[] args = {"-I", bamFile.getCanonicalPath()};
        clp.parseArguments(System.out, args);
        tool.onStartup();
        final SAMFileHeader headerForReads = tool.getHeaderForReads();

        final SamReaderFactory factory = SamReaderFactory.makeDefault()    //read the file directly and compare headers
                .validationStringency(ValidationStringency.SILENT);
        try(SamReader samReader = factory.open(bamFile)) {
            final SAMFileHeader samFileHeader = samReader.getFileHeader();
            Assert.assertEquals(headerForReads, samFileHeader);
        }
        tool.doWork();
        tool.onShutdown();
    }

    @Test
    public void testFeaturesHeader() throws Exception {
        final TestGATKToolWithFeatures tool = new TestGATKToolWithFeatures();
        final CommandLineParser clp = new CommandLineArgumentParser(tool);
        final File vcfFile = new File(TestResources.publicTestDir + "org/broadinstitute/hellbender/engine/feature_data_source_test_with_bigHeader.vcf");
        final String[] args = {"--mask", vcfFile.getCanonicalPath()};
        clp.parseArguments(System.out, args);
        tool.onStartup();
        final Object headerForFeatures = tool.getHeaderForFeatures(tool.mask);
        Assert.assertTrue(headerForFeatures instanceof VCFHeader);
        final VCFHeader  vcfheaderForFeatures = (VCFHeader) headerForFeatures;

        try(final VCFFileReader vcfReader = new VCFFileReader(vcfFile, false)){  //read the file directly and compare headers
            final VCFHeader vcfFileHeader= vcfReader.getFileHeader();
            Assert.assertEquals(vcfheaderForFeatures.getGenotypeSamples(), vcfFileHeader.getGenotypeSamples());
            Assert.assertEquals(vcfheaderForFeatures.getInfoHeaderLines(), vcfFileHeader.getInfoHeaderLines());
            Assert.assertEquals(vcfheaderForFeatures.getFormatHeaderLines(), vcfFileHeader.getFormatHeaderLines());
            Assert.assertEquals(vcfheaderForFeatures.getFilterLines(), vcfFileHeader.getFilterLines());
            Assert.assertEquals(vcfheaderForFeatures.getContigLines(), vcfFileHeader.getContigLines());
            Assert.assertEquals(vcfheaderForFeatures.getOtherHeaderLines(), vcfFileHeader.getOtherHeaderLines());
        }
        tool.doWork();
        tool.onShutdown();
    }

    @Test
    public void testAllowLexicographicallySortedVariantHeader() throws Exception {
        final GATKTool tool = new TestGATKToolWithFeatures();
        final CommandLineParser clp = new CommandLineArgumentParser(tool);
        final File vcfFile = new File(TestResources.publicTestDir + "org/broadinstitute/hellbender/engine/lexicographically_sorted_dict.vcf");
        final String[] args = {"--mask", vcfFile.getCanonicalPath() };
        clp.parseArguments(System.out, args);

        // This method would throw if sequence dictionary validation failed. Here we are testing
        // that it does not throw despite the lexicographically-sorted sequence dictionary in the vcf.
        tool.onStartup();
    }

    @Test(expectedExceptions = UserException.MissingReference.class)
    public void testNonExistentReferenceFile() throws Exception {
        final TestGATKToolWithFeatures tool = new TestGATKToolWithFeatures();
        final CommandLineParser clp = new CommandLineArgumentParser(tool);
        final String[] args = {"--reference", BaseTest.getSafeNonExistentFile("NonExistentReferenceFile.fasta").getAbsolutePath()};
        clp.parseArguments(System.out, args);
        tool.onStartup();
    }

    @Test(expectedExceptions = UserException.IncompatibleSequenceDictionaries.class)
    public void testDisallowLexicographicallySortedVariantHeader_ifClashWithReference() throws Exception {
        final GATKTool tool = new TestGATKToolWithFeatures();
        final CommandLineParser clp = new CommandLineArgumentParser(tool);
        final File vcfFile = new File(TestResources.publicTestDir + "org/broadinstitute/hellbender/engine/lexicographically_sorted_dict.vcf");
        final String[] args = {"--mask", vcfFile.getCanonicalPath(),
                "--reference", TestResources.hg19MiniReference};
        clp.parseArguments(System.out, args);

        // This method throws despite the lexicographically-sorted sequence dictionary in the vcf.
        //This is because the reference sequence dictionary clashes with the one from the VCF.
        tool.onStartup();
    }

    @DataProvider(name="validationStringency")
    public Object[][] validationStringency() {
        return new Object[][]{
                {TestResources.NA12878_chr17_1k_CRAM, TestResources.v37_chr17_1Mb_Reference, 493}
        };
    }

    private void testValidationStringency(
            final String bamFileName,
            final String referenceFileName,
            final String validationStringency, final int count) throws SAMFormatException
    {
        final TestGATKToolValidationStringency tool = new TestGATKToolValidationStringency();
        final CommandLineParser clp = new CommandLineArgumentParser(tool);
        final File bamFile = new File(TestResources.NA12878_chr17_1k_CRAM);
        final File refFile = new File(TestResources.v37_chr17_1Mb_Reference);
        final String[] args = {
                "-I", bamFileName,
                "-R", referenceFileName,
                "-VS", validationStringency
        };

        clp.parseArguments(System.out, args);
        tool.onStartup();
        tool.doWork();
        tool.onShutdown();

        Assert.assertEquals(tool.getCount(), count);
    }

    @Test(dataProvider = "validationStringency", expectedExceptions=SAMFormatException.class)
    public void testReadsValidationStringencyStrict(final String bamFileName, final String referenceFileName, int count) throws Exception {
        testValidationStringency(bamFileName, referenceFileName, "STRICT", count);
    }

    @Test(dataProvider = "validationStringency")
    public void testReadsValidationStringencyLenient(final String bamFileName, final String referenceFileName, int count) throws Exception {
        testValidationStringency(bamFileName, referenceFileName, "LENIENT", count);
    }

    @Test(dataProvider = "validationStringency")
    public void testReadsValidationStringencySilent(final String bamFileName, final String referenceFileName, int count) throws Exception {
        testValidationStringency(bamFileName, referenceFileName, "SILENT", count);
    }

    @Test
    public void testBestSequenceDictionary_fromReads() throws Exception {
        final GATKTool tool = new TestGATKToolWithReads();
        final CommandLineParser clp = new CommandLineArgumentParser(tool);
        final File bamFile = new File(TestResources.publicTestDir + "org/broadinstitute/hellbender/engine/reads_data_source_test1.bam");
        final String[] args = {"-I", bamFile.getCanonicalPath()};
        clp.parseArguments(System.out, args);
        tool.onStartup();
        //read the dict back in and compare to bam dict
        final SAMSequenceDictionary toolDict = tool.getBestAvailableSequenceDictionary();
        try(final SamReader open = SamReaderFactory.makeDefault().open(bamFile)) {
            final SAMSequenceDictionary bamDict = open.getFileHeader().getSequenceDictionary();
            toolDict.assertSameDictionary(bamDict);
            bamDict.assertSameDictionary(toolDict);
            Assert.assertEquals(toolDict, bamDict);
        }
    }

    @Test
    public void testBestSequenceDictionary_fromReadsAndReference() throws Exception {
        final GATKTool tool = new TestGATKToolWithReads();
        final CommandLineParser clp = new CommandLineArgumentParser(tool);
        final File bamFile = new File(TestResources.publicTestDir + "org/broadinstitute/hellbender/engine/reads_data_source_test1.bam");
        final String fastaFile =   TestResources.hg19MiniReference;
        final String[] args = {"-I", bamFile.getCanonicalPath(),
                               "-R", fastaFile};
        clp.parseArguments(System.out, args);
        tool.onStartup();
        //read the dict back in and compare to reference dict
        final SAMSequenceDictionary toolDict = tool.getBestAvailableSequenceDictionary();
        final SAMSequenceDictionary fastaDict = new IndexedFastaSequenceFile(new File(fastaFile)).getSequenceDictionary();
        toolDict.assertSameDictionary(fastaDict);
        fastaDict.assertSameDictionary(toolDict);

        Assert.assertEquals(toolDict, fastaDict);
    }

    @Test
    public void testBestSequenceDictionary_fromVariants() throws Exception {
        final GATKTool tool = new TestGATKToolWithFeatures();
        final CommandLineParser clp = new CommandLineArgumentParser(tool);
        final File vcfFile = new File(TestResources.publicTestDir + "org/broadinstitute/hellbender/engine/feature_data_source_test_withSequenceDict.vcf");
        final String[] args = {"--mask", vcfFile.getCanonicalPath()};
        clp.parseArguments(System.out, args);
        tool.onStartup();
        //read the dict back in and compare to vcf dict
        final SAMSequenceDictionary toolDict = tool.getBestAvailableSequenceDictionary();
        try(final VCFFileReader reader = new VCFFileReader(vcfFile)) {
            final SAMSequenceDictionary vcfDict = reader.getFileHeader().getSequenceDictionary();
            toolDict.assertSameDictionary(vcfDict);
            vcfDict.assertSameDictionary(toolDict);
            Assert.assertEquals(toolDict, vcfDict);
        }
    }

    @Test
    public void testBestSequenceDictionary_fromNothing() throws Exception {
        final GATKTool tool = new TestGATKToolWithNothing();
        final CommandLineParser clp = new CommandLineArgumentParser(tool);
        final String[] args = {};
        clp.parseArguments(System.out, args);
        tool.onStartup();
        //read the dict back in and assert that it's null
        final SAMSequenceDictionary toolDict = tool.getBestAvailableSequenceDictionary();
        Assert.assertNull(toolDict);
    }

    @DataProvider(name="createVCFWriterData")
    public Object[][] createVCFWriterData() {
        return new Object[][]{
                {new File(TestResources.publicTestDir, "org/broadinstitute/hellbender/engine/example_variants.vcf"), ".vcf", ".idx", true, true},
                {new File(TestResources.publicTestDir, "org/broadinstitute/hellbender/engine/example_variants.vcf"), ".vcf", ".idx", false, true},
                {new File(TestResources.publicTestDir, "org/broadinstitute/hellbender/engine/example_variants.vcf"), ".vcf", ".idx", true, false},
                {new File(TestResources.publicTestDir, "org/broadinstitute/hellbender/engine/example_variants.vcf"), ".bcf", ".idx", true, true},
                {new File(TestResources.publicTestDir, "org/broadinstitute/hellbender/engine/example_variants.vcf"), ".bcf", ".idx", false, true},
                {new File(TestResources.publicTestDir, "org/broadinstitute/hellbender/engine/example_variants.vcf"), ".bcf", ".idx", true, false},
                {new File(TestResources.publicTestDir, "org/broadinstitute/hellbender/engine/example_variants.vcf"), ".vcf.bgz", ".tbi", true, true},
                {new File(TestResources.publicTestDir, "org/broadinstitute/hellbender/engine/example_variants.vcf"), ".vcf.gz", ".tbi", false, true},
                {new File(TestResources.publicTestDir, "org/broadinstitute/hellbender/engine/example_variants.vcf"), ".vcf.bgz", ".tbi", true, false}
        };
    }

    @Test(dataProvider = "createVCFWriterData")
    public void testCreateVCFWriterDefaults(
            final File inputFile,           // unused
            final String outputExtension,
            final String indexExtension,
            final boolean createIndex,      // unused
            final boolean createMD5         // unused
    ) throws IOException {

        // create a writer and make sure the default index/md5 params are honored
        final TestGATKToolWithVariants tool = createTestVariantTool(null);

        final File tmpDir = createTempDir("createVCFTest");
        final File outputFile = new File(tmpDir.getAbsolutePath(), "createVCFTest" + outputExtension);

        final VariantContextWriter writer = tool.createVCFWriter(outputFile);
        writer.close();

        final File outFileIndex = new File(outputFile.getAbsolutePath() + indexExtension);
        final File outFileMD5 = new File(outputFile.getAbsolutePath() + ".md5");

        Assert.assertTrue(outputFile.exists(), "No output file was not created");
        Assert.assertTrue(outFileIndex.exists(), "The index file was not created");
        Assert.assertFalse(outFileMD5.exists(), "An md5 file was created and should not have been");
    }

    private TestGATKToolWithVariants createTestVariantTool(final String args[]) {
        final TestGATKToolWithVariants tool = new TestGATKToolWithVariants();
        if (null != args) {
            final CommandLineParser clp = new CommandLineArgumentParser(tool);
            clp.parseArguments(System.out, args);
        }
        return tool;
    }

    @Test(dataProvider = "createVCFWriterData")
    public void testCreateVCFWriterWithOptions(
            final File inputFile,
            final String outputExtension,
            final String indexExtension,
            final boolean createIndex,
            final boolean createMD5) throws IOException {

        // create a writer and make sure the requested index/md5 params are honored
        final TestGATKToolWithVariants tool = new TestGATKToolWithVariants();

        final File outputFile = setupVCFWriter(inputFile, outputExtension, tool, createIndex, createMD5, false);

        final VariantContextWriter writer = tool.createVCFWriter(outputFile);
        writer.close();

        final File outFileIndex = new File(outputFile.getAbsolutePath() + indexExtension);
        final File outFileMD5 = new File(outputFile.getAbsolutePath() + ".md5");

        Assert.assertTrue(outputFile.exists(), "No output file was not created");
        Assert.assertEquals(outFileIndex.exists(), createIndex, "The createIndex argument was not honored");
        Assert.assertEquals(outFileMD5.exists(), createMD5, "The createMD5 argument was not honored");
    }

    @DataProvider(name="createVCFWriterLenientData")
    public Object[][] createVCFWriterLenientData() {
        return new Object[][]{
                {new File(TestResources.publicTestDir, "org/broadinstitute/hellbender/engine/example_variants.vcf"), ".vcf", ".idx", true, true},
                {new File(TestResources.publicTestDir, "org/broadinstitute/hellbender/engine/example_variants.vcf"), ".vcf", ".idx", false, true},
                {new File(TestResources.publicTestDir, "org/broadinstitute/hellbender/engine/example_variants.vcf"), ".vcf", ".idx", true, false}
        };
    }

    @Test(dataProvider = "createVCFWriterLenientData")
    public void testCreateVCFWriterLenientTrue(
            final File inputFile,
            final String outputExtension,
            final String indexExtension,
            final boolean createIndex,
            final boolean createMD5) throws IOException {
        final TestGATKToolWithVariants tool = new TestGATKToolWithVariants();

        // verify lenient==true is honored by writing a bad attribute
        final File outputFile = setupVCFWriter(inputFile, outputExtension, tool, createIndex, createMD5, true);

        try (VariantContextWriter writer = tool.createVCFWriter(outputFile)) {
            writeHeaderAndBadVariant(writer); // write bad attribute succeed with lenient set
        }

        final File outFileIndex = new File(outputFile.getAbsolutePath() + indexExtension);
        final File outFileMD5 = new File(outputFile.getAbsolutePath() + ".md5");

        Assert.assertTrue(outputFile.exists(), "No output file was not created");
        Assert.assertEquals(outFileIndex.exists(), createIndex, "The createIndex argument was not honored");
        Assert.assertEquals(outFileMD5.exists(), createMD5, "The createMD5 argument was not honored");
    }

    @Test(dataProvider = "createVCFWriterLenientData", expectedExceptions = IllegalStateException.class)
    public void testCreateVCFWriterLenientFalse(
            final File inputFile,
            final String outputExtension,
            final String indexExtension, // unused
            final boolean createIndex,
            final boolean createMD5) throws IOException {

        // verify lenient==false is honored by writing a bad attribute
        final TestGATKToolWithVariants tool = new TestGATKToolWithVariants();
        final File outputFile = setupVCFWriter(inputFile, outputExtension, tool, createIndex, createMD5, false);

        try (VariantContextWriter writer = tool.createVCFWriter(outputFile)) {
            writeHeaderAndBadVariant(writer); // throws due to bad attribute
        }
    }

    @CommandLineProgramProperties(
            summary = "TestGATKVariantToolWithNoSequenceDictionary",
            oneLineSummary = "TestGATKVariantToolWithNoSequenceDictionary",
            programGroup = TestProgramGroup.class
    )
    private static final class TestGATKVariantToolWithNoSequenceDictionary extends GATKTool{
        @Argument(fullName="output", shortName="out", doc="Input variants", optional=true)
        public File out;

        @Override
        public SAMSequenceDictionary getBestAvailableSequenceDictionary() {return null;}

        @Override
        public boolean requiresFeatures() {return true;}

        @Override
        public void traverse() {}
    }

    @Test(dataProvider = "createVCFWriterData")
    public void testCreateVCFWriterWithNoSequenceDictionary(
            final File inputFile,
            final String outputExtension,
            final String indexExtension,
            final boolean createIndex,
            final boolean createMD5) throws IOException
    {
        // verify that a null sequence dictionary still results in a file, but with no index
        final TestGATKVariantToolWithNoSequenceDictionary tool = new TestGATKVariantToolWithNoSequenceDictionary();
        final File outputFile = setupVCFWriter(inputFile, outputExtension, tool, createIndex, createMD5, false);

        final VariantContextWriter writer = tool.createVCFWriter(outputFile);
        writer.close();

        final File outFileIndex = new File(outputFile.getAbsolutePath() + indexExtension);
        final File outFileMD5 = new File(outputFile.getAbsolutePath() + ".md5");

        Assert.assertTrue(outputFile.exists(), "No output file was not created");
        Assert.assertEquals(outFileIndex.exists(), false, "An index file should not have been created"); // always false with no seq dictionary
        Assert.assertEquals(outFileMD5.exists(), createMD5, "The createMD5 argument was not honored");
    }

    private File setupVCFWriter(
            final File inputFile,
            final String outputExtension,
            final GATKTool tool,
            final boolean createIndex,
            final boolean createMD5,
            final boolean lenient) throws IOException
    {
        final File tmpDir = createTempDir("createVCFTest");
        final File outputFile = new File(tmpDir.getAbsolutePath(), "createVCFTest" + outputExtension);

        ArgumentsBuilder args = new ArgumentsBuilder();

        args.addInput(inputFile);
        args.addOutput(outputFile);
        args.add("--createOutputVariantIndex"); args.add(Boolean.toString(createIndex));
        args.add("--createOutputVariantMD5"); args.add(Boolean.toString(createMD5));
        if (lenient) {
            args.add("--lenient");
        }

        final CommandLineParser clp = new CommandLineArgumentParser(tool);
        clp.parseArguments(System.out, args.getArgsArray());

        return outputFile;
    }

    @Test
    public void testGetDefaultToolVCFHeaderLines() throws IOException {
        final TestGATKToolWithFeatures tool = new TestGATKToolWithFeatures();
        final File vcfFile = new File(TestResources.publicTestDir + "org/broadinstitute/hellbender/engine/feature_data_source_test_with_bigHeader.vcf");
        final String[] args = {"--mask", vcfFile.getCanonicalPath(), "--" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "true"};
        tool.instanceMain(args);

        Set<VCFHeaderLine> stdHeaderLines = tool.getDefaultToolVCFHeaderLines();
        VCFHeader hdr = new VCFHeader(stdHeaderLines);

        VCFHeaderLine sourceLine = hdr.getOtherHeaderLine("source");
        Assert.assertEquals(sourceLine.getValue(), tool.getClass().getSimpleName());

        VCFIDHeaderLine commandLine = (VCFIDHeaderLine) hdr.getOtherHeaderLine("GATKCommandLine");
        Assert.assertEquals(commandLine.getID(), tool.getClass().getSimpleName());

        String commandLineString = commandLine.toString();
        assertContains(commandLineString,"CommandLine=");
        assertContains(commandLineString,"Version=");
        assertContains(commandLineString,"Date=");
    }

    private void writeHeaderAndBadVariant(final VariantContextWriter writer) {
        final VariantContextBuilder vcBuilder = new VariantContextBuilder(
                "chr1","1", 1, 1, Arrays.asList(Allele.create("A", true)));
        vcBuilder.attribute("fake", new Object());
        final VariantContext vc = vcBuilder.make();
        final VCFHeader vcfHeader = new VCFHeader();
        writer.writeHeader(vcfHeader);
        writer.add(vc);
    }

}
