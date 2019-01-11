package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.*;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.tribble.Feature;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFIDHeaderLine;
import org.broadinstitute.barclay.argparser.*;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.cmdline.GATKPlugin.GATKAnnotationPluginDescriptor;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.TestProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.annotator.Annotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.ClippingRankSumTest;
import org.broadinstitute.hellbender.tools.walkers.annotator.Coverage;
import org.broadinstitute.hellbender.tools.walkers.annotator.StandardAnnotation;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.reference.ReferenceUtils;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import javax.annotation.Nullable;
import java.io.File;
import java.io.IOException;
import java.util.*;

public final class GATKToolUnitTest extends GATKBaseTest {

    public static final String bqsrTestDir = toolsTestDir + "BQSR/";

    public static final String BQSR_WGS_B37_CH20_21_10M_100_CRAM = bqsrTestDir +
            "CEUTrio.HiSeq.WGS.b37.NA12878.20.21.10m-10m100.cram";

    public static final String hg19MicroDictFileName = ReferenceUtils.getFastaDictionaryFileName(hg19MicroReference);
    public static final String hg19MiniDictFileName  = ReferenceUtils.getFastaDictionaryFileName(hg19MiniReference);
    public static final String v37_chr17DictFileName = ReferenceUtils.getFastaDictionaryFileName(v37_chr17_1Mb_Reference);
    public static final String b37_20_21DictFile     = ReferenceUtils.getFastaDictionaryFileName(b37_reference_20_21);

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

        public List<SimpleInterval> getIntervals() {
            return intervalArgumentCollection.getIntervals(getBestAvailableSequenceDictionary());
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

        @Override
        public boolean useVariantAnnotations() {
            return true;
        }
    }

    @CommandLineProgramProperties(
            summary = "TestGATKToolWithDefaultAnnotations",
            oneLineSummary = "TestGATKToolWithDefaultAnnotations",
            programGroup = TestProgramGroup.class
    )
    private static final class TestGATKToolWithDefaultAnnotations extends GATKTool{

        @Override
        public void traverse() {
            //no op
        }

        @Override
        public boolean useVariantAnnotations() {
            return true;
        }

        @Override
        public List<Annotation> getDefaultVariantAnnotations() {
            return Collections.singletonList(new Coverage());
        }
    }

    @CommandLineProgramProperties(
            summary = "TestGATKToolWithDefaultAnnotations",
            oneLineSummary = "TestGATKToolWithDefaultAnnotations",
            programGroup = TestProgramGroup.class
    )
    private static final class TestGATKToolWithDefaultAnnotationGroups extends GATKTool{

        @Override
        public void traverse() {
            //no op
        }

        @Override
        public boolean useVariantAnnotations() {
            return true;
        }

        @Override
        public List<Class<? extends Annotation>> getDefaultVariantAnnotationGroups() {
            return Collections.singletonList(StandardAnnotation.class);
        }
    }

    @DataProvider
    public Object[][] sequenceDictionaryTestValuesCompatible() {

        return new Object[][] {
                { v37_chr17DictFileName, "--input", NA12878_chr17_1k_BAM, null, null },
                { b37_20_21DictFile, "--input", BQSR_WGS_B37_CH20_21_10M_100_CRAM, "--reference", b37_reference_20_21 },
                { hg19MiniDictFileName, "--reference", hg19MicroReference, null, null },
                { hg19MicroDictFileName, "--reference", hg19MiniReference, null, null },
                { v37_chr17DictFileName, "--reference", v37_chr17_1Mb_Reference, null, null },
        };
    }

    @DataProvider
    public Object[][] sequenceDictionaryTestValuesIncompatible() {

        return new Object[][] {
                { v37_chr17DictFileName, "--input", NA12878_20_21_WGS_bam, null, null },
                { b37_20_21DictFile, "--input", BQSR_WGS_B37_CH20_21_10M_100_CRAM, "--reference", v37_chr17_1Mb_Reference },
                { b37_20_21DictFile, "--reference", hg19MiniReference, null, null },
                { v37_chr17DictFileName, "--reference", hg19_chr1_1M_Reference, null, null },
        };
    }

    private void testGATKToolWithSequenceDictionaryHelper(String masterSequenceDictionaryFile,
                                                         String otherSeqArg, String otherSequenceFile,
                                                         String cramRefArg, String cramRefFile) throws Exception {

        final GATKTool tool = new TestGATKToolWithSequenceDictionary();
        final CommandLineParser clp = new CommandLineArgumentParser(tool);

        final String[] args = (cramRefArg == null)
                ? new String[] {"--" + StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME, masterSequenceDictionaryFile, otherSeqArg, otherSequenceFile }
                : new String[] {"--" + StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME, masterSequenceDictionaryFile, otherSeqArg, otherSequenceFile, cramRefArg, cramRefFile };

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

        final String dictFileName = ReferenceUtils.getFastaDictionaryFileName(v37_chr17_1Mb_Reference);
        final SAMSequenceDictionary dict = ReferenceUtils.loadFastaDictionary(new File(dictFileName));

        return new Object[][] {
                { null, NA12878_chr17_1k_BAM, null },
                { dictFileName, NA12878_chr17_1k_BAM, dict },
        };
    }

    @Test(dataProvider = "provideForGetMasterSequenceTest")
    public void testGetMasterSequenceDictionary(String masterSequenceFileName, String inputFileName, SAMSequenceDictionary expectedDict) {
        final GATKTool tool = new TestGATKToolWithSequenceDictionary();
        final CommandLineParser clp = new CommandLineArgumentParser(tool);

        final String[] args = (masterSequenceFileName == null)
                ? new String[] { "--input", inputFileName, }
                : new String[] { "--" + StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME, masterSequenceFileName, "--input", inputFileName, };

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
                { v37_chr17DictFileName, "--input", NA12878_chr17_1k_BAM, null, null, v37_chr17Dict },
                { b37_20_21DictFile, "--input", BQSR_WGS_B37_CH20_21_10M_100_CRAM, "--reference", b37_reference_20_21, b37_20_21Dict },
                { hg19MiniDictFileName, "--reference", hg19MicroReference, null, null, hg19MiniDict },
                { hg19MicroDictFileName, "--reference", hg19MiniReference, null, null, hg19MicroDict },
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
            ? new String[] { "--" + StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME, sequenceFileName, otherSeqArg, otherSequenceFile, }
            : new String[] { "--" + StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME, sequenceFileName, otherSeqArg, otherSequenceFile, cramRefArg, cramRefFile };

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
        final File bamFile = new File(publicTestDir + "org/broadinstitute/hellbender/engine/reads_data_source_test1.bam");
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
    public void testPicardIntervalList() throws Exception {
        final TestGATKToolWithReads tool = new TestGATKToolWithReads();
        final CommandLineParser clp = new CommandLineArgumentParser(tool);
        final File bamFile = new File(publicTestDir + "org/broadinstitute/hellbender/engine/reads_data_source_test1.bam");
        final File intervalsFile = new File(publicTestDir + "picard_intervals.list");
        final String[] args = {
                "-I", bamFile.getCanonicalPath(),
                "-L", intervalsFile.getCanonicalPath()
        };
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

        // ensure that the raw interval argument has not been expanded by Barclay, and that the post-merged
        // intervals list contains 3 intervals (there are 4 in the file; 2 get merged)
        Assert.assertEquals(tool.getIntervals().size(), 3);

        tool.onShutdown();
    }

    @Test
    public void testFeaturesHeader() throws Exception {
        final TestGATKToolWithFeatures tool = new TestGATKToolWithFeatures();
        final CommandLineParser clp = new CommandLineArgumentParser(tool);
        final File vcfFile = new File(publicTestDir + "org/broadinstitute/hellbender/engine/feature_data_source_test_with_bigHeader.vcf");
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
        final File vcfFile = new File(publicTestDir + "org/broadinstitute/hellbender/engine/lexicographically_sorted_dict.vcf");
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
        final String[] args = {"--reference", GATKBaseTest.getSafeNonExistentFile("NonExistentReferenceFile.fasta").getAbsolutePath()};
        clp.parseArguments(System.out, args);
        tool.onStartup();
    }

    @Test(expectedExceptions = UserException.IncompatibleSequenceDictionaries.class)
    public void testDisallowLexicographicallySortedVariantHeader_ifClashWithReference() throws Exception {
        final GATKTool tool = new TestGATKToolWithFeatures();
        final CommandLineParser clp = new CommandLineArgumentParser(tool);
        final File vcfFile = new File(publicTestDir + "org/broadinstitute/hellbender/engine/lexicographically_sorted_dict.vcf");
        final String[] args = {"--mask", vcfFile.getCanonicalPath(),
                "--reference", hg19MiniReference};
        clp.parseArguments(System.out, args);

        // This method throws despite the lexicographically-sorted sequence dictionary in the vcf.
        //This is because the reference sequence dictionary clashes with the one from the VCF.
        tool.onStartup();
    }

    @DataProvider(name="validationStringency")
    public Object[][] validationStringency() {
        return new Object[][]{
                {NA12878_chr17_1k_CRAM, v37_chr17_1Mb_Reference, 493}
        };
    }

    private void testValidationStringency(
            final String bamFileName,
            final String referenceFileName,
            final String validationStringency, final int count) throws SAMFormatException
    {
        final TestGATKToolValidationStringency tool = new TestGATKToolValidationStringency();
        final CommandLineParser clp = new CommandLineArgumentParser(tool);
        final File bamFile = new File(NA12878_chr17_1k_CRAM);
        final File refFile = new File(v37_chr17_1Mb_Reference);
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
        final File bamFile = new File(publicTestDir + "org/broadinstitute/hellbender/engine/reads_data_source_test1.bam");
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
        final File bamFile = new File(publicTestDir + "org/broadinstitute/hellbender/engine/reads_data_source_test1.bam");
        final String fastaFile =   hg19MiniReference;
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
        final File vcfFile = new File(publicTestDir + "org/broadinstitute/hellbender/engine/feature_data_source_test_withSequenceDict.vcf");
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
                {new File(publicTestDir, "org/broadinstitute/hellbender/engine/example_variants.vcf"), ".vcf", ".idx", true, true},
                {new File(publicTestDir, "org/broadinstitute/hellbender/engine/example_variants.vcf"), ".vcf", ".idx", false, true},
                {new File(publicTestDir, "org/broadinstitute/hellbender/engine/example_variants.vcf"), ".vcf", ".idx", true, false},
                {new File(publicTestDir, "org/broadinstitute/hellbender/engine/example_variants.vcf"), ".bcf", ".idx", true, true},
                {new File(publicTestDir, "org/broadinstitute/hellbender/engine/example_variants.vcf"), ".bcf", ".idx", false, true},
                {new File(publicTestDir, "org/broadinstitute/hellbender/engine/example_variants.vcf"), ".bcf", ".idx", true, false},
                {new File(publicTestDir, "org/broadinstitute/hellbender/engine/example_variants.vcf"), ".vcf.bgz", ".tbi", true, true},
                {new File(publicTestDir, "org/broadinstitute/hellbender/engine/example_variants.vcf"), ".vcf.gz", ".tbi", false, true},
                {new File(publicTestDir, "org/broadinstitute/hellbender/engine/example_variants.vcf"), ".vcf.bgz", ".tbi", true, false}
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

    @Test
    public void testMakeEmptyAnnotations() {
        final TestGATKToolWithVariants tool = createTestVariantTool(null);
        Collection<Annotation> annots = tool.makeVariantAnnotations();

        Assert.assertTrue(annots.isEmpty());
    }

    @Test
    public void testGetAllAnnotations() {
        String[] args = {"--"+StandardArgumentDefinitions.ENABLE_ALL_ANNOTATIONS};

        final TestGATKToolWithVariants tool = createTestVariantTool(args);
        Collection<Annotation> annots = tool.makeVariantAnnotations();

        ClassFinder finder = new ClassFinder();
        finder.find(GATKAnnotationPluginDescriptor.pluginPackageName, Annotation.class);

        Set<Class<?>> classes = finder.getConcreteClasses();
        Assert.assertFalse(classes.isEmpty());
        Assert.assertEquals(annots.size(),classes.size());
        for(Class<?> found : classes) {
            Assert.assertTrue(annots.stream().anyMatch(a -> a.getClass()==found));
        }
    }

    @Test
    public void testExcludeAnnotation(){
        String[] args = {"--"+StandardArgumentDefinitions.ENABLE_ALL_ANNOTATIONS, "-AX", "Coverage"};

        final TestGATKToolWithVariants tool = createTestVariantTool(args);
        Collection<Annotation> annots = tool.makeVariantAnnotations();

        // Asserting that the annotation was excluded
        Assert.assertFalse(annots.stream().anyMatch(a -> a.getClass()==Coverage.class));
        Assert.assertFalse(annots.isEmpty());
    }

    @Test
    public void testIncludeAnnotationGroups(){
        String[] args = {"-G", StandardAnnotation.class.getSimpleName()};

        final TestGATKToolWithVariants tool = createTestVariantTool(args);
        Collection<Annotation> annots = tool.makeVariantAnnotations();

        // Asserting that a standard annotation was included but not everything
        ClassFinder finder = new ClassFinder();
        finder.find(GATKAnnotationPluginDescriptor.pluginPackageName, StandardAnnotation.class);

        Set<Class<?>> classes = finder.getConcreteClasses();
        Assert.assertFalse(classes.isEmpty());
        Assert.assertEquals(annots.size(),classes.size());
        for(Class<?> found : classes) {
            Assert.assertTrue(annots.stream().anyMatch(a -> a.getClass()==found));
        }
    }


    @Test
    public void testIncludeAnnotation(){
        String[] args = {"-A", Coverage.class.getSimpleName()};

        final TestGATKToolWithVariants tool = createTestVariantTool(args);
        Collection<Annotation> annots = tool.makeVariantAnnotations();

        // Asserting coverage was added
        Assert.assertTrue(annots.stream().anyMatch(a -> a.getClass()==Coverage.class));
        Assert.assertTrue(annots.size() == 1);
    }

    @Test
    public void testMakeDefaultAnnotations() {
        String[] args = null;

        final TestGATKToolWithDefaultAnnotations tool = createTestVariantTool(new TestGATKToolWithDefaultAnnotations(), args);
        Collection<Annotation> annots = tool.makeVariantAnnotations();

        // Asserting coverage was added by default
        Assert.assertTrue(annots.stream().anyMatch(a -> a.getClass()==Coverage.class));
        Assert.assertTrue(annots.size() == 1);
    }

    @Test
    public void testMakeDefaultAnnotationGroups() {
        String[] args = null;

        final TestGATKToolWithDefaultAnnotationGroups tool = createTestVariantTool(new TestGATKToolWithDefaultAnnotationGroups(), args);
        Collection<Annotation> annots = tool.makeVariantAnnotations();

        ClassFinder finder = new ClassFinder();
        finder.find(GATKAnnotationPluginDescriptor.pluginPackageName, StandardAnnotation.class);

        Set<Class<?>> classes = finder.getConcreteClasses();
        Assert.assertFalse(classes.isEmpty());
        Assert.assertEquals(annots.size(),classes.size());
        for(Class<?> found : classes) {
            Assert.assertTrue(annots.stream().anyMatch(a -> a.getClass()==found));
        }

        Assert.assertFalse(annots.stream().anyMatch(a -> a.getClass()==StandardAnnotation.class));
        Assert.assertFalse(annots.stream().anyMatch(a -> a.getClass()==ClippingRankSumTest.class));
    }

    @Test
    public void testClearDefaultAnnotationsGroups() {
        String[] args = {"--"+StandardArgumentDefinitions.DISABLE_TOOL_DEFAULT_ANNOTATIONS};

        final TestGATKToolWithDefaultAnnotationGroups tool = createTestVariantTool(new TestGATKToolWithDefaultAnnotationGroups(), args);
        Collection<Annotation> annots = tool.makeVariantAnnotations();

        // Asserting that the standard annotation was not included when defaults are disabled
        Assert.assertFalse(annots.stream().anyMatch(a -> a.getClass()==Coverage.class));
        Assert.assertFalse(annots.stream().anyMatch(a -> a.getClass()==StandardAnnotation.class));
        Assert.assertFalse(annots.stream().anyMatch(a -> a.getClass()==ClippingRankSumTest.class));
    }

    @Test
    public void testClearDefaultAnnotations() {
        String[] args = {"--"+StandardArgumentDefinitions.DISABLE_TOOL_DEFAULT_ANNOTATIONS};

        final TestGATKToolWithDefaultAnnotations tool = createTestVariantTool(new TestGATKToolWithDefaultAnnotations(), args);
        Collection<Annotation> annots = tool.makeVariantAnnotations();

        // Asserting that the standard annotation was not included when defaults are disabled
        Assert.assertFalse(annots.stream().anyMatch(a -> a.getClass()==Coverage.class));
        Assert.assertFalse(annots.stream().anyMatch(a -> a.getClass()==StandardAnnotation.class));
        Assert.assertFalse(annots.stream().anyMatch(a -> a.getClass()==ClippingRankSumTest.class));
    }

    @Test
    public void testHelpWithAllPluginDescriptors() {
        // Smoke test to ensure that requesting help from plugin descriptors doesn't crash. Use a tool
        // (TestGATKToolWithVariants) that has both the read filter and annotation plugin descriptors enabled.
        String[] args = {"-h"};
        new TestGATKToolWithVariants().instanceMain(args);
    }

    private TestGATKToolWithVariants createTestVariantTool(final String args[]) {
       return createTestVariantTool(new TestGATKToolWithVariants(), args);
    }

    private <T extends GATKTool> T createTestVariantTool(final T tool, final String args[]) {

        final CommandLineParser clp = tool.getCommandLineParser();
        clp.parseArguments(System.out, args==null? new String[0] : args);

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
                {new File(publicTestDir, "org/broadinstitute/hellbender/engine/example_variants.vcf"), ".vcf", ".idx", true, true},
                {new File(publicTestDir, "org/broadinstitute/hellbender/engine/example_variants.vcf"), ".vcf", ".idx", false, true},
                {new File(publicTestDir, "org/broadinstitute/hellbender/engine/example_variants.vcf"), ".vcf", ".idx", true, false}
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
        args.add("--" + StandardArgumentDefinitions.CREATE_OUTPUT_VARIANT_INDEX_LONG_NAME); args.add(Boolean.toString(createIndex));
        args.add("--" + StandardArgumentDefinitions.CREATE_OUTPUT_VARIANT_MD5_LONG_NAME); args.add(Boolean.toString(createMD5));
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
        final File vcfFile = new File(publicTestDir + "org/broadinstitute/hellbender/engine/feature_data_source_test_with_bigHeader.vcf");
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

    final String baseVariants = packageRootTestDir + "engine/feature_data_source_test.vcf";

    @CommandLineProgramProperties(programGroup = TestProgramGroup.class, oneLineSummary = "GATKTool Intervals Test Walker", summary = "This is a test walker for GATKTool getTraversalIntervals")
    private static class TestIntervalWalker extends GATKTool {
        @Argument(fullName = StandardArgumentDefinitions.VARIANT_LONG_NAME, shortName = StandardArgumentDefinitions.VARIANT_SHORT_NAME, doc = "A VCF file containing variants", common = false, optional = false)
        public String drivingVariantFile;

        public TestIntervalWalker() {

        }

        @Override
        public void traverse() {

        }
    }

    @Test(expectedExceptions = UserException.class)
    public void testSequenceDictionaryRequiredForIntervalQuery() throws Exception {
        //This should have failed because no dictionary is provided
        final TestIntervalWalker tool = new TestIntervalWalker();
        tool.instanceMain(new String[]{
                "-V", baseVariants,
                "-L", "1:21-21"
        });
    }

    @DataProvider(name = "TestGetTraversalIntervalsProvider")
    public Object[][] getTestGetTraversalIntervalsProvider() {
        return new Object[][]{
                {"1:21-21", 1, hg19MiniReference},
                {null, 4, hg19MiniReference},
                {null, null, null}
        };
    }

    @Test(dataProvider = "TestGetTraversalIntervalsProvider")
    public void testGetTraversalIntervals(@Nullable String intervals, Integer expected, String ref) {
        List<String> args = new ArrayList<>(Arrays.asList(
                "-V", baseVariants
        ));

        if (ref != null) {
            args.add("-R");
            args.add(ref);
        }

        if (intervals != null) {
            args.add("-L");
            args.add(intervals);
        }

        final TestIntervalWalker tool = new TestIntervalWalker();
        tool.instanceMain(args.toArray(new String[args.size()]));

        Assert.assertEquals(tool.getTraversalIntervals() == null ? null : tool.getTraversalIntervals().size(), expected);
    }
}
