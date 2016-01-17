package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMFormatException;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.tribble.Feature;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineParser;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Iterator;

public final class GATKToolUnitTest extends BaseTest{

    @CommandLineProgramProperties(
            summary = "TestGATKToolWithReads",
            oneLineSummary = "TestGATKToolWithReads",
            programGroup = ReadProgramGroup.class
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
            summary = "TestGATKToolWithReads",
            oneLineSummary = "TestGATKToolWithReads",
            programGroup = ReadProgramGroup.class
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
            programGroup = ReadProgramGroup.class
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

    @Test
    public void testReadsHeader() throws Exception {
        final GATKTool tool = new TestGATKToolWithReads();
        final CommandLineParser clp = new CommandLineParser(tool);
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
    public void testFeaturesHeader() throws Exception {
        final TestGATKToolWithFeatures tool = new TestGATKToolWithFeatures();
        final CommandLineParser clp = new CommandLineParser(tool);
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
        final TestGATKToolWithFeatures tool = new TestGATKToolWithFeatures();
        final CommandLineParser clp = new CommandLineParser(tool);
        final File vcfFile = new File(publicTestDir + "org/broadinstitute/hellbender/engine/lexicographically_sorted_dict.vcf");
        final String[] args = {"--mask", vcfFile.getCanonicalPath(),
                               "--reference", b37_reference_20_21};
        clp.parseArguments(System.out, args);

        // This method would throw if sequence dictionary validation failed. Here we are testing
        // that it does not throw despite the lexicographically-sorted sequence dictionary in the vcf.
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
        final CommandLineParser clp = new CommandLineParser(tool);
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

}
