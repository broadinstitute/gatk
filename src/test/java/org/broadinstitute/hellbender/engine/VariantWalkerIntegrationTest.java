package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineArgumentParser;
import org.broadinstitute.barclay.argparser.CommandLineParser;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.ReadFilterArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.TestProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.examples.ExampleVariantWalker;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

public final class VariantWalkerIntegrationTest extends CommandLineProgramTest {

    @Override
    public String getTestedClassName() {
        return ExampleVariantWalker.class.getSimpleName();
    }

    @Test(expectedExceptions = UserException.class)
    public void testRequiresIndexForInterval() throws Exception {
        String fileIn = "count_variants_withSequenceDict_noIndex.vcf";
        String moreArgs = "-L 1";
        final File ORIG_FILE = new File(getTestDataDir(), fileIn);
        ArgumentsBuilder ab = new ArgumentsBuilder();
        ab.add("--variant " + ORIG_FILE.getAbsolutePath());
        ab.add(moreArgs);
        this.runCommandLine(ab.getArgsArray());
    }

    @Test(expectedExceptions = UserException.MalformedGenomeLoc.class)
    public void testMissingContigForInterval() throws Exception {
        String fileIn = "count_variants_withSequenceDict.vcf";
        String moreArgs = "-L 25";
        final File ORIG_FILE = new File(getTestDataDir(), fileIn);
        ArgumentsBuilder ab = new ArgumentsBuilder();
        ab.add("--variant " + ORIG_FILE.getAbsolutePath());
        ab.add(moreArgs);
        this.runCommandLine(ab.getArgsArray());
    }

    @Test(expectedExceptions = UserException.class)
    public void testRequiresSequenceDictionaryForInterval() throws Exception {
        String fileIn = "count_variants.vcf";
        String moreArgs = "-L 1";
        final File ORIG_FILE = new File(getTestDataDir(), fileIn);
        ArgumentsBuilder ab = new ArgumentsBuilder();
        ab.add("--variant " + ORIG_FILE.getAbsolutePath());
        ab.add(moreArgs);
        this.runCommandLine(ab.getArgsArray());
    }

    @CommandLineProgramProperties(
            summary = "TestGATKToolWithFeatures",
            oneLineSummary = "TestGATKToolWithFeatures",
            programGroup = TestProgramGroup.class
    )
    private static final class TestGATKToolWithFeatures extends VariantWalker{

        public static final String HAS_BACKING_READ_SOURCE_LONG_NAME = "has-backing-read-source";
        public static final String BACKING_READS_LONG_NAME = "backing-reads";

        @Argument(fullName=HAS_BACKING_READ_SOURCE_LONG_NAME)
        boolean hasBackingReadSource = false;

        @Argument(fullName=BACKING_READS_LONG_NAME, optional=true)
        List<String> backingReads= new ArrayList<>();

        @Override
        public void apply(
                VariantContext variant,
                ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext ) {
            Assert.assertEquals(readsContext.hasBackingDataSource(), hasBackingReadSource);
            if (hasBackingReadSource) {
                Iterator<GATKRead> it = readsContext.iterator();
                while (it.hasNext()) {
                    Assert.assertTrue(backingReads.contains(it.next().getName()));
                }
            }
        }
    }

    @Test
    public void testBestSequenceDictionary_FromVariantReference() throws Exception {
        final GATKTool tool = new TestGATKToolWithFeatures();
        final CommandLineParser clp = new CommandLineArgumentParser(tool);
        final File vcfFile = new File(publicTestDir + "org/broadinstitute/hellbender/engine/example_variants_noSequenceDict.vcf");
        final String[] args = {"-V", vcfFile.getCanonicalPath(), "-R", hg19MiniReference};
        clp.parseArguments(System.out, args);
        tool.onStartup();
        // make sure we DON'T get the seq dictionary from the VCF index, and instead get the one from
        // the reference when its better
        final SAMSequenceDictionary toolDict = tool.getBestAvailableSequenceDictionary();
        Assert.assertFalse(toolDict.getSequences().stream().allMatch(seqRec -> seqRec.getSequenceLength() == 0));

        SAMSequenceDictionary refDict = new ReferenceFileSource(IOUtils.getPath(hg19MiniReference)).getSequenceDictionary();
        toolDict.assertSameDictionary(refDict);
        refDict.assertSameDictionary(toolDict);
        Assert.assertEquals(toolDict, refDict);
    }

    @Test
    public void testBestSequenceDictionary_FromVariantIndex() throws Exception {
        final GATKTool tool = new TestGATKToolWithFeatures();
        final CommandLineParser clp = new CommandLineArgumentParser(tool);
        final File vcfFile = new File(publicTestDir + "org/broadinstitute/hellbender/engine/example_variants_noSequenceDict.vcf");
        final String[] args = {"--variant", vcfFile.getCanonicalPath()};
        clp.parseArguments(System.out, args);
        tool.onStartup();
        // make sure we get the seq dict from the index when there isn't one in the VCF and there is no reference available
        final SAMSequenceDictionary toolDict = tool.getBestAvailableSequenceDictionary();
        Assert.assertTrue(toolDict.getSequences().stream().allMatch(seqRec -> seqRec.getSequenceLength() == 0));
    }

    @Test
    public void testReadFilterOff() throws Exception {
        final GATKTool tool = new TestGATKToolWithFeatures();
        final File vcfFile = new File(publicTestDir + "org/broadinstitute/hellbender/engine/VariantWalkerTest_VariantsWithReads.vcf");
        final File bamFile = new File(publicTestDir + "org/broadinstitute/hellbender/engine/VariantWalkerTest_VariantsWithReads.bam");

        final String[] args = {
                "--variant", vcfFile.getCanonicalPath(),
                "--input", bamFile.getCanonicalPath(),
                "--" + TestGATKToolWithFeatures.HAS_BACKING_READ_SOURCE_LONG_NAME, "true",
                // reads we expect to see when using no filter
                "--" + TestGATKToolWithFeatures.BACKING_READS_LONG_NAME, "a",
                "--" + TestGATKToolWithFeatures.BACKING_READS_LONG_NAME, "d"
        };
        tool.instanceMain(args);
    }

    @Test
    public void testReadFilterOn() throws Exception {
        final GATKTool tool = new TestGATKToolWithFeatures();
        final File vcfFile = new File(publicTestDir + "org/broadinstitute/hellbender/engine/VariantWalkerTest_VariantsWithReads.vcf");
        final File bamFile = new File(publicTestDir + "org/broadinstitute/hellbender/engine/VariantWalkerTest_VariantsWithReads.bam");
        final String[] args = {
                "--variant", vcfFile.getCanonicalPath(),
                "--input", bamFile.getCanonicalPath(),
                "--" + TestGATKToolWithFeatures.HAS_BACKING_READ_SOURCE_LONG_NAME, "true",
                "--" + ReadFilterArgumentDefinitions.READ_FILTER_LONG_NAME, "ReadNameReadFilter",
                "--" + ReadFilterArgumentDefinitions.READ_NAME_LONG_NAME, "d",      // only retain reads named "d"
                "--" + TestGATKToolWithFeatures.BACKING_READS_LONG_NAME, "d",  // name of reads we expect to see
        };
        tool.instanceMain(args);
    }

}
