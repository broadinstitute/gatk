package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.CommandLineParser;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.TestProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.examples.ExampleVariantWalker;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;

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

        @Override
        public boolean requiresFeatures() {
            return true;
        }

        public void apply(
                VariantContext variant,
                ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext ) {
            // no-op
        }
    }

    @Test
    public void testBestSequenceDictionary_FromVariantReference() throws Exception {
        final GATKTool tool = new TestGATKToolWithFeatures();
        final CommandLineParser clp = new CommandLineParser(tool);
        final File vcfFile = new File(publicTestDir + "org/broadinstitute/hellbender/engine/example_variants_noSequenceDict.vcf");
        final String[] args = {"-V", vcfFile.getCanonicalPath(), "-R", hg19MiniReference};
        clp.parseArguments(System.out, args);
        tool.onStartup();
        // make sure we DON'T get the seq dictionary from the VCF index, and instead get the one from
        // the reference when its better
        final SAMSequenceDictionary toolDict = tool.getBestAvailableSequenceDictionary();
        Assert.assertFalse(toolDict.getSequences().stream().allMatch(seqRec -> seqRec.getSequenceLength() == 0));

        SAMSequenceDictionary refDict = new ReferenceFileSource(new File(hg19MiniReference)).getSequenceDictionary();
        toolDict.assertSameDictionary(refDict);
        refDict.assertSameDictionary(toolDict);
        Assert.assertEquals(toolDict, refDict);
    }

    @Test
    public void testBestSequenceDictionary_FromVariantIndex() throws Exception {
        final GATKTool tool = new TestGATKToolWithFeatures();
        final CommandLineParser clp = new CommandLineParser(tool);
        final File vcfFile = new File(publicTestDir + "org/broadinstitute/hellbender/engine/example_variants_noSequenceDict.vcf");
        final String[] args = {"--variant", vcfFile.getCanonicalPath()};
        clp.parseArguments(System.out, args);
        tool.onStartup();
        // make sure we get the seq dict from the index when there isn't one in the VCF and there is no reference available
        final SAMSequenceDictionary toolDict = tool.getBestAvailableSequenceDictionary();
        Assert.assertTrue(toolDict.getSequences().stream().allMatch(seqRec -> seqRec.getSequenceLength() == 0));
    }

}
