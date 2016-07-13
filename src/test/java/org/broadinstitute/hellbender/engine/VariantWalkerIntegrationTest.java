package org.broadinstitute.hellbender.engine;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.examples.ExampleVariantWalker;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
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

}
