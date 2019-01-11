package org.broadinstitute.hellbender.tools.walkers.bqsr;

import htsjdk.samtools.ValidationStringency;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.Main;
import org.broadinstitute.hellbender.tools.validation.CompareBaseQualities;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;

//Note: has to be in this package to have access to Main.instanceMain
public final class BothStepsOfBQSRIntegrationTest extends CommandLineProgramTest {

    @Test //Tests that BQSR with and without indel recalibration produces the same bam file
    public void testIndelRemoval() throws Exception {
        //Note; using a small 27M file for speed
        final File bamIn = new File(NA12878_20_21_WGS_bam);
        final String interval = "20";

        final File recalOutWithIndels = baseRecalibrator(bamIn, interval, false);
        final File bamOutWithIndels = applyBQSR(bamIn, interval, recalOutWithIndels, false);

        final File recalOutWithoutIndels = baseRecalibrator(bamIn, interval, true);
        final File bamOutWithoutIndels = applyBQSR(bamIn, interval, recalOutWithoutIndels, true);

        final ArgumentsBuilder args1 = new ArgumentsBuilder();
        args1.addArgument("VALIDATION_STRINGENCY", ValidationStringency.SILENT.name());
        args1.addPositionalArgument(bamOutWithIndels.getAbsolutePath());
        args1.addPositionalArgument(bamOutWithoutIndels.getAbsolutePath());
        final Object result = new Main().instanceMain(makeCommandLineArgs(args1.getArgsList(), CompareBaseQualities.class.getSimpleName()));
        Assert.assertEquals(result, 0);
    }

    private File applyBQSR(final File bamIn, final String interval, final File recalOut, final boolean skipIndels) {
        final File bamOut = GATKBaseTest.createTempFile("applyBQSR." + skipIndels, ".bam");
        final ArgumentsBuilder args1 = new ArgumentsBuilder();
        args1.addInput(bamIn);
        args1.addFileArgument("bqsr", recalOut);
        args1.addArgument("L", interval);
        args1.addOutput(bamOut);
        new Main().instanceMain(makeCommandLineArgs(args1.getArgsList(), ApplyBQSR.class.getSimpleName()));
        return bamOut;
    }

    private File baseRecalibrator(final File bamIn, final String interval, final boolean skipIndels) {
        final File recalOut = GATKBaseTest.createTempFile("baseRecalibrator." + skipIndels, ".recal");

        final ArgumentsBuilder args1 = new ArgumentsBuilder();
        args1.addInput(bamIn);
        args1.addOutput(recalOut);
        args1.addArgument("L", interval);
        args1.addFileArgument("known-sites", new File(dbsnp_138_b37_20_21_vcf));
        args1.addReference(new File(b37_reference_20_21));
        args1.addBooleanArgument("indels", !skipIndels);
        new Main().instanceMain(makeCommandLineArgs(args1.getArgsList(), BaseRecalibrator.class.getSimpleName()));
        return recalOut;
    }


}
