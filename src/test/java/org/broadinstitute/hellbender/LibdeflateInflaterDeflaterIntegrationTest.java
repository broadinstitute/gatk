package org.broadinstitute.hellbender;

import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.BlockGunzipper;
import htsjdk.samtools.util.zip.DeflaterFactory;
import htsjdk.samtools.util.zip.InflaterFactory;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.testutils.SamAssertionUtils;
import org.broadinstitute.hellbender.utils.NativeUtils;
import org.testng.Assert;
import org.testng.SkipException;
import org.testng.annotations.AfterClass;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

/**
 * Tests that GATK tools pick up HTSJDK's libdeflate-backed inflater and deflater by default, and that the
 * --use-jdk-inflater / --use-jdk-deflater arguments switch to the JDK implementations.
 */
public class LibdeflateInflaterDeflaterIntegrationTest extends CommandLineProgramTest {

    private static final String INPUT_FILE = "CEUTrio.HiSeq.WGS.b37.NA12878.20.21.tiny.md.bam";

    @Override
    public String getTestedToolName() {
        return "PrintReads";
    }

    private static boolean isLibdeflateSupported() {
        return (NativeUtils.runningOnLinux() || NativeUtils.runningOnMac()) && !NativeUtils.runningOnPPCArchitecture();
    }

    @Test
    public void libdeflateIsUsedByDefault() throws IOException {
        if (!isLibdeflateSupported()) {
            throw new SkipException("libdeflate is not available on this platform");
        }
        runPrintReads(false, false);
        Assert.assertEquals(CommandLineProgram.getDeflaterName(), "LibdeflateDeflater");
        Assert.assertEquals(CommandLineProgram.getInflaterName(), "LibdeflateInflater");
    }

    @Test
    public void jdkDeflaterFlagSelectsJdkDeflater() throws IOException {
        runPrintReads(false, true);
        Assert.assertEquals(CommandLineProgram.getDeflaterName(), "Deflater");
    }

    @Test
    public void jdkInflaterFlagSelectsJdkInflater() throws IOException {
        runPrintReads(true, false);
        Assert.assertEquals(CommandLineProgram.getInflaterName(), "Inflater");
    }

    @Test
    public void bothJdkFlagsSelectJdkImplementations() throws IOException {
        runPrintReads(true, true);
        Assert.assertEquals(CommandLineProgram.getDeflaterName(), "Deflater");
        Assert.assertEquals(CommandLineProgram.getInflaterName(), "Inflater");
    }

    /** Runs PrintReads with the given flags, asserts the round-tripped BAM matches the input, and leaves the
     *  factories installed by that run in place so the caller can inspect them. */
    private void runPrintReads(final boolean useJdkInflater, final boolean useJdkDeflater) throws IOException {
        final File origBam = new File(largeFileTestDir, INPUT_FILE);
        final File outFile = GATKBaseTest.createTempFile(INPUT_FILE, ".bam");

        final ArrayList<String> args = new ArrayList<>();
        args.add("--input"); args.add(origBam.getAbsolutePath());
        args.add("--output"); args.add(outFile.getAbsolutePath());
        args.add("--use-jdk-inflater"); args.add(String.valueOf(useJdkInflater));
        args.add("--use-jdk-deflater"); args.add(String.valueOf(useJdkDeflater));

        runCommandLine(args);
        SamAssertionUtils.assertSamsEqual(outFile, origBam);
    }

    @AfterClass
    public void restoreDefaultFactories() {
        // Other tests in the JVM share the static factories; put HTSJDK's defaults back.
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());
        BlockGunzipper.setDefaultInflaterFactory(new InflaterFactory());
    }
}
