package org.broadinstitute.hellbender;

import com.intel.gkl.compression.IntelDeflater;
import com.intel.gkl.compression.IntelDeflaterFactory;
import com.intel.gkl.compression.IntelInflater;
import com.intel.gkl.compression.IntelInflaterFactory;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.BlockGunzipper;
import htsjdk.samtools.util.zip.DeflaterFactory;
import htsjdk.samtools.util.zip.InflaterFactory;
import org.broadinstitute.hellbender.utils.NativeUtils;
import org.broadinstitute.hellbender.utils.RandomDNA;
import org.broadinstitute.hellbender.testutils.SamAssertionUtils;
import org.testng.Assert;
import org.testng.SkipException;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.zip.DataFormatException;
import java.util.zip.Deflater;
import java.util.zip.Inflater;

/**
 * Test GATK integration with the Intel inflater and deflater
 */
public class IntelInflaterDeflaterIntegrationTest extends CommandLineProgramTest {

    private static final String INPUT_FILE = "CEUTrio.HiSeq.WGS.b37.NA12878.20.21.tiny.md.bam";

    @Override
    public String getTestedToolName() {
        return "PrintReads";
    }

    private boolean isIntelInflaterDeflaterSupported() {
        return (NativeUtils.runningOnLinux() || NativeUtils.runningOnMac()) && !NativeUtils.runningOnPPCArchitecture();
    }

    @Test
    public void testIntelInflaterIsAvailable() {
        if ( ! NativeUtils.runningOnLinux()  && ! NativeUtils.runningOnMac()) {
            throw new SkipException("IntelInflater is not available on this platform");
        }

        if ( NativeUtils.runningOnPPCArchitecture() ) {
            throw new SkipException("IntelInflater is not available for this architecture");
        }

        Assert.assertTrue(new IntelInflater().load(null), "IntelInflater shared library was not loaded. " +
                "This could be due to a configuration error, or your system might not support it.");
    }

    @Test
    public void testIntelDeflaterIsAvailable() {
        if ( ! NativeUtils.runningOnLinux()  && ! NativeUtils.runningOnMac()) {
            throw new SkipException("IntelDeflater is not available on this platform");
        }

        if ( NativeUtils.runningOnPPCArchitecture() ) {
            throw new SkipException("IntelDeflater is not available for this architecture");
        }

        Assert.assertTrue(new IntelDeflater().load(null), "IntelDeflater shared library was not loaded. " +
                "This could be due to a configuration error, or your system might not support it.");
    }

    @DataProvider(name = "JdkFlags")
    public Object[][] makeJdkFlags() {
        final List<Object[]> tests = new ArrayList<>();

        // flags used by PrintReads: {--use_jdk_inflater, --use_jdk_deflater}
        tests.add(new Object[]{false, true});
        tests.add(new Object[]{true, false});
        tests.add(new Object[]{false, false});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "JdkFlags")
    public void testIntelInflaterDeflaterWithPrintReads(final boolean use_jdk_inflater, final boolean use_jdk_deflater) throws Exception {
        if (!isIntelInflaterDeflaterSupported()) {
            throw new SkipException("IntelInflater/IntelDeflater not available on this platform");
        }

        final File ORIG_BAM = new File(largeFileTestDir, INPUT_FILE);
        final File outFile = GATKBaseTest.createTempFile(INPUT_FILE, ".bam");

        final ArrayList<String> args = new ArrayList<>();
        args.add("--input"); args.add(ORIG_BAM.getAbsolutePath());
        args.add("--output"); args.add(outFile.getAbsolutePath());
        args.add("--use-jdk-inflater"); args.add(String.valueOf(use_jdk_inflater));
        args.add("--use-jdk-deflater"); args.add(String.valueOf(use_jdk_deflater));

        // store current default factories, so they can be restored later
        InflaterFactory currentInflaterFactory = BlockGunzipper.getDefaultInflaterFactory();
        DeflaterFactory currentDeflaterFactory = BlockCompressedOutputStream.getDefaultDeflaterFactory();

        // set default factories to jdk version
        // because PrintReads cannot change the factory to Jdk if it was already set to Intel
        BlockGunzipper.setDefaultInflaterFactory(new InflaterFactory());
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());

        // run PrintReads
        runCommandLine(args);

        // restore default factories
        BlockGunzipper.setDefaultInflaterFactory(currentInflaterFactory);
        BlockCompressedOutputStream.setDefaultDeflaterFactory(currentDeflaterFactory);

        // validate input and output files are the same
        SamAssertionUtils.assertSamsEqual(outFile, ORIG_BAM);
    }

    @Test
    public void deflateInflateWithIntel() throws DataFormatException {
        if (!isIntelInflaterDeflaterSupported()) {
            throw new SkipException("IntelInflater/IntelDeflater not available on this platform");
        }

        // create buffers and random input
        final int LEN = 64 * 1024;
        final byte[] input = new RandomDNA().nextBases(LEN);
        final byte[] compressed = new byte[2 * LEN];
        final byte[] result = new byte[LEN];

        final IntelInflaterFactory intelInflaterFactory = new IntelInflaterFactory();
        final IntelDeflaterFactory intelDeflaterFactory = new IntelDeflaterFactory();
        Assert.assertTrue(intelInflaterFactory.usingIntelInflater());
        Assert.assertTrue(intelDeflaterFactory.usingIntelDeflater());

        for (int i = 0; i < 10; i++) {
            // create deflater with compression level i
            final Deflater deflater = intelDeflaterFactory.makeDeflater(i, true);

            // setup deflater
            deflater.setInput(input);
            deflater.finish();

            // compress data
            int compressedBytes = 0;

            // buffer for compressed data is 2x the size of the data we are compressing
            // so this loop should always finish in one iteration
            while (!deflater.finished()) {
                compressedBytes = deflater.deflate(compressed, 0, compressed.length);
            }
            deflater.end();

            // decompress and check output == input
            Inflater inflater = intelInflaterFactory.makeInflater(true);

            inflater.setInput(compressed, 0, compressedBytes);
            inflater.inflate(result);
            inflater.end();

            Assert.assertEquals(input, result);

            // clear compressed and result buffers for next iteration
            Arrays.fill(compressed, (byte)0);
            Arrays.fill(result, (byte)0);
        }
    }
}
