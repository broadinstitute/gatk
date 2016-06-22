package org.broadinstitute.hellbender;

import com.intel.gkl.compression.IntelDeflater;
import com.intel.gkl.compression.IntelDeflaterFactory;
import htsjdk.samtools.*;
import htsjdk.samtools.util.zip.DeflaterFactory;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.NativeUtils;
import org.broadinstitute.hellbender.utils.RandomDNA;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.SkipException;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;
import java.util.zip.Deflater;
import java.util.zip.Inflater;

/**
 * Test that it's possible to load libIntelDeflater
 */
public class IntelDeflaterIntegrationTest extends BaseTest {

    private final static Logger log = LogManager.getLogger(IntelDeflaterIntegrationTest.class);
    private static final String INPUT_FILE = publicTestDir + "CEUTrio.HiSeq.WGS.b37.NA12878.20.21.bam";


    private boolean isIntelDeflaterSupported() {
        return (NativeUtils.runningOnLinux() || NativeUtils.runningOnMac()) && !NativeUtils.runningOnPPCArchitecture();
    }

    @Test
    public void testIntelDeflaterIsAvailable(){
        if ( ! NativeUtils.runningOnLinux()  && ! NativeUtils.runningOnMac()) {
            throw new SkipException("IntelDeflater not available on this platform");
        }

        if ( NativeUtils.runningOnPPCArchitecture() ) {
            throw new SkipException("IntelDeflater not available for this architecture");
        }

        Assert.assertTrue(new IntelDeflater().load(), "libIntelDeflater.so was not loaded. " +
                "This could be due to a configuration error, or your system might not support it.");
    }

    @Test
    public void deflateInflateWithIntelDeflater() {
        if (!isIntelDeflaterSupported()) {
            throw new SkipException("IntelDeflater not available on this platform");
        }

        // create buffers and random input
        final int LEN = 64 * 1024;
        final byte[] input = new RandomDNA().nextBases(LEN);
        final byte[] compressed = new byte[2 * LEN];
        final byte[] result = new byte[LEN];

        final IntelDeflaterFactory intelDeflaterFactory = new IntelDeflaterFactory();

        for (int i = 0; i < 10; i++) {
            // create deflater with compression level i
            final Deflater deflater = intelDeflaterFactory.makeDeflater(i, true);
            Assert.assertTrue(intelDeflaterFactory.usingIntelDeflater());

            // setup deflater
            deflater.reset();
            deflater.setInput(input);
            deflater.finish();

            // compress data
            int compressedBytes = 0;
            while (!deflater.finished()) {
                compressedBytes = deflater.deflate(compressed, 0, compressed.length);
            }
            deflater.end();

            // log results
            log.info("%d bytes compressed to %d bytes : %2.2f%% compression\n",
                     LEN, compressedBytes, 100.0 - 100.0 * compressedBytes / LEN);

            // decompress and check output == input
            Inflater inflater = new Inflater(true);
            try {
                inflater.setInput(compressed, 0, compressedBytes);
                inflater.inflate(result);
                inflater.end();
            } catch (java.util.zip.DataFormatException e) {
                e.printStackTrace();
            }

            Assert.assertEquals(input, result);
        }
    }

    @Test(enabled = false)
    public void bamWriteReadWithIntelDeflater() throws IOException {
        if (!isIntelDeflaterSupported()) {
            throw new SkipException("IntelDeflater not available on this platform");
        }

        // create input and output files
        final File inputFile = new File(INPUT_FILE);
        final File outputFile = File.createTempFile("output", ".bam");
        outputFile.deleteOnExit();

        // setup SAM reader
        SamReaderFactory readerFactory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
        readerFactory = readerFactory.enable(SamReaderFactory.Option.EAGERLY_DECODE);

        // create deflater factory for IntelDeflater
        final DeflaterFactory intelDeflaterFactory = new IntelDeflaterFactory();

        log.info("input filesize = " + inputFile.length());

        // loop through all compression levels
        for (int compressionLevel = 0; compressionLevel < 10; compressionLevel++) {
            long totalRecords = 0;
            try (final SamReader reader = readerFactory.open(inputFile)) {
                // create input SAM reader and output SAM writer
                final SAMFileHeader header = reader.getFileHeader();
                final SAMFileWriterFactory writerFactory = new SAMFileWriterFactory();
                writerFactory.setCompressionLevel(compressionLevel);
// TODO: remove comment when HTSJDK is updated
//                writerFactory.setDeflaterFactory(intelDeflaterFactory);
                final SAMFileWriter writer = writerFactory.makeBAMWriter(header, true, outputFile);

                // read input and write to output
                long totalTime = 0;
                for (final SAMRecord record : reader) {
                    final long start = System.currentTimeMillis();
                    writer.addAlignment(record);
                    totalTime += System.currentTimeMillis() - start;
                    totalRecords++;
                }

                writer.close();

                log.info(String.format("PROFILE %d %.3f %d",
                        compressionLevel, totalTime / 1000.0, outputFile.length()));
            }

            // check that input reads match output reads
            log.info("Checking generated output. Total records = " + totalRecords);
            final SamReader expectedFile = readerFactory.open(inputFile);
            final SamReader generatedFile = readerFactory.open(outputFile);
            Iterator<SAMRecord> generatedIterator = generatedFile.iterator();

            for (final SAMRecord expected : expectedFile) {
                SAMRecord generated = generatedIterator.next();
                assert(expected.toString().equals(generated.toString()));
            }
        }
    }
}
