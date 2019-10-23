package org.broadinstitute.hellbender.tools.walkers.longreads;

import htsjdk.samtools.*;
import htsjdk.samtools.util.CloseableIterator;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.text.XReadLines;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class ZMWAwareBamSharderIntegrationTest extends CommandLineProgramTest {

    public static final String SMALL_TEST_BAM = "src/test/resources/chrM_and_chr20_subset.remaining_sorted.bam";
    public static final String CUSTOM_GATK_CONFIG = "GATKConfig_bamindexer.properties";

    @Test
    public void testSmallBam() {
        Assert.assertFalse(Defaults.USE_ASYNC_IO_WRITE_FOR_SAMTOOLS);

        final File outputBam = createTempFile("ZMWAwareBamSharderIntegrationTest_testSmallBam", ".bam");
        final File outputIndex = createTempFile("ZMWAwareBamSharderIntegrationTest_testSmallBam", ".index");
        final int TARGET_SHARD_SIZE = 100;

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addArgument("gatk-config-file", CUSTOM_GATK_CONFIG);
        args.addInput(new File(SMALL_TEST_BAM));
        args.addOutput(outputBam);
        args.addArgument("output-index", outputIndex.getAbsolutePath());
        args.addArgument("target-shard-size", Integer.toString(TARGET_SHARD_SIZE));

        runCommandLine(args);

        try ( final XReadLines indexReader = new XReadLines(outputIndex) ) {
            final List<String> indexContents = indexReader.readLines();

            for ( final String indexLine : indexContents ) {
                final String[] tokens = indexLine.split("\\t", -1);
                Assert.assertEquals(tokens.length, 3);

                final String readName = tokens[0];
                final long chunkStart = Long.parseLong(tokens[1]);
                final long chunkEnd = Long.parseLong(tokens[2]);

                final Chunk chunk = new Chunk(chunkStart, chunkEnd);
                final BAMFileSpan span = new BAMFileSpan(chunk);

                final SamReader.PrimitiveSamReaderToSamReaderAdapter bamReader = (SamReader.PrimitiveSamReaderToSamReaderAdapter)SamReaderFactory.makeDefault().open(outputBam);
                final CloseableIterator<SAMRecord> iterator = bamReader.iterator(span);

                final List<SAMRecord> queryResults = new ArrayList<>();
                while ( iterator.hasNext() ) {
                    queryResults.add(iterator.next());
                }

                System.out.println("Num reads in shard = " + queryResults.size());
                Assert.assertEquals(queryResults.get(0).getReadName(), readName);

                iterator.close();
                bamReader.close();
            }

        } catch ( IOException e ) {
            throw new UserException("oops", e);
        }
    }
}

