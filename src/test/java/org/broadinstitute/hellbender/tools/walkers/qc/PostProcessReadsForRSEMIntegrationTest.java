package org.broadinstitute.hellbender.tools.walkers.qc;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.engine.ReadsPathDataSource;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

public class PostProcessReadsForRSEMIntegrationTest extends CommandLineProgramTest {
    @Test
    public void test(){
        final ArgumentsBuilder args = new ArgumentsBuilder();
        final String testTranscriptomeSam = publicTestDir + "transcriptome_query_sorted_abbr_header.bam";
        final File output = createTempFile("output", "bam");

        args.addInput(testTranscriptomeSam);
        args.addOutput(output.getAbsolutePath());
        args.add("XL", publicTestDir + "gencode_v19_MT_transcriptome_abbr_header.interval_list");
        runCommandLine(args, PostProcessReadsForRSEM.class.getSimpleName());

        final ReadsPathDataSource outputReadSource = new ReadsPathDataSource(output.toPath());
        final Iterator<GATKRead> outputSamIterator = outputReadSource.iterator();

        final List<String> inputReadNames = Arrays.asList("SL-NVZ:HNJ2HDRXY211220:HNJ2HDRXY:2:2101:10004:11804",
                "SL-NVZ:HNJ2HDRXY211220:HNJ2HDRXY:2:2101:10004:27367",
                "SL-NVZ:HNJ2HDRXY211220:HNJ2HDRXY:2:2101:10004:30906",
                "SL-NVZ:HNJ2HDRXY211220:HNJ2HDRXY:2:2101:10004:36761",
                "SL-NVZ:HNJ2HDRXY211220:HNJ2HDRXY:2:2101:10004:9079",
                "SL-NVZ:HNJ2HDRXY211220:HNJ2HDRXY:2:2101:10276:1689"); // aligns to a transcript in MT
        final int[] expectedReadPairCount = new int[] { 4, 5, 2, 1, 7, 0 };
        String currentReadName = inputReadNames.get(0);
        int i = 0;
        int j = 0;
        int count = 0;

        while (outputSamIterator.hasNext()){
            final GATKRead outputRead = outputSamIterator.next();
            // Check that read1 and read2 alternate
            Assert.assertTrue(i % 2 == 0 ? outputRead.isFirstOfPair() : ! outputRead.isFirstOfPair());
            if (! outputRead.getName().equals(currentReadName)){
                Assert.assertEquals(count, 2*expectedReadPairCount[j++]);
                currentReadName = outputRead.getName();
                count = 1;
                i++;
                continue;
            }

            count++;
            i++;
        }

    }
}