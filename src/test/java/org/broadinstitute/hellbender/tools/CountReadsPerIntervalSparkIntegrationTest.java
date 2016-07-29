package org.broadinstitute.hellbender.tools;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.testng.annotations.Test;

import java.io.File;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.util.List;

import static org.testng.Assert.assertEquals;
import static org.testng.AssertJUnit.assertFalse;

public final class CountReadsPerIntervalSparkIntegrationTest extends CommandLineProgramTest {

    @Test
    public void test() throws Exception {
        List<String> outNoShuffle = run(false);
        List<String> outWithShuffle = run(true);
        assertFalse(outNoShuffle.isEmpty());
        assertEquals(outNoShuffle, outWithShuffle);
    }

    @Test
    public void testWithShuffle() throws Exception {
        final File bam = new File(publicTestDir + "large/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.bam");
        final File out = File.createTempFile("out", ".txt");
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--input");
        args.add(bam.getAbsolutePath());
        args.add("--output");
        args.add(out.getAbsolutePath());
        args.add("--shuffle");
        this.runCommandLine(args.getArgsArray());

        Files.readAllLines(out.toPath(), Charset.forName("UTF-8"));
    }

    public List<String> run(boolean shuffle) throws Exception {
        final File bam = new File(publicTestDir + "large/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.bam");
        final File out = File.createTempFile("out", ".txt");
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--input");
        args.add(bam.getAbsolutePath());
        args.add("--output");
        args.add(out.getAbsolutePath());
        if (shuffle) {
            args.add("--shuffle");
        }
        this.runCommandLine(args.getArgsArray());
        return Files.readAllLines(out.toPath(), Charset.forName("UTF-8"));
    }


}