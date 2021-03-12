package org.broadinstitute.hellbender.tools;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

public class CreateReadsBundleJsonIntegrationTest extends CommandLineProgramTest {

//    final String bam = "gs://hellbender/test/resources/benchmark/CEUTrio.HiSeq.WEx.b37.NA12892.bam";
//    final String bai = "gs://hellbender/test/resources/benchmark/CEUTrio.HiSeq.WEx.b37.NA12892.bam.bai";
//    final String cram = "gs://hellbender/test/resources/large/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.cram";
//    final String crai = "gs://hellbender/test/resources/large/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.cram.bai";

    private static final File TEST_LARGE_DATA_DIR = new File("src/test/resources/large/");
//src/test/resources/large/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.bam
//src/test/resources/large/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.bam.bai
//src/test/resources/large/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.cram
//src/test/resources/large/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.cram.bai

    @Test
    //public void testCreateReadsBundle(final String[] testArgs) {
    public void testCreateReadsBundle() throws IOException {
        final Path outputJSON = createTempPath("testCreateReadsBundle", ".json");
        final Path inputBAM = Paths.get("src/test/resources/large/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.bam");
        final Path inputBAI = Paths.get("src/test/resources/large/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.bam.bai");

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addInput(inputBAM);
        args.add(StandardArgumentDefinitions.READ_INDEX_LONG_NAME, inputBAI);
        args.addOutput(outputJSON);

        final Object res = runCommandLine(args);
        Assert.assertNull(res);

        final String outputContents = readEntireFileAsString(new GATKPath(outputJSON.toUri().toString()));
        //{"name":"htsbundle","INDEX":"src/test/resources/large/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.bam.bai","version":"1.0","READS":"src/test/resources/large/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.bam"}
        //Assert.assertEquals(outputContents, "");
        System.out.println(outputContents);
    }

    public String readEntireFileAsString(final GATKPath gatkPath) throws IOException {
        return new String(Files.readAllBytes(gatkPath.toPath()));
    }
}
