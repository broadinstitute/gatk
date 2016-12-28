package org.broadinstitute.hellbender.tools.walkers.mutect;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.testng.Assert;
import org.testng.annotations.Test;
import java.io.File;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Stream;

import static org.testng.Assert.*;

/**
 * Created by tsato on 12/28/16.
 */
public class CountFalsePositivesTest extends CommandLineProgramTest {

    /***
     * Make sure that the tool runs to completion and prints out the correct number of lines
     */
    @Test
    public void testEvaluateMutect2() throws Exception {
        final String dreamDir =  publicTestDir + "org/broadinstitute/hellbender/tools/mutect/dream/vcfs/";
        List<String> vcfList = Arrays.asList(dreamDir + "sample_1.vcf", dreamDir + "sample_2.vcf", dreamDir + "sample_3.vcf", dreamDir + "sample_4.vcf");
        final File vcfListFile = createTempFile("vcfList", ".txt");
        final File output = createTempFile("output", ".txt");
        Path file = Paths.get(vcfListFile.getAbsolutePath());
        Files.write(file, vcfList, Charset.forName("UTF-8"));
        final String[] args = {
                "-V", vcfListFile.toString(),
                "-O", output.toString()
        };
        runCommandLine(args);
        Stream<String> outputLines = Files.lines(Paths.get(output.getAbsolutePath()));
        Assert.assertEquals(outputLines.count(), 9L);

    }

}