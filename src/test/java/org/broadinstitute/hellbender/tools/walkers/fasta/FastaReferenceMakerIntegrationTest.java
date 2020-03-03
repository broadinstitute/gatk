package org.broadinstitute.hellbender.tools.walkers.fasta;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.BaseTest;
import org.broadinstitute.hellbender.testutils.FastaTestUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.List;

public class FastaReferenceMakerIntegrationTest extends CommandLineProgramTest {

    @DataProvider
    public Object[][] getFastaParameters(){
        return new Object[][]{
                {Arrays.asList("1:10,000,100-10,000,500", "1:10,100,000-10,101,000", "1:10,900,000-10,900,001"), getTestFile("reference_only.fasta")},
                {Arrays.asList("1:10,000,100-10,000,200", "1:10,000,201-10,000,301"), getTestFile("reference_only_contiguous_same_contig.fasta")},
                {Arrays.asList("1:10,000,100-10,000,200", "2:10,000,201-10,000,301"), getTestFile("reference_only_contiguous_diff_contigs.fasta")}
        };
    }

    @Test(dataProvider = "getFastaParameters")
    public void runMakeFastaTest(List<String> intervals, File expected) throws IOException {
        ArgumentsBuilder args = new ArgumentsBuilder();
        final File out = BaseTest.createTempFile("subset", ".fasta");
        args.addReference(new File(b37Reference))
                .addOutput(out);
        intervals.forEach(interval -> args.add("L", interval));
        runCommandLine(args);

        FastaTestUtils.assertFastaFilesContainTheSameSequence(out.toPath(), expected.toPath());
    }


    @DataProvider
    public Object[][] getBasesPerLine(){
        return new Object[][]{
                {13}, {60}, {100}
        };
    }

    @Test(dataProvider = "getBasesPerLine")
    public void testBasePerLine(int basesPerLine) throws IOException {
        ArgumentsBuilder args = new ArgumentsBuilder();
        final File out = BaseTest.createTempFile("subset", ".fasta");
        final File hg19mini = new File(hg19MiniReference);
        args.addReference(hg19mini)
                .addOutput(out)
                .add(FastaReferenceMaker.LINE_WIDTH_LONG_NAME, String.valueOf(basesPerLine));

        runCommandLine(args);

        FastaTestUtils.assertFastaFilesContainTheSameSequence(out.toPath(), hg19mini.toPath());
        final List<String> lines = Files.readAllLines(out.toPath());
        for(int i = 0; i < lines.size(); i++){
            String line = lines.get(i);
            if(!line.startsWith(">")){
                Assert.assertTrue(line.length() == basesPerLine || i == lines.size()-1 || lines.get(i+1).startsWith(">"));
            }
        }

    }
}