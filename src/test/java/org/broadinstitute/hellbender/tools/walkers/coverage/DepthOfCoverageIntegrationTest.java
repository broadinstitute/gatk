package org.broadinstitute.hellbender.tools.walkers.coverage;

import com.opencsv.CSVReader;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

public class DepthOfCoverageIntegrationTest extends CommandLineProgramTest {

    private String getDoCExtensionFromFile(File docOutputFile, String basename) {
        String[] split = docOutputFile.getAbsolutePath().split(basename);
        if (split.length == 2) {
            return split[1];
        } else if (split.length == 1) {
            return "";
        }
        Assert.fail("There was a problem loading DoCExtension: "+basename);
        return null;
    }

    // Files were taken from GATK3 outputs being run with equivalent arguments
    private File getExpectedDataDir() {
        return getTestFile( "expected/");
    }

    @Test
    // Baseline test asserting that histogram window boundaries work and that coverage counts are correct when all base quality scores are included.
    public void testBaseOutputNoFiltering() throws IOException {
        final String expectedBaseName = "depthofcoveragenofiltering";
        final File baseOutputFile = createTempDir("depthofcoveragenofiltering");
        final File output = IOUtils.createTempFileInDirectory( "depthofcoveragenofiltering", ".csv", baseOutputFile);

        String cmd = "-R "+hg38Reference+" " +
                "-I "+largeFileTestDir + "multiSampleSubsetted.bam " +
                "-L "+ getTestFile("artificial.target_region.interval_list ")+
                "--min-base-quality 0 --include-deletions --print-base-counts -pt readgroup -pt sample -pt platform -pt library --output-format CSV --summary-coverage-threshold 10 --summary-coverage-threshold 15 --summary-coverage-threshold 20 --summary-coverage-threshold 25";
        cmd += " -O "+output.getAbsolutePath();
        runCommandLine(cmd.split(" "));

        File[] actualFiles = baseOutputFile.listFiles();

        compareOutputDirectories(expectedBaseName, output.getName(), actualFiles);
    }

    @Test
    // Asserting that the per-base statistics are reasonable when run without an interval list file
    public void testCoverageBehaviorWhenProvidedNoIntervalFile() throws IOException {
        final String expectedBaseName = "testCoverageBehaviorWhenProvidedNoIntervalFile";
        final File baseOutputFile = createTempDir("testCoverageBehaviorWhenProvidedNoIntervalFile");
        final File output = IOUtils.createTempFileInDirectory( "testCoverageBehaviorWhenProvidedNoIntervalFile", ".csv", baseOutputFile);

        String cmd = "-R "+hg38Reference+" " +
        "-I "+largeFileTestDir+"multiSampleSubsetted.bam"+
        " -O "+output.getAbsolutePath() + " -L chr1:1570000-1571000 --omit-interval-statistics --omit-locus-table --omit-per-sample-statistics";
        runCommandLine(cmd.split(" "));
        File[] actualFiles = baseOutputFile.listFiles();

        compareOutputDirectories(expectedBaseName, output.getName(), actualFiles);
    }

    @Test
    // Coverage drops to nothing because '--max-base-quality 4' is set, resulting in very poor output coverage
    public void testNoCoverageDueToFiltering() throws IOException {
        final String expectedBaseName = "depthofcoveragenwithiltering";
        final File baseOutputFile = createTempDir("testNoCoverageDueToFiltering");
        final File output = IOUtils.createTempFileInDirectory( "testNoCoverageDueToFiltering", ".csv", baseOutputFile);

        String cmd = "-R "+hg38Reference+" " +
                "-I "+largeFileTestDir + "multiSampleSubsetted.bam " +
                "-L "+ getTestFile("artificial.target_region.interval_list ")+
                "--min-base-quality 5 --max-base-quality 4 --include-deletions --print-base-counts -pt readgroup -pt sample -pt library -pt platform --output-format CSV";
        cmd += " -O "+output.getAbsolutePath();
        runCommandLine(cmd.split(" "));

        File[] actualFiles = baseOutputFile.listFiles();

        compareOutputDirectories(expectedBaseName, output.getName(), actualFiles);
    }

    @Test
    // NOTE, the gene list file was not generated with GATK3 due to different gene list merging behavior
    public void testGeneListAllGenesCompletelyCoveredByIntervals() throws IOException {
        final String expectedBaseName = "testGeneListDataAllCovered";
        final File baseOutputFile = createTempDir("testGeneListAllGenesCompletelyCoveredByIntervals");
        final File output = IOUtils.createTempFileInDirectory( "testGeneListAllGenesCompletelyCoveredByIntervals", ".csv", baseOutputFile);
        output.delete();

        String[] cmd = new String[]{ "-R",hg38Reference,
                "-I",largeFileTestDir+"multiSampleSubsetted.bam",
                "-L",getTestFile("artificial.gene_target.interval_list").getAbsolutePath(),
                "--calculate-coverage-over-genes",getTestFile("refGene_CDK11B.refseq").getAbsolutePath(),
                "--min-base-quality","0","--include-deletions","-pt","sample","--output-format","CSV","--omit-depth-output-at-each-base","--omit-per-sample-statistics",
                "-O",output.getAbsolutePath()};
        runCommandLine(cmd);

        File[] actualFiles = baseOutputFile.listFiles();

        compareOutputDirectories(expectedBaseName, output.getName(), actualFiles);
    }

    @Test
    // This is asserting that the new locus overlapping behavior works without merging
    public void testIntervalListOverlappingUniqueBehavior() throws IOException {
        final String expectedBaseName = "testIntervalListOverlappingUniqueBehavior"; //interval_summary files were not generated with gatk3 as they behave by interval merging
        final File baseOutputFile = createTempDir("testNoCoverageDueToFiltering");
        final File output = IOUtils.createTempFileInDirectory( "testNoCoverageDueToFiltering", ".csv", baseOutputFile);
        output.delete();

        ArgumentsBuilder cmd = new ArgumentsBuilder();
        cmd.addReference(new File(hg38Reference));
        cmd.addInput(new File(largeFileTestDir+"multiSampleSubsetted.bam"));
        cmd.add("min-base-quality","0");
        cmd.addRaw("--include-deletions");
        cmd.add("partition-type", "sample");
        cmd.addRaw("--omit-depth-output-at-each-base");
        cmd.addOutput(output);

        // We construct overlapping intervals to demonstrate the new interval merging behavior
        for (int i = 1656275; i < 1656975; i+=10) {
            cmd.addInterval(new SimpleInterval("chr1", i, i+100));
        }
        runCommandLine(cmd);
        File[] actualFiles = baseOutputFile.listFiles();

        compareOutputDirectories(expectedBaseName, output.getName(), actualFiles);
    }


    @Test
    // NOTE, the gene list file was not generated with GATK3 due to different gene list merging behavior
    public void testGeneListDataPartallylCovered() throws IOException {
        final File baseOutputFile = createTempDir("testNoCoverageDueToFiltering");
        final File output = IOUtils.createTempFileInDirectory( "testNoCoverageDueToFiltering", ".csv", baseOutputFile);

        String[] cmd = new String[]{ "-R",hg38Reference,
                "-I",largeFileTestDir+"multiSampleSubsetted.bam",
                "-L", "chr1:1656275-1677440", // this should only completely span two of the genes in the list
                "--calculate-coverage-over-genes",getTestFile("refGene_CDK11B.refseq").getAbsolutePath(),
                "--min-base-quality","0","--include-deletions","-pt","sample","--output-format","CSV","--omit-depth-output-at-each-base","--omit-per-sample-statistics","--omit-genes-not-entirely-covered-by-traversal",
                "-O",output.getAbsolutePath()};
        runCommandLine(cmd);

        try (final CSVReader actualParser = new CSVReader(new FileReader(output.getAbsolutePath()+".sample_gene_summary"))) {
            actualParser.readNext(); // Skip header line
            List<String[]> lines = actualParser.readAll();
            Assert.assertEquals(lines.size(), 2);
            // The only two genes completely encompassed by this interval span
            Assert.assertEquals(lines.get(0)[0], "SLC35E2");
            Assert.assertEquals(lines.get(1)[0], "SLC35E2");
            Assert.assertNotEquals(lines.get(0), lines.get(1));
        }
    }

    @Test
    // NOTE, the gene list file was not generated with GATK3 due to different gene list merging behavior
    public void testGeneListDataPartallylCoveredIntronGap() throws IOException {
        final File baseOutputFile = createTempDir("testNoCoverageDueToFiltering");
        final File output = IOUtils.createTempFileInDirectory( "testNoCoverageDueToFiltering", ".csv", baseOutputFile);

        String[] cmd = new String[]{ "-R",hg38Reference,
                "-I",largeFileTestDir+"multiSampleSubsetted.bam",
                "-L", "chr1:1656275-1666275", "-L", "chr1:1666280-1677440", // this has a gap that aligns with introns for both genes
                "--calculate-coverage-over-genes",getTestFile("refGene_CDK11B.refseq").getAbsolutePath(),
                "--min-base-quality","0","--include-deletions","-pt","sample","--output-format","CSV","--omit-depth-output-at-each-base","--omit-per-sample-statistics","--omit-genes-not-entirely-covered-by-traversal",
                "-O",output.getAbsolutePath()};
        runCommandLine(cmd);

        try (final CSVReader actualParser = new CSVReader(new FileReader(output.getAbsolutePath()+".sample_gene_summary"))) {
            actualParser.readNext(); // Skip header line
            List<String[]> lines = actualParser.readAll();
            Assert.assertEquals(lines.size(), 2);
            // The only two genes completely encompassed by this interval span
            Assert.assertEquals(lines.get(0)[0], "SLC35E2");
            Assert.assertEquals(lines.get(1)[0], "SLC35E2");
            Assert.assertNotEquals(lines.get(0), lines.get(1));
        }
    }

    private void compareOutputDirectories(final String expectedBaseName, final String actualFileBaseName, final File[] actual) throws IOException {
        List<File> expectedFiles = Arrays.stream(Objects.requireNonNull(getExpectedDataDir().listFiles())).filter(f -> f.getName().contains(expectedBaseName)).sorted().collect(Collectors.toList());
        List<File> actualFiles = Arrays.stream(actual).sorted().collect(Collectors.toList());

        // First check that all the expected files correspond to actual files
        List<String> expectedFilenames = expectedFiles.stream().map(f -> getDoCExtensionFromFile(f, expectedBaseName)).sorted().collect(Collectors.toList());
        List<String> actualFilenames = actualFiles.stream().map(f -> getDoCExtensionFromFile(f, actualFileBaseName)).sorted().collect(Collectors.toList());
        Assert.assertEquals(actualFilenames, expectedFilenames);

        // Now assert that the outputs exactly match with the expected outputs
        for (int i = 0; i < actual.length; i++) {
            // Need this line for comparison to gatk3 in which the interval_summary file has non-deterministic header/sample line ordering
            compareCSVByRowHeaderOrderAgnostic(actualFiles.get(i), expectedFiles.get(i), ',', actualFileBaseName);
        }
    }

    // Helper class that allows us to get around mismatches sample line ordering that may mismatch in results files
    private void compareCSVByRowHeaderOrderAgnostic( final File actualFile, final File expectedFile, char separator, String actualFileBaseName) {
        try (final CSVReader actualParser = new CSVReader(new FileReader(actualFile), separator);
             final CSVReader expectedParser = new CSVReader(new FileReader(expectedFile), separator);
        ) {
            List<String> actualHeader = Arrays.asList(actualParser.readNext());
            List<String> expectedHeader = Arrays.asList(expectedParser.readNext());
            // Deal with archaic standard of delimiting header with a different character from the body
            if (expectedHeader.size() == 1) {
                expectedHeader = new ArrayList<>(Arrays.asList(expectedHeader.get(0).split("[\t,]")));
            }
            //Assert headers match paying no heed to ordering
            Assert.assertEquals(actualHeader.stream().sorted().iterator(), expectedHeader.stream().sorted().iterator());


            Iterator<String[]> actualLines = StreamSupport.stream(actualParser.spliterator(),false).sorted(Comparator.comparing(l -> l[0]+","+l[1])).iterator();
            Iterator<String[]> expectedLines = StreamSupport.stream(expectedParser.spliterator(),false).sorted(Comparator.comparing(l -> l[0]+","+l[1])).iterator();
            // Special case to handle scrambled output header for the base coverage in gatk3
            if (getDoCExtensionFromFile(actualFile, actualFileBaseName).equals("")) {
                // Just assert all the correct counts are present for each line because the order is wrong in the expected files
                int linecount = 0;
                while (actualLines.hasNext() && expectedLines.hasNext()) {
                    List<String> actualLine = Arrays.asList(actualLines.next());
                    List<String> expectedLine = Arrays.asList(expectedLines.next());
                    Assert.assertEquals(actualLine.stream().sorted().iterator(), expectedLine.stream().sorted().iterator(), "Base output files disagreed at line "+linecount);
                    linecount++;
                }
                if (actualLines.hasNext()) {
                    Assert.fail("Finished parsing expected file and found lines remaining in the actual file  "+actualFile.getName());
                } else if (expectedLines.hasNext()) {
                    Assert.fail("Reached end of parsing expected file lines and found lines remaining "+expectedFile.getName());
                }

            } else {
                compareLinesByColumn(actualFile, expectedFile, actualHeader, expectedHeader, actualLines, expectedLines);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }


    private void compareLinesByColumn(File actualFile, File expectedFile, List<String> actualHeader, List<String> expectedHeader, Iterator<String[]> actualLines, Iterator<String[]> expectedLines) {
        while (actualLines.hasNext() && expectedLines.hasNext()) {
            // Assert that the text value is the same for each header line
            String[] actual = actualLines.next();
            String[] expected = expectedLines.next();
            for (int a = 0, e = expectedHeader.indexOf(actualHeader.get(0));
                 a < actualHeader.size();
                 a++, e = a < actualHeader.size() ? expectedHeader.indexOf(actualHeader.get(a)) : 0) {
                // Treat the 'Total' line specially as gatk3 didn't assign the right number of columns to that row
                if (actual[0].equals("Total")) {
                    // gatk3 produced incomplete lines for the total line before, assert that the data fields are valid.
                    Assert.assertEquals(Arrays.copyOf(actual, 3) , Arrays.copyOf(expected, 3));
                } else {
                    Assert.assertEquals(actual[a], expected[e], "Actual and expected fields differed at header line : "+ actualHeader.get(a) +", for files " + actualFile + " and " + expectedFile);
                }
            }
        }
        if (actualLines.hasNext()) {
            Assert.fail("Finished parsing expected file and found lines remaining in the actual file  "+actualFile.getName());
        } else if (expectedLines.hasNext()) {
            Assert.fail("Reached end of parsing expected file lines and found lines remaining "+expectedFile.getName());
        }
    }


}