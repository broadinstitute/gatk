package org.broadinstitute.hellbender.tools.gvs.ingest;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import static htsjdk.samtools.ValidationStringency.STRICT;

public class CreateVariantIngestFilesIntegrationTest extends CommandLineProgramTest {

    // *manual_nocalls.vcf.gz modified manually to introduce a no-call (./.) GT at chr20 61417
    // no-call to test marking ./. as GQ 0 in PET and abstaining from putting variant into VET
    // original GT at chr20 61417 = 0/1
    final String input_vcf_file = "NA12878.haplotypeCalls.reblocked.chr20.100k.manual_nocall.vcf.gz";
    final String interval_list_file = "wgs_calling_regions.hg38.chr20.100k.interval_list";
    final String sample_map_file = "test_sample_map.tsv";

    // TSV column order for VET, matching VetFieldEnum.values() and the existing golden TSV header
    private static final List<String> VET_TSV_COLUMNS =
            Arrays.stream(VetFieldEnum.values()).map(Enum::name).collect(Collectors.toList());

    // TSV column order for ref_ranges, matching RefRangesTsvWriter's hard-coded header order
    private static final List<String> REF_RANGES_TSV_COLUMNS =
            List.of("location", "sample_id", "length", "state");

    @Test
    public void testCreateDefaultVariantIngestFiles() throws Exception {
        testCreateVariantIngestFiles(
                getToolTestDataDir() + "expected.ref_ranges_001_NA12878.tsv",
                getToolTestDataDir() + "expected.vet_001_NA12878.tsv",
                new String[]{});
    }

    @Test
    public void testIgnoreGQ0CreateVariantIngestFiles() throws Exception {
        testCreateVariantIngestFiles(
                getToolTestDataDir() + "expected.ignore_gq0.ref_ranges_001_NA12878.tsv",
                getToolTestDataDir() + "expected.vet_001_NA12878.tsv",
                new String[]{"--ref-block-gq-to-ignore ZERO"});
    }

    @Test
    public void testIgnoreGQ40AndAboveCreateVariantIngestFiles() throws Exception {
        testCreateVariantIngestFiles(
                getToolTestDataDir() + "expected.ignore_gq40_and_above.ref_ranges_001_NA12878.tsv",
                getToolTestDataDir() + "expected.vet_001_NA12878.tsv",
                new String[]{"--ref-block-gq-to-ignore FORTY", "--ignore-above-gq-threshold true"});
    }

    public void testCreateVariantIngestFiles(final String expectedRefRangesTsv,
                                             final String expectedVetTsv,
                                             final String[] additionalArgs) throws Exception {

        final File outputDir = createTempDir("output_dir");
        final File convertedDir = createTempDir("converted_");

        final ArgumentsBuilder args = new ArgumentsBuilder(additionalArgs);
        args
                .add("V", getToolTestDataDir() + input_vcf_file)
                .add("L", getToolTestDataDir() + interval_list_file)
                .add("interval-set-rule", "INTERSECTION")
                .add("output-type", "PARQUET")
                .add("SNM", getToolTestDataDir() + sample_map_file)
                .add("enable-reference-ranges", true)
                .add("ref-version", "38")
                .add("output-directory", outputDir)
        ;
        runCommandLine(args);

        // Find the generated Parquet files by their filename prefix
        final File vetParquetFile = findParquetFileByPrefix(outputDir, "vet_");
        final File refRangesParquetFile = findParquetFileByPrefix(outputDir, "ref_ranges_");

        // PARQUET mode also emits a sample_chromosome_ploidy file; verify it was created
        Assert.assertNotNull(findParquetFileByPrefix(outputDir, "sample_chromosome_ploidy_"),
                "Expected a sample_chromosome_ploidy Parquet file to be created");

        // Convert each Parquet file to TSV so we can compare against the existing golden TSV files
        final File convertedVetTsv = new File(convertedDir, "vet.tsv");
        final File convertedRefRangesTsv = new File(convertedDir, "ref_ranges.tsv");
        ParquetToTsvConverter.convert(vetParquetFile, VET_TSV_COLUMNS, convertedVetTsv);
        ParquetToTsvConverter.convert(refRangesParquetFile, REF_RANGES_TSV_COLUMNS, convertedRefRangesTsv);

        // Compare converted TSV contents against the golden TSV files.
        // List order: ref_ranges first, then vet — matching the golden-file list order.
        final List<File> convertedFiles = Arrays.asList(convertedRefRangesTsv, convertedVetTsv);
        final List<String> expectedOutputFiles = Arrays.asList(expectedRefRangesTsv, expectedVetTsv);
        IntegrationTestSpec.assertMatchingFiles(convertedFiles, expectedOutputFiles, false, STRICT);
    }

    private File findParquetFileByPrefix(final File dir, final String prefix) {
        final File[] matches = dir.listFiles((d, name) -> name.startsWith(prefix) && name.endsWith(".parquet"));
        Assert.assertNotNull(matches, "Could not list files in output directory");
        Assert.assertEquals(matches.length, 1,
                "Expected exactly one Parquet file with prefix '" + prefix + "' in " + dir);
        return matches[0];
    }
}
