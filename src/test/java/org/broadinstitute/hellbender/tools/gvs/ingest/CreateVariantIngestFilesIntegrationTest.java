package org.broadinstitute.hellbender.tools.gvs.ingest;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.tools.gvs.common.IngestConstants;
import org.broadinstitute.hellbender.tools.gvs.common.IngestUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Objects;
import java.util.Set;
import java.util.stream.Collectors;

import static htsjdk.samtools.ValidationStringency.STRICT;

public class CreateVariantIngestFilesIntegrationTest extends CommandLineProgramTest {

    // *manual_nocalls.vcf.gz modified manually to introduce a no-call (./.) GT at chr20 61417
    // no-call to test marking ./. as GQ 0 in PET and abstaining from putting variant into VET
    // original GT at chr20 61417 = 0/1
    private static final String INPUT_VCF_FILE = "NA12878.haplotypeCalls.reblocked.chr20.100k.manual_nocall.vcf.gz";
    final String interval_list_file = "wgs_calling_regions.hg38.chr20.100k.interval_list";
    final String sample_map_file = "test_sample_map.tsv";

    // Sample ID as recorded in test_sample_map.tsv
    private static final long SAMPLE_ID = 2L;
    // Superpartition (table number) derived the same way the production code does
    private static final String TABLE_NUMBER =
            String.format("%03d", IngestUtils.getTableNumber(SAMPLE_ID, IngestConstants.partitionPerTable));

    // Expected Parquet filename prefixes, matching the naming logic in VetCreator / RefRangesCreator /
    // SamplePloidyCreator (note: ploidy filenames omit the table number)
    private static final String VET_PARQUET_PREFIX         = "vet_"                      + TABLE_NUMBER + "_" + SAMPLE_ID + "_";
    private static final String REF_RANGES_PARQUET_PREFIX  = "ref_ranges_"               + TABLE_NUMBER + "_" + SAMPLE_ID + "_";
    private static final String PLOIDY_PARQUET_PREFIX      = "sample_chromosome_ploidy_"               + SAMPLE_ID + "_";

    // Exact expected Parquet output filenames
    private static final String VET_PARQUET_FILENAME        = VET_PARQUET_PREFIX        + INPUT_VCF_FILE + ".parquet";
    private static final String REF_RANGES_PARQUET_FILENAME = REF_RANGES_PARQUET_PREFIX + INPUT_VCF_FILE + ".parquet";
    private static final String PLOIDY_PARQUET_FILENAME     = PLOIDY_PARQUET_PREFIX     + INPUT_VCF_FILE + ".parquet";

    // TSV column order for VET, matching VetFieldEnum.values() and the existing golden TSV header
    private static final List<String> VET_TSV_COLUMNS =
            Arrays.stream(VetFieldEnum.values()).map(Enum::name).collect(Collectors.toList());

    // Column order for converting ref_ranges Parquet to TSV. The golden files use location-first
    // ordering (location, sample_id, length, state), which differs from the Parquet schema's
    // sample_id-first field order; ParquetToTsvConverter uses this list to reorder on read.
    private static final List<String> REF_RANGES_TSV_COLUMNS =
            List.of("location", "sample_id", "length", "state");

    // TSV column order for sample_chromosome_ploidy, matching the Parquet schema field order
    private static final List<String> PLOIDY_TSV_COLUMNS =
            List.of("sample_id", "chromosome", "ploidy");

    // All three test cases use the same input sample on the same chromosome, so the ploidy
    // output is identical regardless of which GQ bands are filtered.
    private static final String EXPECTED_PLOIDY_TSV =
            "expected.sample_chromosome_ploidy_NA12878.tsv";

    @Test
    public void testCreateDefaultVariantIngestFiles() throws Exception {
        testCreateVariantIngestFiles(
                getToolTestDataDir() + "expected.ref_ranges_001_NA12878.tsv",
                getToolTestDataDir() + "expected.vet_001_NA12878.tsv",
                getToolTestDataDir() + EXPECTED_PLOIDY_TSV,
                new String[]{});
    }

    @Test
    public void testIgnoreGQ0CreateVariantIngestFiles() throws Exception {
        testCreateVariantIngestFiles(
                getToolTestDataDir() + "expected.ignore_gq0.ref_ranges_001_NA12878.tsv",
                getToolTestDataDir() + "expected.vet_001_NA12878.tsv",
                getToolTestDataDir() + EXPECTED_PLOIDY_TSV,
                new String[]{"--ref-block-gq-to-ignore ZERO"});
    }

    @Test
    public void testIgnoreGQ40AndAboveCreateVariantIngestFiles() throws Exception {
        testCreateVariantIngestFiles(
                getToolTestDataDir() + "expected.ignore_gq40_and_above.ref_ranges_001_NA12878.tsv",
                getToolTestDataDir() + "expected.vet_001_NA12878.tsv",
                getToolTestDataDir() + EXPECTED_PLOIDY_TSV,
                new String[]{"--ref-block-gq-to-ignore FORTY", "--ignore-above-gq-threshold true"});
    }

    public void testCreateVariantIngestFiles(final String expectedRefRangesTsv,
                                             final String expectedVetTsv,
                                             final String expectedPloidyTsv,
                                             final String[] additionalArgs) throws Exception {

        final File outputDir = createTempDir("output_dir");
        final File convertedDir = createTempDir("converted_");

        final ArgumentsBuilder args = new ArgumentsBuilder(additionalArgs);
        args
                .add("V", getToolTestDataDir() + INPUT_VCF_FILE)
                .add("L", getToolTestDataDir() + interval_list_file)
                .add("interval-set-rule", "INTERSECTION")
                .add("output-type", "PARQUET")
                .add("SNM", getToolTestDataDir() + sample_map_file)
                .add("enable-reference-ranges", true)
                .add("ref-version", "38")
                .add("output-directory", outputDir)
        ;
        runCommandLine(args);

        // Locate the generated Parquet files by their exact expected names
        final File vetParquetFile        = new File(outputDir, VET_PARQUET_FILENAME);
        final File refRangesParquetFile  = new File(outputDir, REF_RANGES_PARQUET_FILENAME);
        final File ploidyParquetFile     = new File(outputDir, PLOIDY_PARQUET_FILENAME);
        Assert.assertTrue(vetParquetFile.exists(),       "Expected VET Parquet file not found: "        + VET_PARQUET_FILENAME);
        Assert.assertTrue(refRangesParquetFile.exists(), "Expected ref_ranges Parquet file not found: " + REF_RANGES_PARQUET_FILENAME);
        Assert.assertTrue(ploidyParquetFile.exists(),    "Expected ploidy Parquet file not found: "     + PLOIDY_PARQUET_FILENAME);

        // Assert that the output directory contains exactly the three expected Parquet files and nothing else,
        // to catch regressions that accidentally emit extra outputs (e.g. duplicate partitions, stray files).
        // CRC sidecar files (.<name>.crc) are an expected Hadoop local-filesystem artifact and are excluded.
        final Set<String> actualOutputFileNames = Arrays.stream(Objects.requireNonNull(outputDir.listFiles()))
                .map(File::getName)
                .filter(name -> !name.endsWith(".crc"))
                .collect(Collectors.toSet());
        final Set<String> expectedOutputFileNames = new HashSet<>(Arrays.asList(
                VET_PARQUET_FILENAME, REF_RANGES_PARQUET_FILENAME, PLOIDY_PARQUET_FILENAME));
        Assert.assertEquals(actualOutputFileNames, expectedOutputFileNames,
                "Output directory should contain exactly the three expected Parquet files");

        // Convert each Parquet file to TSV so we can compare against the existing golden TSV files
        final File convertedVetTsv = new File(convertedDir, "vet.tsv");
        final File convertedRefRangesTsv = new File(convertedDir, "ref_ranges.tsv");
        final File convertedPloidyTsv = new File(convertedDir, "sample_chromosome_ploidy.tsv");
        ParquetToTsvConverter.convert(vetParquetFile, VET_TSV_COLUMNS, convertedVetTsv);
        ParquetToTsvConverter.convert(refRangesParquetFile, REF_RANGES_TSV_COLUMNS, convertedRefRangesTsv);
        ParquetToTsvConverter.convert(ploidyParquetFile, PLOIDY_TSV_COLUMNS, convertedPloidyTsv);

        // Compare converted TSV contents against the golden TSV files.
        // List order: ref_ranges, then vet, then ploidy.
        final List<File> convertedFiles = Arrays.asList(convertedRefRangesTsv, convertedVetTsv, convertedPloidyTsv);
        final List<String> expectedOutputFiles = Arrays.asList(expectedRefRangesTsv, expectedVetTsv, expectedPloidyTsv);
        IntegrationTestSpec.assertMatchingFiles(convertedFiles, expectedOutputFiles, false, STRICT);
    }

    @Test
    public void testCreateVrsAlleleParquetFiles() throws Exception {
        // This test requires a local hg38 reference. Set the system property
        // "gvs.test.hg38.reference" to the path of Homo_sapiens_assembly38.fasta,
        // or place the file at the default location below.
        final String refPath = System.getProperty(
                "gvs.test.hg38.reference",
                "/Users/hatcher/Projects/local_ref/Homo_sapiens_assembly38.fasta");
        final File refFile = new File(refPath);
        if (!refFile.exists()) {
            throw new org.testng.SkipException(
                    "Skipping testCreateVrsAlleleParquetFiles: reference not found at " + refPath
                    + ". Set -Dgvs.test.hg38.reference=<path> to enable.");
        }

        final File outputDir = createTempDir("output_dir_vrs");
        final File convertedDir = createTempDir("converted_vrs_");

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args
                .add("R", refPath)
                .add("V", getToolTestDataDir() + INPUT_VCF_FILE)
                .add("L", getToolTestDataDir() + interval_list_file)
                .add("interval-set-rule", "INTERSECTION")
                .add("output-type", "PARQUET")
                .add("SNM", getToolTestDataDir() + sample_map_file)
                .add("enable-reference-ranges", true)
                .add("enable-vrs-ids", true)
                .add("ref-version", "38")
                .add("output-directory", outputDir)
        ;
        runCommandLine(args);

        final List<File> allOutputFiles = Arrays.stream(Objects.requireNonNull(outputDir.listFiles()))
                .filter(f -> !f.getName().endsWith(".crc"))
                .sorted()
                .collect(Collectors.toList());

        // Verify a vrs_allele TSV file was created alongside the standard outputs
        final List<String> fileNames = allOutputFiles.stream()
                .map(File::getName)
                .collect(Collectors.toList());
        Assert.assertTrue(fileNames.stream().anyMatch(n -> n.startsWith("vrs_allele_")),
                "Expected a vrs_allele_*.parquet file in output, got: " + fileNames);
        Assert.assertTrue(fileNames.stream().anyMatch(n -> n.startsWith("vet_")),
                "Expected a vet_*.parquet file in output");
        Assert.assertTrue(fileNames.stream().anyMatch(n -> n.startsWith("ref_ranges_")),
                "Expected a ref_ranges_*.parquet file in output");

        // Convert the VRS allele Parquet to TSV, then check its header and that it has at least one data row
        final File vrsFile = allOutputFiles.stream()
                .filter(f -> f.getName().startsWith("vrs_allele_"))
                .findFirst()
                .orElseThrow(() -> new AssertionError("vrs_allele file not found"));
        final File vrsTsv = new File(convertedDir, "vrs_allele.tsv");
        ParquetToTsvConverter.convert(vrsFile, List.of("vrs_allele_id", "vrs_location_id", "refget_accession",
                "ref_genome_coord_start", "ref_genome_coord_end", "state_type"), vrsTsv);
        final List<String> vrsLines = Files.readAllLines(vrsTsv.toPath());
        Assert.assertTrue(vrsLines.size() > 1,
                "Expected at least one data row in vrs_allele TSV, got: " + vrsLines.size() + " lines");

        // Verify the header contains all expected column names
        final String header = vrsLines.get(0);
        for (String col : new String[]{"vrs_allele_id", "vrs_location_id", "refget_accession",
                "ref_genome_coord_start", "ref_genome_coord_end", "state_type"}) {
            Assert.assertTrue(header.contains(col), "Header missing column: " + col);
        }

        // Verify each data row has a properly formatted vrs_allele_id
        for (int i = 1; i < vrsLines.size(); i++) {
            final String line = vrsLines.get(i);
            Assert.assertTrue(line.startsWith("ga4gh:VA."),
                    "Row " + i + " vrs_allele_id should start with ga4gh:VA., got: " + line);
        }
    }
}
