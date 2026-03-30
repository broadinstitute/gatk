package org.broadinstitute.hellbender.tools.gvs.ingest;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Objects;
import java.util.stream.Collectors;

import static htsjdk.samtools.ValidationStringency.STRICT;

public class CreateVariantIngestFilesIntegrationTest extends CommandLineProgramTest {

    // *manual_nocalls.vcf.gz modified manually to introduce a no-call (./.) GT at chr20 61417
    // no-call to test marking ./. as GQ 0 in PET and abstaining from putting variant into VET
    // original GT at chr20 61417 = 0/1
    final String input_vcf_file = "NA12878.haplotypeCalls.reblocked.chr20.100k.manual_nocall.vcf.gz";
    final String interval_list_file = "wgs_calling_regions.hg38.chr20.100k.interval_list";
    final String sample_map_file = "test_sample_map.tsv";
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
        final List<String> expectedOutputFiles = Arrays.asList(
                expectedRefRangesTsv,
                expectedVetTsv
        );

        final ArgumentsBuilder args = new ArgumentsBuilder(additionalArgs);
        args
                .add("V", getToolTestDataDir() + input_vcf_file)
                .add("L", getToolTestDataDir() + interval_list_file)
                .add("interval-set-rule", "INTERSECTION")
                .add("output-type", "TSV")
                .add("SNM", getToolTestDataDir() + sample_map_file)
                .add("enable-reference-ranges", true)
                .add("ref-version", "38")
                .add("output-directory", outputDir)
        ;
        runCommandLine(args);

        final List<File> allOutputFiles = Arrays.asList(
                Objects.requireNonNull(outputDir.listFiles())
        );
        Collections.sort(allOutputFiles);

        final List<File> expectedFileNames = Arrays.asList(
                new File(outputDir + "/ref_ranges_001_" + input_vcf_file + ".tsv"),
                new File(outputDir + "/vet_001_" + input_vcf_file + ".tsv")
        );

        // test contents of files match
        IntegrationTestSpec.assertMatchingFiles(allOutputFiles, expectedOutputFiles, false, STRICT);

        // test file names are expected
        Assert.assertEquals(allOutputFiles, expectedFileNames);
    }

    @Test
    public void testCreateVrsAlleleTsvFiles() throws Exception {
        // This test requires a local hg38 reference. Set the system property
        // "gvs.test.hg38.reference" to the path of Homo_sapiens_assembly38.fasta,
        // or place the file at the default location below.
        final String refPath = System.getProperty(
                "gvs.test.hg38.reference",
                "/Users/hatcher/Projects/local_ref/Homo_sapiens_assembly38.fasta");
        final File refFile = new File(refPath);
        if (!refFile.exists()) {
            throw new org.testng.SkipException(
                    "Skipping testCreateVrsAlleleTsvFiles: reference not found at " + refPath
                    + ". Set -Dgvs.test.hg38.reference=<path> to enable.");
        }

        final File outputDir = createTempDir("output_dir_vrs");

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args
                .add("R", refPath)
                .add("V", getToolTestDataDir() + input_vcf_file)
                .add("L", getToolTestDataDir() + interval_list_file)
                .add("interval-set-rule", "INTERSECTION")
                .add("output-type", "TSV")
                .add("SNM", getToolTestDataDir() + sample_map_file)
                .add("enable-reference-ranges", true)
                .add("enable-vrs-ids", true)
                .add("ref-version", "38")
                .add("output-directory", outputDir)
        ;
        runCommandLine(args);

        final List<File> allOutputFiles = Arrays.asList(
                Objects.requireNonNull(outputDir.listFiles())
        );
        Collections.sort(allOutputFiles);

        // Verify a vrs_allele TSV file was created alongside the standard outputs
        final List<String> fileNames = allOutputFiles.stream()
                .map(File::getName)
                .collect(Collectors.toList());
        Assert.assertTrue(fileNames.stream().anyMatch(n -> n.startsWith("vrs_allele_")),
                "Expected a vrs_allele_*.tsv file in output, got: " + fileNames);
        Assert.assertTrue(fileNames.stream().anyMatch(n -> n.startsWith("vet_")),
                "Expected a vet_*.tsv file in output");
        Assert.assertTrue(fileNames.stream().anyMatch(n -> n.startsWith("ref_ranges_")),
                "Expected a ref_ranges_*.tsv file in output");

        // Verify the VRS allele TSV has the expected header and at least one data row
        final File vrsFile = allOutputFiles.stream()
                .filter(f -> f.getName().startsWith("vrs_allele_"))
                .findFirst()
                .orElseThrow(() -> new AssertionError("vrs_allele file not found"));
        final List<String> vrsLines = Files.readAllLines(vrsFile.toPath());
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
