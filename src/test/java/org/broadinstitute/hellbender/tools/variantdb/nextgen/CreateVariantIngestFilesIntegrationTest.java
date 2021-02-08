package org.broadinstitute.hellbender.tools.variantdb.nextgen;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.File;
import java.util.*;

import static htsjdk.samtools.ValidationStringency.STRICT;

public class CreateVariantIngestFilesIntegrationTest extends CommandLineProgramTest {

    @Test
    public void testPetTsvIsCreated() throws Exception {
        final String input_vcf_file = "NA12878.haplotypeCalls.reblocked.chr20.100k.vcf.gz";
        final String interval_list_file = "wgs_calling_regions.hg38.chr20.100k.interval_list";
        final String sample_map_file = "test_sample_map.tsv";
        final File outputDir = createTempDir("output_dir");
        final List<String> expectedOutputFiles = new ArrayList<>(Arrays.asList(
                getToolTestDataDir() + "expected.metadata_001_NA12878.tsv",
                getToolTestDataDir() + "expected.pet_001_NA12878.tsv",
                getToolTestDataDir() + "expected.vet_001_NA12878.tsv"));

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args
                .add("V", getToolTestDataDir() + input_vcf_file)
                .add("L", getToolTestDataDir() + interval_list_file)
                .add("interval-set-rule", "INTERSECTION")
                .add("mode", "GENOMES")
                .add("output-type", "TSV")
                .add("ref-block-gq-to-ignore", "SIXTY")
                .add("SNM", getToolTestDataDir() + sample_map_file)
                .add("ref-version", "38")
                .add("output-directory", outputDir)
        ;
        runCommandLine(args);

        final List<File> allOutputFiles = Arrays.asList(
                Objects.requireNonNull(outputDir.listFiles())
        );
        Collections.sort(allOutputFiles);

        IntegrationTestSpec.assertMatchingFiles(allOutputFiles, expectedOutputFiles, false, STRICT);
    }
}
