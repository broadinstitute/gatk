package org.broadinstitute.hellbender.tools.tumorheterogeneity;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.ExomeStandardArgumentDefinitions;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.List;

/**
 * Integration tests for {@link TumorHeterogeneity}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class TumorHeterogeneityIntegrationTest extends CommandLineProgramTest {
    @Test
    public void testTumorHeterogeneity1Clone() {
        final List<File> ACNV_SEGMENT_FILES = Arrays.asList(
//                new File("/home/slee/working/ipython/purity-ploidy/integration-test/1_clone/seed-0/purity-0.0/total_segments-log2cr_sd-0.001-maf_sd-0.001.acnv.seg"),
//                new File("/home/slee/working/ipython/purity-ploidy/integration-test/1_clone/seed-0/purity-0.1/total_segments-log2cr_sd-0.001-maf_sd-0.001.acnv.seg"),
//                new File("/home/slee/working/ipython/purity-ploidy/integration-test/1_clone/seed-0/purity-0.2/total_segments-log2cr_sd-0.001-maf_sd-0.001.acnv.seg"),
//                new File("/home/slee/working/ipython/purity-ploidy/integration-test/1_clone/seed-0/purity-0.3/total_segments-log2cr_sd-0.001-maf_sd-0.001.acnv.seg"),
//                new File("/home/slee/working/ipython/purity-ploidy/integration-test/1_clone/seed-0/purity-0.4/total_segments-log2cr_sd-0.001-maf_sd-0.001.acnv.seg"),
//                new File("/home/slee/working/ipython/purity-ploidy/integration-test/1_clone/seed-0/purity-0.5/total_segments-log2cr_sd-0.001-maf_sd-0.001.acnv.seg"),
//                new File("/home/slee/working/ipython/purity-ploidy/integration-test/1_clone/seed-0/purity-0.6/total_segments-log2cr_sd-0.001-maf_sd-0.001.acnv.seg"),
//                new File("/home/slee/working/ipython/purity-ploidy/integration-test/1_clone/seed-0/purity-0.7/total_segments-log2cr_sd-0.001-maf_sd-0.001.acnv.seg"),
//                new File("/home/slee/working/ipython/purity-ploidy/integration-test/1_clone/seed-0/purity-0.8/total_segments-log2cr_sd-0.001-maf_sd-0.001.acnv.seg"),
//                new File("/home/slee/working/ipython/purity-ploidy/integration-test/1_clone/seed-0/purity-0.9/total_segments-log2cr_sd-0.001-maf_sd-0.001.acnv.seg"),
//                new File("/home/slee/working/ipython/purity-ploidy/integration-test/1_clone/seed-0/purity-1.0/total_segments-log2cr_sd-0.001-maf_sd-0.001.acnv.seg")
//                new File("/home/slee/working/ipython/purity-ploidy/integration-test/1_clone/seed-1/purity-0.0/total_segments-log2cr_sd-0.001-maf_sd-0.001.acnv.seg"),
//                new File("/home/slee/working/ipython/purity-ploidy/integration-test/1_clone/seed-1/purity-0.1/total_segments-log2cr_sd-0.001-maf_sd-0.001.acnv.seg"),
                new File("/home/slee/working/ipython/purity-ploidy/integration-test/1_clone/seed-1/purity-0.2/total_segments-log2cr_sd-0.001-maf_sd-0.001.acnv.seg"),
//                new File("/home/slee/working/ipython/purity-ploidy/integration-test/1_clone/seed-1/purity-0.3/total_segments-log2cr_sd-0.001-maf_sd-0.001.acnv.seg"),
//                new File("/home/slee/working/ipython/purity-ploidy/integration-test/1_clone/seed-1/purity-0.4/total_segments-log2cr_sd-0.001-maf_sd-0.001.acnv.seg"),
//                new File("/home/slee/working/ipython/purity-ploidy/integration-test/1_clone/seed-1/purity-0.5/total_segments-log2cr_sd-0.001-maf_sd-0.001.acnv.seg"),
//                new File("/home/slee/working/ipython/purity-ploidy/integration-test/1_clone/seed-1/purity-0.6/total_segments-log2cr_sd-0.001-maf_sd-0.001.acnv.seg"),
//                new File("/home/slee/working/ipython/purity-ploidy/integration-test/1_clone/seed-1/purity-0.7/total_segments-log2cr_sd-0.001-maf_sd-0.001.acnv.seg"),
//                new File("/home/slee/working/ipython/purity-ploidy/integration-test/1_clone/seed-1/purity-0.8/total_segments-log2cr_sd-0.001-maf_sd-0.001.acnv.seg"),
//                new File("/home/slee/working/ipython/purity-ploidy/integration-test/1_clone/seed-1/purity-0.9/total_segments-log2cr_sd-0.001-maf_sd-0.001.acnv.seg"),
                new File("/home/slee/working/ipython/purity-ploidy/integration-test/1_clone/seed-1/purity-1.0/total_segments-log2cr_sd-0.001-maf_sd-0.001.acnv.seg")
//                new File("/home/slee/working/ipython/purity-ploidy/integration-test/1_clone/seed-1-cr-only/purity-0.0/total_segments-log2cr_sd-0.001-maf_sd-0.001.acnv.seg"),
//                new File("/home/slee/working/ipython/purity-ploidy/integration-test/1_clone/seed-1-cr-only/purity-0.1/total_segments-log2cr_sd-0.001-maf_sd-0.001.acnv.seg"),
//                new File("/home/slee/working/ipython/purity-ploidy/integration-test/1_clone/seed-1-cr-only/purity-0.2/total_segments-log2cr_sd-0.001-maf_sd-0.001.acnv.seg"),
//                new File("/home/slee/working/ipython/purity-ploidy/integration-test/1_clone/seed-1-cr-only/purity-0.3/total_segments-log2cr_sd-0.001-maf_sd-0.001.acnv.seg"),
//                new File("/home/slee/working/ipython/purity-ploidy/integration-test/1_clone/seed-1-cr-only/purity-0.4/total_segments-log2cr_sd-0.001-maf_sd-0.001.acnv.seg"),
//                new File("/home/slee/working/ipython/purity-ploidy/integration-test/1_clone/seed-1-cr-only/purity-0.5/total_segments-log2cr_sd-0.001-maf_sd-0.001.acnv.seg"),
//                new File("/home/slee/working/ipython/purity-ploidy/integration-test/1_clone/seed-1-cr-only/purity-0.6/total_segments-log2cr_sd-0.001-maf_sd-0.001.acnv.seg"),
//                new File("/home/slee/working/ipython/purity-ploidy/integration-test/1_clone/seed-1-cr-only/purity-0.7/total_segments-log2cr_sd-0.001-maf_sd-0.001.acnv.seg"),
//                new File("/home/slee/working/ipython/purity-ploidy/integration-test/1_clone/seed-1-cr-only/purity-0.8/total_segments-log2cr_sd-0.001-maf_sd-0.001.acnv.seg"),
//                new File("/home/slee/working/ipython/purity-ploidy/integration-test/1_clone/seed-1-cr-only/purity-0.9/total_segments-log2cr_sd-0.001-maf_sd-0.001.acnv.seg"),
//                new File("/home/slee/working/ipython/purity-ploidy/integration-test/1_clone/seed-1-cr-only/purity-1.0/total_segments-log2cr_sd-0.001-maf_sd-0.001.acnv.seg")
        );

        for (final File ACNV_SEGMENT_FILE : ACNV_SEGMENT_FILES) {
            final String OUTPUT_PREFIX = ACNV_SEGMENT_FILE.getAbsolutePath().replace(".acnv.seg", "");
            final String[] arguments = {
                    "--" + ExomeStandardArgumentDefinitions.SEGMENT_FILE_LONG_NAME, ACNV_SEGMENT_FILE.getAbsolutePath(),
                    "--" + TumorHeterogeneity.OUTPUT_PREFIX_LONG_NAME, OUTPUT_PREFIX,
                    "--" + TumorHeterogeneity.MAX_ALLELIC_COPY_NUMBER_CLONAL_LONG_NAME, "4",
                    "--" + TumorHeterogeneity.MAX_ALLELIC_COPY_NUMBER_LONG_NAME, "4",
                    "--" + TumorHeterogeneity.NUM_WALKERS_CLONAL_LONG_NAME, "100",
                    "--" + TumorHeterogeneity.NUM_SAMPLES_CLONAL_LONG_NAME, "500",
                    "--" + TumorHeterogeneity.NUM_BURN_IN_CLONAL_LONG_NAME, "250",
                    "--" + TumorHeterogeneity.NUM_SAMPLES_LONG_NAME, "0",
                    "--" + TumorHeterogeneity.NUM_BURN_IN_LONG_NAME, "0",
                    "--" + TumorHeterogeneity.COPY_RATIO_NORMALIZATION_PRIOR_ALPHA_LONG_NAME, "1E1",
                    "--" + TumorHeterogeneity.COPY_RATIO_NORMALIZATION_PRIOR_BETA_LONG_NAME, "1E1",
                    "--" + TumorHeterogeneity.COPY_RATIO_NOISE_CONSTANT_PRIOR_ALPHA_LONG_NAME, "1",
                    "--" + TumorHeterogeneity.COPY_RATIO_NOISE_CONSTANT_PRIOR_BETA_LONG_NAME, "1E3",
                    "--" + TumorHeterogeneity.PLOIDY_MISMATCH_PENALTY_LONG_NAME, "1E6",
                    "--" + TumorHeterogeneity.MODE_PURITY_BIN_SIZE_LONG_NAME, "0.025",
                    "--" + TumorHeterogeneity.MODE_PLOIDY_BIN_SIZE_LONG_NAME, "0.025",
                    "--verbosity", "INFO"
            };
            runCommandLine(arguments);
        }
    }

    @Test
    public void testTumorHeterogeneityWGS() {
        final List<File> ACNV_SEGMENT_FILES = Arrays.asList(
//                new File("/home/slee/working/ipython/purity-ploidy/integration-test/wgs/TCGA-05-4396-01A-21D-1855-08-sim-final.seg"),
//                new File("/home/slee/working/ipython/purity-ploidy/integration-test/wgs/TCGA-05-4432-01A-01D-1931-08-sim-final.seg"),
//                new File("/home/slee/working/ipython/purity-ploidy/integration-test/wgs/TCGA-38-4628-01A-01D-1931-08-sim-final.seg"),
                new File("/home/slee/working/ipython/purity-ploidy/integration-test/wgs/TCGA-44-2659-01A-01D-1931-08-sim-final.seg"),
                new File("/home/slee/working/ipython/purity-ploidy/integration-test/wgs/TCGA-78-7156-01A-11D-2036-08-sim-final.seg"),
                new File("/home/slee/working/ipython/purity-ploidy/integration-test/wgs/TCGA-78-7158-01A-11D-2036-08-sim-final.seg")
        );

        for (final File ACNV_SEGMENT_FILE : ACNV_SEGMENT_FILES) {
            final String OUTPUT_PREFIX = ACNV_SEGMENT_FILE.getAbsolutePath().replace(".acnv.seg", "");
            final String[] arguments = {
                    "--" + ExomeStandardArgumentDefinitions.SEGMENT_FILE_LONG_NAME, ACNV_SEGMENT_FILE.getAbsolutePath(),
                    "--" + TumorHeterogeneity.OUTPUT_PREFIX_LONG_NAME, OUTPUT_PREFIX,
                    "--" + TumorHeterogeneity.MAX_ALLELIC_COPY_NUMBER_CLONAL_LONG_NAME, "6",
                    "--" + TumorHeterogeneity.MAX_ALLELIC_COPY_NUMBER_LONG_NAME, "6",
                    "--" + TumorHeterogeneity.NUM_WALKERS_CLONAL_LONG_NAME, "100",
                    "--" + TumorHeterogeneity.NUM_SAMPLES_CLONAL_LONG_NAME, "500",
                    "--" + TumorHeterogeneity.NUM_BURN_IN_CLONAL_LONG_NAME, "250",
                    "--" + TumorHeterogeneity.NUM_SAMPLES_LONG_NAME, "0",
                    "--" + TumorHeterogeneity.NUM_BURN_IN_LONG_NAME, "0",
                    "--" + TumorHeterogeneity.COPY_RATIO_NORMALIZATION_PRIOR_ALPHA_LONG_NAME, "1E1",
                    "--" + TumorHeterogeneity.COPY_RATIO_NORMALIZATION_PRIOR_BETA_LONG_NAME, "1E1",
                    "--" + TumorHeterogeneity.COPY_RATIO_NOISE_CONSTANT_PRIOR_ALPHA_LONG_NAME, "1",
                    "--" + TumorHeterogeneity.COPY_RATIO_NOISE_CONSTANT_PRIOR_BETA_LONG_NAME, "1E2",
                    "--" + TumorHeterogeneity.PLOIDY_MISMATCH_PENALTY_LONG_NAME, "1E2",
                    "--" + TumorHeterogeneity.MODE_PURITY_BIN_SIZE_LONG_NAME, "0.025",
                    "--" + TumorHeterogeneity.MODE_PLOIDY_BIN_SIZE_LONG_NAME, "0.025",
                    "--verbosity", "INFO"
            };
            runCommandLine(arguments);
        }
    }

    @Test
    public void testTumorHeterogeneityWES() {
        final List<File> ACNV_SEGMENT_FILES = Arrays.asList(
//                new File("/home/slee/working/ipython/purity-ploidy/integration-test/wes/TCGA-05-4396-01A-21D-1855-08-sim-final.seg"),
//                new File("/home/slee/working/ipython/purity-ploidy/integration-test/wes/TCGA-05-4432-01A-01D-1265-08-sim-final.seg"),
//                new File("/home/slee/working/ipython/purity-ploidy/integration-test/wes/TCGA-38-4628-01A-01D-1265-08-sim-final.seg"),
                new File("/home/slee/working/ipython/purity-ploidy/integration-test/wes/TCGA-44-2659-01A-01D-0969-08-sim-final.seg"),
                new File("/home/slee/working/ipython/purity-ploidy/integration-test/wes/TCGA-78-7156-01A-11D-2036-08-sim-final.seg"),
                new File("/home/slee/working/ipython/purity-ploidy/integration-test/wes/TCGA-78-7158-01A-11D-2036-08-sim-final.seg")
        );

        for (final File ACNV_SEGMENT_FILE : ACNV_SEGMENT_FILES) {
            final String OUTPUT_PREFIX = ACNV_SEGMENT_FILE.getAbsolutePath().replace(".acnv.seg", "");
            final String[] arguments = {
                    "--" + ExomeStandardArgumentDefinitions.SEGMENT_FILE_LONG_NAME, ACNV_SEGMENT_FILE.getAbsolutePath(),
                    "--" + TumorHeterogeneity.OUTPUT_PREFIX_LONG_NAME, OUTPUT_PREFIX,
                    "--" + TumorHeterogeneity.MAX_ALLELIC_COPY_NUMBER_CLONAL_LONG_NAME, "6",
                    "--" + TumorHeterogeneity.MAX_ALLELIC_COPY_NUMBER_LONG_NAME, "6",
                    "--" + TumorHeterogeneity.MAX_NUM_POPULATIONS_LONG_NAME, "3",
                    "--" + TumorHeterogeneity.NUM_WALKERS_CLONAL_LONG_NAME, "100",
                    "--" + TumorHeterogeneity.NUM_SAMPLES_CLONAL_LONG_NAME, "500",
                    "--" + TumorHeterogeneity.NUM_BURN_IN_CLONAL_LONG_NAME, "250",
                    "--" + TumorHeterogeneity.NUM_WALKERS_LONG_NAME, "100",
                    "--" + TumorHeterogeneity.NUM_SAMPLES_LONG_NAME, "200",
                    "--" + TumorHeterogeneity.NUM_BURN_IN_LONG_NAME, "100",
                    "--" + TumorHeterogeneity.COPY_RATIO_NORMALIZATION_PRIOR_ALPHA_LONG_NAME, "1E1",
                    "--" + TumorHeterogeneity.COPY_RATIO_NORMALIZATION_PRIOR_BETA_LONG_NAME, "1E1",
                    "--" + TumorHeterogeneity.COPY_RATIO_NOISE_CONSTANT_PRIOR_ALPHA_LONG_NAME, "1",
                    "--" + TumorHeterogeneity.COPY_RATIO_NOISE_CONSTANT_PRIOR_BETA_LONG_NAME, "1E2",
                    "--" + TumorHeterogeneity.PLOIDY_MISMATCH_PENALTY_LONG_NAME, "1E3",
                    "--" + TumorHeterogeneity.SUBCLONE_VARIANCE_PENALTY_LONG_NAME, "1E3",
                    "--" + TumorHeterogeneity.MODE_PURITY_BIN_SIZE_LONG_NAME, "0.025",
                    "--" + TumorHeterogeneity.MODE_PLOIDY_BIN_SIZE_LONG_NAME, "0.025",
                    "--verbosity", "INFO"
            };
            runCommandLine(arguments);
        }
    }

    @Test
    public void testTumorHeterogeneityHCC() {
        final List<File> ACNV_SEGMENT_FILES = Arrays.asList(
                new File("/home/slee/working/ipython/purity-ploidy/purity-series/0-0-SM-74NEG-sim-final.seg"),
                new File("/home/slee/working/ipython/purity-ploidy/purity-series/1-10-SM-74P2T-sim-final.seg"),
                new File("/home/slee/working/ipython/purity-ploidy/purity-series/2-30-SM-74P35-sim-final.seg"),
                new File("/home/slee/working/ipython/purity-ploidy/purity-series/3-40-SM-74P3J-sim-final.seg"),
                new File("/home/slee/working/ipython/purity-ploidy/purity-series/4-60-SM-74P3M-sim-final.seg"),
                new File("/home/slee/working/ipython/purity-ploidy/purity-series/5-70-SM-74P3K-sim-final.seg"),
                new File("/home/slee/working/ipython/purity-ploidy/purity-series/6-80-SM-74P51-sim-final.seg"),
                new File("/home/slee/working/ipython/purity-ploidy/purity-series/7-90-SM-74P56-sim-final.seg"),
                new File("/home/slee/working/ipython/purity-ploidy/purity-series/8-100-SM-74P4M-sim-final.seg")
        );

        for (final File ACNV_SEGMENT_FILE : ACNV_SEGMENT_FILES) {
            final String OUTPUT_PREFIX = ACNV_SEGMENT_FILE.getAbsolutePath().replace(".acnv.seg", "");
            final String[] arguments = {
                    "--" + ExomeStandardArgumentDefinitions.SEGMENT_FILE_LONG_NAME, ACNV_SEGMENT_FILE.getAbsolutePath(),
                    "--" + TumorHeterogeneity.OUTPUT_PREFIX_LONG_NAME, OUTPUT_PREFIX,
                    "--" + TumorHeterogeneity.MAX_ALLELIC_COPY_NUMBER_CLONAL_LONG_NAME, "6",
                    "--" + TumorHeterogeneity.MAX_ALLELIC_COPY_NUMBER_LONG_NAME, "6",
                    "--" + TumorHeterogeneity.MAX_NUM_POPULATIONS_LONG_NAME, "3",
                    "--" + TumorHeterogeneity.NUM_WALKERS_CLONAL_LONG_NAME, "50",
                    "--" + TumorHeterogeneity.NUM_SAMPLES_CLONAL_LONG_NAME, "200",
                    "--" + TumorHeterogeneity.NUM_BURN_IN_CLONAL_LONG_NAME, "100",
                    "--" + TumorHeterogeneity.NUM_WALKERS_LONG_NAME, "500",
                    "--" + TumorHeterogeneity.NUM_SAMPLES_LONG_NAME, "0",
                    "--" + TumorHeterogeneity.NUM_BURN_IN_LONG_NAME, "0",
                    "--" + TumorHeterogeneity.CONCENTRATION_PRIOR_ALPHA_LONG_NAME, "1",
                    "--" + TumorHeterogeneity.CONCENTRATION_PRIOR_BETA_LONG_NAME, "1E1",
                    "--" + TumorHeterogeneity.COPY_RATIO_NORMALIZATION_PRIOR_ALPHA_LONG_NAME, "1",
                    "--" + TumorHeterogeneity.COPY_RATIO_NORMALIZATION_PRIOR_BETA_LONG_NAME, "1",
                    "--" + TumorHeterogeneity.COPY_RATIO_NOISE_CONSTANT_PRIOR_ALPHA_LONG_NAME, "1",
                    "--" + TumorHeterogeneity.COPY_RATIO_NOISE_CONSTANT_PRIOR_BETA_LONG_NAME, "1E2",
                    "--" + TumorHeterogeneity.PLOIDY_MISMATCH_PENALTY_LONG_NAME, "1E6",
                    "--" + TumorHeterogeneity.SUBCLONE_VARIANCE_PENALTY_LONG_NAME, "1E-6",
                    "--" + TumorHeterogeneity.MODE_PURITY_BIN_SIZE_LONG_NAME, "0.025",
                    "--" + TumorHeterogeneity.MODE_PLOIDY_BIN_SIZE_LONG_NAME, "0.025",
                    "--verbosity", "INFO"
            };
            runCommandLine(arguments);
        }
    }

    @Test
    public void testTumorHeterogeneity2Clone() {
        final List<File> ACNV_SEGMENT_FILES = Arrays.asList(
                new File("/home/slee/working/ipython/purity-ploidy/integration-test/2_clone/seed-3/purity-0.2/total_segments-log2cr_sd-0.001-maf_sd-0.001.acnv.seg"),
                new File("/home/slee/working/ipython/purity-ploidy/integration-test/2_clone/seed-3/purity-0.4/total_segments-log2cr_sd-0.001-maf_sd-0.001.acnv.seg"),
                new File("/home/slee/working/ipython/purity-ploidy/integration-test/2_clone/seed-3/purity-0.6/total_segments-log2cr_sd-0.001-maf_sd-0.001.acnv.seg"),
                new File("/home/slee/working/ipython/purity-ploidy/integration-test/2_clone/seed-3/purity-0.8/total_segments-log2cr_sd-0.001-maf_sd-0.001.acnv.seg"),
                new File("/home/slee/working/ipython/purity-ploidy/integration-test/2_clone/seed-3/purity-1.0/total_segments-log2cr_sd-0.001-maf_sd-0.001.acnv.seg")
        );

        for (final File ACNV_SEGMENT_FILE : ACNV_SEGMENT_FILES) {
            final String OUTPUT_PREFIX = ACNV_SEGMENT_FILE.getAbsolutePath().replace(".acnv.seg", "");
            final String[] arguments = {
                    "--" + ExomeStandardArgumentDefinitions.SEGMENT_FILE_LONG_NAME, ACNV_SEGMENT_FILE.getAbsolutePath(),
                    "--" + TumorHeterogeneity.OUTPUT_PREFIX_LONG_NAME, OUTPUT_PREFIX,
                    "--" + TumorHeterogeneity.MAX_ALLELIC_COPY_NUMBER_CLONAL_LONG_NAME, "4",
                    "--" + TumorHeterogeneity.MAX_ALLELIC_COPY_NUMBER_LONG_NAME, "4",
                    "--" + TumorHeterogeneity.NUM_WALKERS_CLONAL_LONG_NAME, "50",
                    "--" + TumorHeterogeneity.NUM_SAMPLES_CLONAL_LONG_NAME, "100",
                    "--" + TumorHeterogeneity.NUM_BURN_IN_CLONAL_LONG_NAME, "50",
                    "--" + TumorHeterogeneity.NUM_WALKERS_LONG_NAME, "500",
                    "--" + TumorHeterogeneity.NUM_SAMPLES_LONG_NAME, "100",
                    "--" + TumorHeterogeneity.NUM_BURN_IN_LONG_NAME, "50",
                    "--" + TumorHeterogeneity.CONCENTRATION_PRIOR_ALPHA_LONG_NAME, "1",
                    "--" + TumorHeterogeneity.CONCENTRATION_PRIOR_BETA_LONG_NAME, "1E1",
                    "--" + TumorHeterogeneity.COPY_RATIO_NORMALIZATION_PRIOR_ALPHA_LONG_NAME, "1",
                    "--" + TumorHeterogeneity.COPY_RATIO_NORMALIZATION_PRIOR_BETA_LONG_NAME, "1",
                    "--" + TumorHeterogeneity.COPY_RATIO_NOISE_CONSTANT_PRIOR_ALPHA_LONG_NAME, "1",
                    "--" + TumorHeterogeneity.COPY_RATIO_NOISE_CONSTANT_PRIOR_BETA_LONG_NAME, "1E2",
                    "--" + TumorHeterogeneity.PLOIDY_MISMATCH_PENALTY_LONG_NAME, "1E6",
                    "--" + TumorHeterogeneity.SUBCLONE_VARIANCE_PENALTY_LONG_NAME, "1E-6",
                    "--" + TumorHeterogeneity.MODE_PURITY_BIN_SIZE_LONG_NAME, "0.025",
                    "--" + TumorHeterogeneity.MODE_PLOIDY_BIN_SIZE_LONG_NAME, "0.025",
                    "--verbosity", "DEBUG"
            };
            runCommandLine(arguments);
        }
    }
}