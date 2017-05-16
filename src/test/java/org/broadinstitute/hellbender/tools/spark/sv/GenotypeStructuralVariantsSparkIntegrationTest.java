package org.broadinstitute.hellbender.tools.spark.sv;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.IntervalArgumentCollection;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by valentin on 4/21/17.
 */
public class GenotypeStructuralVariantsSparkIntegrationTest extends CommandLineProgramTest {

    private static final long serialVersionUID = 1L;

    private static final File INPUT_VCF = new File(publicTestDir, "/test-files/sv_evidence_sam_for_variants_input_20_21.vcf.gz");
    private static final File INPUT_BAM = new File(publicTestDir, "/test-files/sv_evidence_for_variants_input.bam");
    private static final File INPUT_FASTQ_DIR = new File(publicTestDir, "/test-files/sv_evidence_for_variants_fastqs");
    private static final File INPUT_INTERVALS = new File(publicTestDir, "/test-files/sv_evidence_for_variants_intervals.list");

    @Test
    public void test() {
        final File outputVcf = createTempFile("output", ".vcf.gz");
        outputVcf.delete();
        final File outputVcfIndex = new File(outputVcf.getPath() + ".tbi");
        outputVcfIndex.delete();
        outputVcfIndex.deleteOnExit();
        final List<String> args = new ArrayList<String>() {
            private static final long serialVersionUID = 1L;
            {
             add("-" + StandardArgumentDefinitions.VARIANT_SHORT_NAME);
             add(INPUT_VCF.getAbsolutePath());
             add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
             add(outputVcf.getAbsolutePath());
             add("-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME);
             add(b37_reference_20_21);
             add("-" + StandardArgumentDefinitions.INPUT_SHORT_NAME);
             add( INPUT_BAM.getAbsolutePath());
             add("-" + GenotypeStructuralVariantsSpark.FASTQ_FILE_DIR_SHORT_NAME);
             add( INPUT_FASTQ_DIR.getAbsolutePath());
             add("-L");
             add( INPUT_INTERVALS.getAbsolutePath());
            }
        };
        runCommandLine(args);
        Assert.assertTrue(outputVcf.isFile());
    }
}

