package org.broadinstitute.hellbender.tools.walkers.sv;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.testng.annotations.Test;

import java.io.File;

import static org.testng.Assert.*;

public class SVAnnotateIntegrationTest extends CommandLineProgramTest {
    final File inputVCF = new File(getToolTestDataDir() + "integration.vcf.gz");
    private final String largeFileDirectory = GATKBaseTest.largeFileTestDir + "SVAnnotate/";
    final File gtf = new File(largeFileDirectory + "MANE.selected.GRCh38.v0.95.select_ensembl_genomic.gtf");
    final File noncodingElements = new File(largeFileDirectory + "noncoding.selected.hg38.bed.gz");

    @Test
    public void testBasicOptions() {
        final File output = createTempFile("annotated",".vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addOutput(output)
                .addVCF(inputVCF)
                .add("proteinCodingGTF", gtf)
                .add("nonCodingBed", noncodingElements);

        runCommandLine(args, SVAnnotate.class.getSimpleName());

    }
}