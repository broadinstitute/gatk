package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.Main;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Optional;

import static org.testng.Assert.*;

public class GetNormalArtifactDataIntegrationTest extends CommandLineProgramTest {

    private static final String DREAM_BAMS_DIR = largeFileTestDir + "mutect/dream_synthetic_bams/";
    private static final File DREAM_4_NORMAL = new File(DREAM_BAMS_DIR, "normal_4.bam");
    private static final File DREAM_4_TUMOR = new File(DREAM_BAMS_DIR, "tumor_4.bam");

    @Test
    public void test() {
        Utils.resetRandomGenerator();
        final File tumor = DREAM_4_TUMOR;
        final File normal = DREAM_4_NORMAL;
        final String normalSample = getSampleName(normal);
        final File output = createTempFile("output", ".table");

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addInput(tumor)
                .addInput(normal)
                .addOutput(output)
                .addReference(b37Reference)
                .addInterval("20")
                .add(M2ArgumentCollection.NORMAL_SAMPLE_LONG_NAME, normalSample);

        runCommandLine(args);

        final List<NormalArtifactRecord> records = NormalArtifactRecord.readFromFile(output);

        final File gather = createTempFile("gather", ".table");
        final ArgumentsBuilder args2 = new ArgumentsBuilder()
                .addInput(output)
                .addInput(output)
                .addOutput(gather);

        runCommandLine(args2, GatherNormalArtifactData.class.getSimpleName());
    }

    private String getSampleName(final File bam)  {
        try {
            final File nameFile = createTempFile("sample_name", ".txt");
            new Main().instanceMain(makeCommandLineArgs(Arrays.asList("-I", bam.getAbsolutePath(), "-O", nameFile.getAbsolutePath(), "-encode"), "GetSampleName"));
            return Files.readAllLines(nameFile.toPath()).get(0);
        } catch (final IOException ex) {
            throw new IllegalArgumentException(ex);
        }
    }

}