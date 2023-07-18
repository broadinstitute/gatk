package org.broadinstitute.hellbender.tools.reference;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.alignment.MummerExecutor;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public class CompareReferencesUnitTest extends GATKBaseTest {

    @Test
    public void testGenerateFastaForSequence() throws IOException {
        GATKPath ref = new GATKPath(getToolTestDataDir() + "hg19mini.fasta");
        File expectedOutput = new File(getToolTestDataDir() + "1.fasta");
        String sequenceName = "1";
        File output = createTempFile("example_chr1", ".fasta");

        CompareReferences.generateFastaForSequence(ref, sequenceName, new GATKPath(output.toString()));

        IntegrationTestSpec.assertEqualTextFiles(output, expectedOutput);
    }

}