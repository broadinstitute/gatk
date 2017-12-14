package org.broadinstitute.hellbender.engine;

import org.broadinstitute.hellbender.engine.datasources.ReferenceFileSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.GATKBaseTest;

import org.testng.annotations.Test;

import java.io.IOException;

public class ReferenceFileSourceUnitTest extends GATKBaseTest {

    @Test(expectedExceptions = UserException.MissingReference.class)
    public void testMissingReferenceFile() throws IOException {
        new ReferenceFileSource(GATKBaseTest.getSafeNonExistentFile("NonExistentReference.fasta").getAbsolutePath());
    }

}
