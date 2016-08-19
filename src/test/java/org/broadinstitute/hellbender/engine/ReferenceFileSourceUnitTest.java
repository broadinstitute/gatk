package org.broadinstitute.hellbender.engine;

import org.broadinstitute.hellbender.engine.datasources.ReferenceFileSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.test.BaseTest;

import org.testng.annotations.Test;

import java.io.IOException;

public class ReferenceFileSourceUnitTest extends BaseTest {

    @Test(expectedExceptions = UserException.MissingReference.class)
    public void testMissingReferenceFile() throws IOException {
        new ReferenceFileSource(BaseTest.getSafeNonExistentFile("NonExistentReference.fasta").getAbsolutePath());
    }

}
