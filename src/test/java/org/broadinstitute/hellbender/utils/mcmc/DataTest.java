package org.broadinstitute.hellbender.utils.mcmc;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Collections;

/**
 * Unit tests for {@link Data}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class DataTest extends BaseTest {
    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testEmptyDatasetException() {
        new Data<>("empty", Collections.emptyList());
    }

    @Test(expectedExceptions = UserException.class)
    public void testBadInputFileException() {
        new Data<>("badInputFile", new File(publicTestDir + "BAD_INPUT_FILE"), Integer::parseInt);
    }
}