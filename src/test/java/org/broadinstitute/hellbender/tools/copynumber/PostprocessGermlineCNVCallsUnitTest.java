package org.broadinstitute.hellbender.tools.copynumber;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.python.PythonUnitTestRunner;
import org.testng.annotations.Test;

import java.util.Collections;

/**
 * Note that this class runs tests from the gcnvkernel module, so if changes are made to the python tests, the python
 * package must be rebuilt to reflect those changes, e.g.:
 * ./gradlew pythonPackageArchive
 * pip install --upgrade $(gatk_repo_path)/build/gatkPythonPackageArchive.zip
 */
public class PostprocessGermlineCNVCallsUnitTest extends GATKBaseTest {
    @Test(groups = "python")
    public void testPythonVCFReading() {
        final String testFile = "gcnvkernel.io.test_io_vcf_parsing";

        final PythonUnitTestRunner runner = new PythonUnitTestRunner();
        runner.setup(Collections.singletonList("gcnvkernel"), Collections.emptyList(), Collections.emptyList());
        runner.runTest(testFile);
    }

    @Test(groups = "python")
    public void testSegmentationEngine() {
        final String testFile = "gcnvkernel.postprocess.test_viterbiSegmentationEngine";

        final PythonUnitTestRunner runner = new PythonUnitTestRunner();
        runner.setup(Collections.singletonList("gcnvkernel"), Collections.emptyList(), Collections.emptyList());
        runner.runTest(testFile);
    }

}