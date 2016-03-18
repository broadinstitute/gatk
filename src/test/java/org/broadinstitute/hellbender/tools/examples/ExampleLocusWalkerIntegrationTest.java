/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2015 Daniel Gómez-Sánchez
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

package org.broadinstitute.hellbender.tools.examples;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.Arrays;

/**
 * @author Daniel Gómez-Sánchez (magicDGS)
 */
public class ExampleLocusWalkerIntegrationTest extends CommandLineProgramTest {
	private static final String TEST_DATA_DIRECTORY = publicTestDir + "org/broadinstitute/hellbender/engine/";
	private static final String TEST_OUTPUT_DIRECTORY = publicTestDir + "org/broadinstitute/hellbender/tools/examples/";

	@Test
	public void testExampleVariantWalker() throws IOException {
		IntegrationTestSpec testSpec = new IntegrationTestSpec(
			" -L 1" +
				" -R " + hg19MiniReference +
				" -I " + TEST_DATA_DIRECTORY + "reads_data_source_test1.bam" +
				" -V " + TEST_DATA_DIRECTORY + "feature_data_source_test.vcf" +
				" -O %s",
			Arrays.asList(TEST_OUTPUT_DIRECTORY + "expected_ExampleLocusWalkerIntegrationTest_output.txt")
		);
		testSpec.executeTest("testExampleLocusWalker", this);
	}

}
