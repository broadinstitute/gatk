package org.broadinstitute.hellbender.tools.picard.analysis.directed;

import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.AfterTest;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileReader;
import java.util.*;

/**
 * @author jgentry@broadinstitute.org
 */

public final class CollectRrbsMetricsTest extends CommandLineProgramTest {
	public static final File CHR_M_SAM = new File(getTestDataDir(), "picard/analysis/directed/CollectRrbsMetrics/chrMReads.sam");
	public static final File CHR_M_REFERENCE = new File(getTestDataDir(), "picard/analysis/directed/CollectRrbsMetrics/chrM.reference.fasta");

	private File rootTestDir;

	@BeforeTest
	private void setUp() throws Exception {
		rootTestDir = BaseTest.createTempFile("crmt.", ".tmp");
		Assert.assertTrue(rootTestDir.delete());
		Assert.assertTrue(rootTestDir.mkdir());
	}

	@AfterTest
	private void tearDown() {
		IOUtil.deleteDirectoryTree(rootTestDir);
	}

	@Test
	public void chrMReads() throws Exception {
		final MetricsFile<RrbsSummaryMetrics, ?> metricsFile = getSummaryFile(
                CHR_M_SAM.getAbsolutePath(), CHR_M_REFERENCE.getAbsolutePath(), rootTestDir + "/READ_TEST", new ArrayList<>());
		final RrbsSummaryMetrics metrics = metricsFile.getMetrics().get(0);
		Assert.assertEquals(metrics.READS_ALIGNED.intValue(), 5);
		Assert.assertEquals(metrics.NON_CPG_BASES.intValue(), 15);
		Assert.assertEquals(metrics.NON_CPG_CONVERTED_BASES.intValue(), 11);
		Assert.assertEquals(metrics.PCT_NON_CPG_BASES_CONVERTED, 0.733333);
		Assert.assertEquals(metrics.CPG_BASES_SEEN.intValue(), 5);
		Assert.assertEquals(metrics.CPG_BASES_CONVERTED.intValue(), 1);
		Assert.assertEquals(metrics.PCT_CPG_BASES_CONVERTED, 0.2);
		Assert.assertEquals(metrics.MEAN_CPG_COVERAGE, 1.666667);
		Assert.assertEquals(metrics.MEDIAN_CPG_COVERAGE.intValue(), 2);
		Assert.assertEquals(metrics.READS_WITH_NO_CPG.intValue(), 1);
		Assert.assertEquals(metrics.READS_IGNORED_SHORT.intValue(), 1);
		Assert.assertEquals(metrics.READS_IGNORED_MISMATCHES.intValue(), 1);
	}

	private MetricsFile<RrbsSummaryMetrics, ?> getSummaryFile(final String input, final String reference, final String prefix,
															  final List<String> sequences) throws Exception {
		final List<String> argList = new ArrayList<>();
		argList.add("--input");
        argList.add(input);
		argList.add("--METRICS_FILE_PREFIX");
        argList.add(prefix);
		argList.add("--reference");
        argList.add(reference);
		for (final String sequence : sequences) {
			argList.add("--SEQUENCE_NAMES");
            argList.add(sequence);
		}

        runCommandLine(argList);

		final MetricsFile<RrbsSummaryMetrics, ?> retVal = new MetricsFile<RrbsSummaryMetrics, Integer>();
		retVal.read(new FileReader(prefix + ".rrbs_summary_metrics"));
		return retVal;
	}


}
