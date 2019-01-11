package org.broadinstitute.hellbender.tools.walkers.validation;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.engine.AbstractConcordanceWalker;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Created by Takuto Sato on 1/31/17.
 */
public class ConcordanceIntegrationTest extends CommandLineProgramTest{
    final double epsilon = 1e-3;

    @Test
    public void testSimple() throws Exception {
        final String testDir =  toolsTestDir + "concordance/";
        final File evalVcf = new File(testDir, "gatk4-dream3-mini.vcf");
        final File truthVcf = new File(testDir, "gatk3-dream3-mini.vcf");
        final File summary = createTempFile("summary", ".txt");

        final String[] args = {
                "--" + AbstractConcordanceWalker.EVAL_VARIANTS_LONG_NAME, evalVcf.toString(),
                "--" + AbstractConcordanceWalker.TRUTH_VARIANTS_LONG_NAME, truthVcf.toString(),
                "--" + Concordance.SUMMARY_LONG_NAME, summary.toString(),
        };
        runCommandLine(args);

        ConcordanceSummaryRecord.Reader reader = new ConcordanceSummaryRecord.Reader(summary);
        ConcordanceSummaryRecord snpRecord = reader.readRecord();
        ConcordanceSummaryRecord indelRecord = reader.readRecord();

        Assert.assertEquals(snpRecord.getVariantType(), VariantContext.Type.SNP);
        Assert.assertEquals(indelRecord.getVariantType(), VariantContext.Type.INDEL);

        // numbers verified by manual inspection
        Assert.assertEquals(snpRecord.getTruePositives(), 1);
        Assert.assertEquals(indelRecord.getTruePositives(), 4);
        Assert.assertEquals(snpRecord.getFalsePositives(), 0);
        Assert.assertEquals(indelRecord.getFalsePositives(), 3);
        Assert.assertEquals(snpRecord.getFalseNegatives(), 3);
        Assert.assertEquals(indelRecord.getFalseNegatives(), 2);
        Assert.assertEquals(snpRecord.getSensitivity(), 1.0/4, epsilon);
        Assert.assertEquals(snpRecord.getPrecision(), 1.0, epsilon);

    }

    // Test going from an integer chromosome (22) to a character chromosome (X)
    @Test
    public void testt22X() throws Exception {
        final String testDir = toolsTestDir + "concordance/";
        final File truthVcf = new File(testDir, "gatk3-dream3-x.vcf");
        final File evalVcf = new File(testDir, "gatk4-dream3-x.vcf");
        final File tpfp = createTempFile("tpfp", ".vcf");
        final File tpfn = createTempFile("tpfn", ".vcf");
        final File ftnfn = createTempFile("ftnfn", ".vcf");
        final File summary = createTempFile("summary", ".txt");

        final String[] args = {
                "--" + AbstractConcordanceWalker.EVAL_VARIANTS_LONG_NAME, evalVcf.toString(),
                "--" + AbstractConcordanceWalker.TRUTH_VARIANTS_LONG_NAME, truthVcf.toString(),
                "--" + Concordance.SUMMARY_LONG_NAME, summary.toString(),
                "--" + Concordance.TRUE_POSITIVES_AND_FALSE_NEGATIVES_LONG_NAME, tpfn.getAbsolutePath(),
                "--" + Concordance.TRUE_POSITIVES_AND_FALSE_POSITIVES_LONG_NAME, tpfp.getAbsolutePath(),
                "--" + Concordance.FILTERED_TRUE_NEGATIVES_AND_FALSE_NEGATIVES_LONG_NAME, ftnfn.getAbsolutePath()
        };
        runCommandLine(args);
        ConcordanceSummaryRecord.Reader summaryReader = new ConcordanceSummaryRecord.Reader(summary);
        ConcordanceSummaryRecord snpRecord = summaryReader.readRecord();
        ConcordanceSummaryRecord indelRecord = summaryReader.readRecord();

        // Takuto painstakingly counted the correct values
        Assert.assertEquals(snpRecord.getTruePositives(), 2);
        Assert.assertEquals(indelRecord.getTruePositives(), 2);

        Assert.assertEquals(snpRecord.getFalsePositives(), 1);
        Assert.assertEquals(indelRecord.getFalsePositives(), 5);

        Assert.assertEquals(snpRecord.getFalseNegatives(), 6);
        Assert.assertEquals(indelRecord.getFalseNegatives(), 0);

        // Test the output vcfs
        final List<VariantContext> truePositivesAndFalseNegatives = VariantContextTestUtils.streamVcf(tpfn).collect(Collectors.toList());
        final List<VariantContext> truePositivesAndFalsePositives = VariantContextTestUtils.streamVcf(tpfp).collect(Collectors.toList());

        Assert.assertEquals(truePositivesAndFalseNegatives.size(), 10);
        Assert.assertEquals(truePositivesAndFalsePositives.size(), 10);

        final VariantContext firstVc = truePositivesAndFalsePositives.get(0);

        Assert.assertEquals(firstVc.getContig(), "22");
        Assert.assertEquals(firstVc.getStart(), 50357818);
        Assert.assertEquals(firstVc.getAttribute(Concordance.TRUTH_STATUS_VCF_ATTRIBUTE), ConcordanceState.TRUE_POSITIVE.getAbbreviation());
        Assert.assertEquals(firstVc.getType(), VariantContext.Type.SNP);
    }

    // Truth and eval are the same file
    // The truth file only contains the PASS sites; should get 100% sensitivity and specificity.
    @Test
    public void testPerfectMatchVcf() throws Exception {
        final String testDir = toolsTestDir + "concordance/";
        final File truthVcf = new File(testDir, "same-truth.vcf");
        final File evalVcf = new File(testDir, "same-eval.vcf");
        final File summary = createTempFile("summary", ".txt");

        final String[] args = {
                "--" + AbstractConcordanceWalker.EVAL_VARIANTS_LONG_NAME, evalVcf.toString(),
                "--" + AbstractConcordanceWalker.TRUTH_VARIANTS_LONG_NAME, truthVcf.toString(),
                "--" + Concordance.SUMMARY_LONG_NAME, summary.toString(),
        };

        runCommandLine(args);
        ConcordanceSummaryRecord.Reader reader = new ConcordanceSummaryRecord.Reader(summary);
        ConcordanceSummaryRecord snpRecord = reader.readRecord();
        ConcordanceSummaryRecord indelRecord = reader.readRecord();

        Assert.assertEquals(snpRecord.getTruePositives() + indelRecord.getTruePositives(), 284);
        Assert.assertEquals(snpRecord.getFalsePositives(), 0);
        Assert.assertEquals(indelRecord.getFalsePositives(), 0);
        Assert.assertEquals(snpRecord.getFalseNegatives(), 0);
        Assert.assertEquals(indelRecord.getFalseNegatives(), 0);
        Assert.assertEquals(snpRecord.getSensitivity(), 1.0, epsilon);
        Assert.assertEquals(snpRecord.getPrecision(), 1.0, epsilon);
        Assert.assertEquals(indelRecord.getSensitivity(), 1.0, epsilon);
        Assert.assertEquals(indelRecord.getPrecision(), 1.0, epsilon);
    }

    @Test
    public void testFilterAnalysis() throws Exception {
        final String testDir = toolsTestDir + "concordance/";
        final File truthVcf = new File(testDir, "filter-analysis-truth.vcf");
        final File evalVcf = new File(testDir, "filter-analysis-eval.vcf");
        final File summary = createTempFile("summary", ".txt");
        final File filterAnalysis = createTempFile("filter-analysis", ".txt");

        final String[] args = {
                "--" + AbstractConcordanceWalker.EVAL_VARIANTS_LONG_NAME, evalVcf.getAbsolutePath(),
                "--" + AbstractConcordanceWalker.TRUTH_VARIANTS_LONG_NAME, truthVcf.getAbsolutePath(),
                "--" + Concordance.SUMMARY_LONG_NAME, summary.getAbsolutePath(),
                "--" + Concordance.FILTER_ANALYSIS_LONG_NAME, filterAnalysis.getAbsolutePath()
        };

        runCommandLine(args);
        final List<FilterAnalysisRecord> filterAnalysisRecords = FilterAnalysisRecord.readFromFile(filterAnalysis);
        Assert.assertEquals(filterAnalysisRecords.size(), 3);

        // this ensures that the first record is average_filter, then bad_filter, then good_filter
        Collections.sort(filterAnalysisRecords, Comparator.comparing(FilterAnalysisRecord::getFilter));
        final FilterAnalysisRecord average = filterAnalysisRecords.get(0);
        final FilterAnalysisRecord bad = filterAnalysisRecords.get(1);
        final FilterAnalysisRecord good = filterAnalysisRecords.get(2);

        Assert.assertEquals(average.getFalseNegativeCount(),1);
        Assert.assertEquals(average.getUniqueFalseNegativeCount(),0);
        Assert.assertEquals(average.getTrueNegativeCount(), 2);
        Assert.assertEquals(average.getUniqueTrueNegativeCount(), 1);

        Assert.assertEquals(bad.getFalseNegativeCount(),4);
        Assert.assertEquals(bad.getUniqueFalseNegativeCount(),3);
        Assert.assertEquals(bad.getTrueNegativeCount(), 0);
        Assert.assertEquals(bad.getUniqueTrueNegativeCount(), 0);

        Assert.assertEquals(good.getFalseNegativeCount(),0);
        Assert.assertEquals(good.getUniqueFalseNegativeCount(),0);
        Assert.assertEquals(good.getTrueNegativeCount(), 5);
        Assert.assertEquals(good.getUniqueTrueNegativeCount(), 4);



    }

    @Test
    public void testDreamSensitivity() throws Exception {
        final String testDir = toolsTestDir + "concordance/";
        final File evalVcf = new File(testDir, "dream3-chr21.vcf");
        final File truthVcf = new File(testDir, "dream3-truth-minus-SV-chr21.vcf");
        final File summary = createTempFile("summary", ".txt");

        final String[] args = {
                "--" + AbstractConcordanceWalker.EVAL_VARIANTS_LONG_NAME, evalVcf.toString(),
                "--" + AbstractConcordanceWalker.TRUTH_VARIANTS_LONG_NAME, truthVcf.toString(),
                "--" + Concordance.SUMMARY_LONG_NAME, summary.toString(),
        };

        runCommandLine(args);

        ConcordanceSummaryRecord.Reader reader = new ConcordanceSummaryRecord.Reader(summary);
        ConcordanceSummaryRecord snpRecord = reader.readRecord();
        ConcordanceSummaryRecord indelRecord = reader.readRecord();

        // The expected sensitivities come from dream challenge's official python script
        // We don't expect them to match exactly so give it a generous delta
        // Currently we tend to overestimate sensitivities TODO: investigate
        Assert.assertEquals(snpRecord.getSensitivity(), 0.879, 0.06); // with python script (masked, chr21) you get 0.879
        Assert.assertEquals(indelRecord.getSensitivity(), 0.791, 0.12); // with python script (masked) you get 0.791

        // count the number of INDEL variants in the truth using awk. it should match # TP + # FN
        // grep -v ^# dream3-truth-minus-MSK-chr21.vcf | grep -v SVTYPE | awk 'length($4) != length($5) {print $0 }' | wc -l
        // 94
        Assert.assertEquals(indelRecord.getTruePositives() + indelRecord.getFalseNegatives(), 94);

        // analogously for SNPs:
        // grep -v ^# dream3-truth-minus-MSK-chr21.vcf | grep -v SVTYPE | awk 'length($4) == length($5) {print $0 }' | wc -l
        // 78
        Assert.assertEquals(snpRecord.getTruePositives() + snpRecord.getFalseNegatives(), 78);
    }
}