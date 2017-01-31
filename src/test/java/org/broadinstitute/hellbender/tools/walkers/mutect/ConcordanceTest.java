package org.broadinstitute.hellbender.tools.walkers.mutect;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

/**
 * Created by tsato on 1/31/17.
 */
public class ConcordanceTest extends CommandLineProgramTest{
    final double epsilon = 1e-3;

    @Test
    public void testSimple() throws Exception {
        final String testDir =  publicTestDir + "org/broadinstitute/hellbender/tools/concordance/";
        final File evalVcf = new File(testDir + "gatk4-dream3-mini.vcf");
        final File truthVcf = new File(testDir + "gatk3-dream3-mini.vcf");
        final File variantTable = createTempFile("table", ".txt");
        final File summary = createTempFile("summary", ".txt");

        final String[] args = {
                "--eval", evalVcf.toString(),
                "--truth", truthVcf.toString(),
                "--table",  variantTable.toString(),
                "--summary", summary.toString(),
                "-tumor", "IS3.snv.indel.sv"
        };
        runCommandLine(args);

        SummaryTableReader reader = new SummaryTableReader(summary);
        SummaryRecord record = reader.readRecord();
        Assert.assertEquals(record.getTruePositives(), 5);
        Assert.assertEquals(record.getFalsePositives(), 3);
        Assert.assertEquals(record.getFalseNegatives(), 5);
        Assert.assertTrue(Math.abs(record.getSensitivity() - 5.0/10) < epsilon);
        Assert.assertTrue(Math.abs(record.getPrecision() - 5.0/8) < epsilon);
    }

    // Test going from an integer chromosome (22) to a character chromosome (X)
    @Test
    public void test22X() throws Exception {
        final String testDir = publicTestDir + "org/broadinstitute/hellbender/tools/concordance/";
        final File truthVcf = new File(testDir + "gatk3-dream3-x.vcf");
        final File evalVcf = new File(testDir + "gatk4-dream3-x.vcf");
        final File table = createTempFile("table", ".txt");
        final File summary = createTempFile("summary", ".txt");

        final String[] args = {
                "--eval", evalVcf.toString(),
                "--truth", truthVcf.toString(),
                "--table", table.toString(),
                "--summary", summary.toString(),
                "-tumor", "IS3.snv.indel.sv"
        };
        runCommandLine(args);
        SummaryTableReader reader = new SummaryTableReader(summary);
        SummaryRecord record = reader.readRecord();

        Assert.assertEquals(record.getTruePositives(), 4);
        Assert.assertEquals(record.getFalsePositives(), 6);
        Assert.assertEquals(record.getFalseNegatives(), 6);
        Assert.assertTrue(Math.abs(record.getSensitivity() - 4.0/10) < epsilon);
        Assert.assertTrue(Math.abs(record.getPrecision() - 4.0/10) < epsilon);
    }

    @Test
    public void testPerfectMatchVcf() throws Exception {
        // Truth and eval come from the same vcf file.
        // The truth file only contains the PASS sites. Should get 100% sensitivity and specificity.
        final String testDir = publicTestDir + "org/broadinstitute/hellbender/tools/concordance/";
        final File truthVcf = new File(testDir + "same-truth.vcf");
        final File evalVcf = new File(testDir + "same-eval.vcf");
        final File table = createTempFile("table", ".txt");
        final File summary = createTempFile("summary", ".txt");

        final String[] args = {
                "--eval", evalVcf.toString(),
                "--truth", truthVcf.toString(),
                "--table", table.toString(),
                "--summary", summary.toString(),
                "-tumor", "TUMOR"
        };

        runCommandLine(args);
        SummaryTableReader reader = new SummaryTableReader(summary);
        SummaryRecord record = reader.readRecord();

        Assert.assertEquals(record.getTruePositives(), 284);
        Assert.assertEquals(record.getFalsePositives(), 0);
        Assert.assertEquals(record.getFalseNegatives(), 0);
        Assert.assertTrue(Math.abs(record.getSensitivity() - 1.0) < epsilon);
        Assert.assertTrue(Math.abs(record.getPrecision() - 1.0) < epsilon);
    }

    private class SummaryTableReader extends TableReader<SummaryRecord> {
        private SummaryTableReader(final File summaryTable) throws IOException {
            super(summaryTable);
        }

        @Override
        protected SummaryRecord createRecord(final DataLine dataLine) {
            final long truePositives = Long.parseLong(dataLine.get(SummaryRecord.TRUE_POSITIVE_COLUMN_NAME));
            final long falsePositives = Long.parseLong(dataLine.get(SummaryRecord.FALSE_POSITIVE_COLUMN_NAME));
            final long falseNegatives = Long.parseLong(dataLine.get(SummaryRecord.FALSE_NEGATIVE_COLUMN_NAME));

            return new SummaryRecord(truePositives, falsePositives, falseNegatives);
        }
    }

    // TODO: test multi allelic
    // TODO: test high confidence intervals
    // TODO; test table writers
    // TODO: test xy chromosomes
    // TODO: test MT chromosomes

}