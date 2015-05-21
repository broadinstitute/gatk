package org.broadinstitute.hellbender.dev;

import com.google.api.services.genomics.Genomics;
import com.google.api.services.genomics.model.Read;
import com.google.cloud.genomics.dataflow.readers.bam.ReadConverter;
import htsjdk.samtools.*;
import htsjdk.samtools.util.Locatable;
import htsjdk.tribble.Feature;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.CommandLineParser;
import org.broadinstitute.hellbender.dev.pipelines.bqsr.CalibrationTablesBuilder;
import org.broadinstitute.hellbender.dev.tools.walkers.bqsr.BaseRecalibratorUprooted;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.ApplyBQSR;
import org.broadinstitute.hellbender.tools.IntegrationTestSpec;
import org.broadinstitute.hellbender.tools.recalibration.RecalibrationArgumentCollection;
import org.broadinstitute.hellbender.tools.recalibration.RecalibrationTables;
import org.broadinstitute.hellbender.tools.recalibration.covariates.Covariate;
import org.broadinstitute.hellbender.tools.recalibration.covariates.StandardCovariateList;
import org.broadinstitute.hellbender.tools.walkers.bqsr.BaseRecalibrator;
import org.broadinstitute.hellbender.tools.walkers.bqsr.RecalibrationEngine;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import org.testng.asserts.Assertion;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Tracking BaserRecalibratorIntegrationTest, since we're testing the same functionality.
 */
public final class CalibrationTablesBuilderIntegrationTest extends CommandLineProgramTest {

    private static class BQSRTest {
        final String reference;
        final String bam;
        final String knownSites;
        final String args;
        final String expectedFileName;

        private BQSRTest(String reference, String bam, String knownSites, String args, String expectedFileName) {
            this.reference = reference;
            this.bam = bam;
            this.knownSites = knownSites;
            this.args = args;
            this.expectedFileName = expectedFileName;
        }

        public String getCommandLine() {
            return  " -R " + reference +
                    " -I " + bam +
                    " " + args +
                    (knownSites.isEmpty() ? "": " -knownSites " + knownSites) +
                    " --RECAL_TABLE_FILE %s" +
                    " -sortAllCols";
        }


        public List<String> getKnownSites() {
            ArrayList<String> ret = new ArrayList<>();
            if (!knownSites.isEmpty()) {
                ret.add(knownSites);
            }
            String[] split = Utils.escapeExpressions(args);
            for (int i=0; i<split.length; i++) {
                if (split[i].equalsIgnoreCase("-knownSites")) {
                    ret.add(split[i+1]);
                }
            }
            return ret;
        }

        @Override
        public String toString() {
            return String.format("BQSR(bam='%s', args='%s')", bam, args);
        }
    }

    private String getResourceDir(){
        return getTestDataDir() + "/" + "BQSR" + "/";
    }

    @DataProvider(name = "BQSRTest")
    public Object[][] createBQSRTestData() {
        final String hg18Reference = publicTestDir + "human_g1k_v37.chr17_1Mb.fasta";
        final String b36Reference = getResourceDir() + "human_b36_both.chr1_1k.fasta";
        final String HiSeqBam = getResourceDir() + "NA12878.chr17_69k_70k.dictFix.bam";
        final String dbSNPb37 =  getResourceDir() + "dbsnp_132.b37.excluding_sites_after_129.chr17_69k_70k.vcf";
        final String origQualsBam = getResourceDir() + "originalQuals.1kg.chr1.1-1K.1RG.dictFix.bam";
        final String dbSNPb36 = getResourceDir() + "dbsnp_132.b36.excluding_sites_after_129.chr1_1k.vcf";

        final String moreSites = getResourceDir() + "bqsr.fakeSitesForTesting.b37.chr17.vcf"; //for testing 2 input files

        return new Object[][]{
                {new BQSRTest(hg18Reference, HiSeqBam, dbSNPb37, "", getResourceDir() + "expected.NA12878.chr17_69k_70k.txt")},
                {new BQSRTest(hg18Reference, HiSeqBam, dbSNPb37, "-knownSites " + moreSites, getResourceDir() + "expected.NA12878.chr17_69k_70k.2inputs.txt")},
                {new BQSRTest(hg18Reference, HiSeqBam, dbSNPb37, "--indels_context_size 4", getResourceDir() + "expected.NA12878.chr17_69k_70k.indels_context_size4.txt")},
                {new BQSRTest(hg18Reference, HiSeqBam, dbSNPb37, "--low_quality_tail 5", getResourceDir() + "expected.NA12878.chr17_69k_70k.low_quality_tail5.txt")},
                {new BQSRTest(hg18Reference, HiSeqBam, dbSNPb37, "--quantizing_levels 6", getResourceDir() + "expected.NA12878.chr17_69k_70k.quantizing_levels6.txt")},
                {new BQSRTest(hg18Reference, HiSeqBam, dbSNPb37, "--mismatches_context_size 4", getResourceDir() + "expected.NA12878.chr17_69k_70k.mismatches_context_size4.txt")},
                {new BQSRTest(b36Reference, origQualsBam, dbSNPb36, "-OQ", getResourceDir() + "expected.originalQuals.1kg.chr1.1-1K.1RG.dictFix.OQ.txt")},
        };
    }
    @Test(dataProvider = "BQSRTest")
    public void testBQSR(BQSRTest params) throws IOException {
        String outputDest = "test-table-pre.txt";

        // 1. Grab the needed inputs.
        String toolCmdLine = String.format(params.getCommandLine(), outputDest);
        SamReader reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(params.bam));
        SAMFileHeader header = reader.getFileHeader();
        List<Read> input = getReads(reader);
        List<SimpleInterval> knownSites = getKnownSites(params);

        // 2. Run the computation.
        CalibrationTablesBuilder calibrationTablesBuilder = new CalibrationTablesBuilder(header, params.reference, toolCmdLine);
        calibrationTablesBuilder.add(input, knownSites);
        calibrationTablesBuilder.done();
        RecalibrationTables output = calibrationTablesBuilder.getRecalibrationTables();
        StandardCovariateList requestedCovariates = calibrationTablesBuilder.getRequestedCovariates();

        // 3. This is where we'd merge the inputs if we had several.
        RecalibrationEngine.finalizeRecalibrationTables(output);

        // 4. Generate report.
        BaseRecalibratorUprooted baseRecalibratorUprooted = BaseRecalibratorUprooted.fromCommandLine(header, toolCmdLine, System.out);
        baseRecalibratorUprooted.checkClientArguments();
        baseRecalibratorUprooted.onTraversalStart(null);
        baseRecalibratorUprooted.saveReport(output, requestedCovariates);

        // 5. Compare against expected report.
        final File gotTablePre = new File(outputDest);
        final File expectedTablePre = new File(params.expectedFileName);
        IntegrationTestSpec.compareTextFiles(gotTablePre, expectedTablePre);
    }

    @Test
    public void testBQSRFailWithoutDBSNP() throws IOException {
        final String resourceDir =  getTestDataDir() + "/" + "BQSR" + "/";

        final String hg18Reference = publicTestDir + "human_g1k_v37.chr17_1Mb.fasta";
        final String HiSeqBam = resourceDir + "NA12878.chr17_69k_70k.dictFix.bam";

        final String  NO_DBSNP = "";
        final String  NO_ARGS = "";
        final BQSRTest params = new BQSRTest(hg18Reference, HiSeqBam, NO_DBSNP, NO_ARGS, resourceDir + "expected.NA12878.chr17_69k_70k.txt");

        try {
            testBQSR(params);
            Assert.fail("Should have thrown a CommandLineException");
        } catch (UserException.CommandLineException expected) {
            // OK
        }
    }


    @Test
    public void testBQSRFailWithIncompatibleReference() throws IOException {
        final String resourceDir =  getTestDataDir() + "/" + "BQSR" + "/";

        final String HiSeqBam_Hg18 = resourceDir + "HiSeq.1mb.1RG.2k_lines.bam";

        final String  NO_ARGS = "";
        final BQSRTest params = new BQSRTest(hg19MiniReference, HiSeqBam_Hg18, hg19_chr1_1M_dbSNP, NO_ARGS, resourceDir + "expected.txt");
        try {
            testBQSR(params);
            Assert.fail("Should have thrown MissingContigInSequenceDictionary");
        } catch (UserException.MissingContigInSequenceDictionary expected) {
            // OK
        }
    }


    private List<Read> getReads(SamReader reader) {
        List<Read> ret = new ArrayList<>();
        ReadFilter readFilter = BaseRecalibratorUprooted.readFilter();
        for (SAMRecord sr : reader) {
            if (!readFilter.test(sr)) continue;
            try {
                Read e = ReadConverter.makeRead(sr);
                ret.add(e);
            } catch (SAMException x) {
                System.out.println("Skipping read " + sr.getReadName() + " because we can't convert it.");
            } catch (NullPointerException y) {
                System.out.println("Skipping read " + sr.getReadName() + " because we can't convert it. (null?)");
            }
        }
        return ret;
    }

    private List<SimpleInterval> getKnownSites(BQSRTest params) {
        List<SimpleInterval> ret = new ArrayList<>();
        for (String knownSites : params.getKnownSites()) {
            System.out.println("Reading the known sites");
            FeatureDataSource<VariantContext> source = new FeatureDataSource<>(new File(knownSites), new VCFCodec(), "KnownIntervals");
            for (VariantContext foo : source) {
                int start = foo.getStart();
                int end = foo.getEnd();
                String contig = foo.getContig();
                ret.add(new SimpleInterval(contig, start, end));
            }
        }
        return ret;
    }
}
