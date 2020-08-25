package org.broadinstitute.hellbender.tools.walkers.validation;

import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.Histogram;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.engine.AbstractConcordanceWalker;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.sam.util.Pair;

import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Created by Michael Gatzen on 9/9/20.
 */
public class ReferenceBlockConcordanceIntegrationTest extends CommandLineProgramTest{

    private static final String CONCORDANCE_TEST_DIR = toolsTestDir + "concordance/";
    private static final String HAPLOTYPECALLER_TEST_DIR = toolsTestDir + "haplotypecaller/";

    @Test
    public void testIdentical() throws Exception {
        // TODO Is it fine to use a test file from a different tool or should I copy it?
        final File truthVcf = new File(HAPLOTYPECALLER_TEST_DIR, "expected.testGVCFMode.gatk4.g.vcf");
        final File evalVcf = new File(HAPLOTYPECALLER_TEST_DIR, "expected.testGVCFMode.gatk4.g.vcf");
        final Path truthBlockHistogramFile = createTempPath("truth_block_histogram", ".tsv");
        final Path evalBlockHistogramFile = createTempPath("eval_block_histogram", ".tsv");
        final Path confidenceConcordanceHistogramFile = createTempPath("confidence_concordance_histogram", ".tsv");

        final String[] args = {
                "--" + AbstractConcordanceWalker.EVAL_VARIANTS_LONG_NAME, evalVcf.toString(),
                "--" + AbstractConcordanceWalker.TRUTH_VARIANTS_LONG_NAME, truthVcf.toString(),
                "--" + ReferenceBlockConcordance.TRUTH_BLOCK_HISTOGRAM_LONG_NAME, truthBlockHistogramFile.toString(),
                "--" + ReferenceBlockConcordance.EVAL_BLOCK_HISTOGRAM_LONG_NAME, evalBlockHistogramFile.toString(),
                "--" + ReferenceBlockConcordance.CONFIDENCE_CONCORDANCE_HISTOGRAM_LONG_NAME, confidenceConcordanceHistogramFile.toString(),
        };
        runCommandLine(args);

        MetricsFile<?, String> truthBlockMetrics = new MetricsFile<>();
        truthBlockMetrics.read(new FileReader(truthBlockHistogramFile.toFile()));
        MetricsFile<?, String> evalBlockMetrics = new MetricsFile<>();
        evalBlockMetrics.read(new FileReader(evalBlockHistogramFile.toFile()));

        Assert.assertEquals(truthBlockMetrics.getNumHistograms(), 1);
        Assert.assertEquals(evalBlockMetrics.getNumHistograms(), 1);

        Histogram<String> truthBlockHistogram = truthBlockMetrics.getHistogram();
        Histogram<String> evalBlockHistogram = evalBlockMetrics.getHistogram();

        // Got this number by counting the <NON_REF> alt alleles in the test GVCF file
        Assert.assertEquals(truthBlockHistogram.getSumOfValues(), 1034);
        Assert.assertEquals(evalBlockHistogram.getSumOfValues(), 1034);

        // Check block histograms both ways, in case one histogram has more entries than the other
        truthBlockHistogram.values().forEach(bin -> {
            Assert.assertTrue(evalBlockHistogram.containsKey(bin.getId()));
            Assert.assertEquals(bin.getValue(), evalBlockHistogram.get(bin.getId()).getValue());
        });

        evalBlockHistogram.values().forEach(bin -> {
            Assert.assertTrue(truthBlockHistogram.containsKey(bin.getId()));
            Assert.assertEquals(bin.getValue(), truthBlockHistogram.get(bin.getId()).getValue());
        });

        MetricsFile<?, String> confidenceConcordanceMetrics = new MetricsFile<>();
        confidenceConcordanceMetrics.read(new FileReader(confidenceConcordanceHistogramFile.toFile()));
        Assert.assertEquals(confidenceConcordanceMetrics.getNumHistograms(), 1);
        Histogram<String> confidenceConcordanceHistogram = confidenceConcordanceMetrics.getHistogram();

        confidenceConcordanceHistogram.values().forEach(bin -> {
            String[] confidenceValues = bin.getId().split(",");
            Assert.assertEquals(confidenceValues[0], confidenceValues[1]);
        });
    }

    private Pair<File, File> writeTestGVCFs(List<TestReferenceBlockConcordanceVariant> truthVariants, List<TestReferenceBlockConcordanceVariant> evalVariants) throws Exception {
        File truthFile = createTempFile("truth", ".gvcf");
        FileWriter writer = new FileWriter(truthFile);
        writer.write("##fileformat=VCFv4.2\n");
        writer.write("##contig=<ID=test_contig,length=1000>\n");
        writer.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTESTSAMPLE\n");
        for (TestReferenceBlockConcordanceVariant variant : truthVariants) {
            writer.write(String.format("test_contig\t%s\t.\tA\t%s\t%s\t.\t%s\tGT:GQ\t%s:%s\n",
                    variant.start,
                    variant.getAltAllele(),
                    variant.getAltAllele().equals("<NON_REF>") ? "." : "1",
                    variant.getAltAllele().equals("<NON_REF>") ? String.format("END=%s", variant.getStop()) : ".",
                    variant.getAltAllele().equals("<NON_REF>") ? "0/0" : "0/1",
                    variant.getConfidence()));
        }
        writer.close();

        File evalFile = createTempFile("eval", ".gvcf");
        writer = new FileWriter(evalFile);
        writer.write("##fileformat=VCFv4.2\n");
        writer.write("##contig=<ID=test_contig,length=1000>\n");
        writer.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTESTSAMPLE\n");
        for (TestReferenceBlockConcordanceVariant variant : evalVariants) {
            writer.write(String.format("test_contig\t%s\t.\tA\t%s\t%s\t.\t%s\tGT:GQ\t%s:%s\n",
                    variant.start,
                    variant.getAltAllele(),
                    variant.getAltAllele().equals("<NON_REF>") ? "." : "1",
                    variant.getAltAllele().equals("<NON_REF>") ? String.format("END=%s", variant.getStop()) : ".",
                    variant.getAltAllele().equals("<NON_REF>") ? "0/0" : "0/1",
                    variant.getConfidence()));
        }
        writer.close();

        return new Pair<>(truthFile, evalFile);
    }

    @DataProvider
    public Object[][] provideSyntheticGVCFs() {
        return new Object[][] {
                // No non_ref blocks
                {
                        // Truth variants
                        Arrays.asList(
                                new TestReferenceBlockConcordanceVariant("C", 1, 1, 99)
                        ),
                        // Eval variants
                        Arrays.asList(
                                new TestReferenceBlockConcordanceVariant("<NON_REF>", 1, 1, 99)
                        ),
                        // Truth block histogram
                        Stream.of(new Object[][] {

                        }).collect(Collectors.toMap(data -> (String) data[0], data -> (Integer) data[1])),
                        // Eval block histogram
                        Stream.of(new Object[][] {
                                { "1,99", 1}
                        }).collect(Collectors.toMap(data -> (String) data[0], data -> (Integer) data[1])),
                        // Confidence concordance
                        Stream.of(new Object[][] {

                        }).collect(Collectors.toMap(data -> (String) data[0], data -> (Integer) data[1])),
                },
                // Two non_ref blocks
                {
                        // Truth variants
                        Arrays.asList(
                                new TestReferenceBlockConcordanceVariant("<NON_REF>", 1, 1, 99)
                        ),
                        // Eval variants
                        Arrays.asList(
                                new TestReferenceBlockConcordanceVariant("<NON_REF>", 1, 1, 99)
                        ),
                        // Truth block histogram
                        Stream.of(new Object[][] {
                                { "1,99", 1}
                        }).collect(Collectors.toMap(data -> (String) data[0], data -> (Integer) data[1])),
                        // Eval block histogram
                        Stream.of(new Object[][] {
                                { "1,99", 1}
                        }).collect(Collectors.toMap(data -> (String) data[0], data -> (Integer) data[1])),
                        // Confidence concordance
                        Stream.of(new Object[][] {
                                { "99,99", 1}
                        }).collect(Collectors.toMap(data -> (String) data[0], data -> (Integer) data[1])),
                },
                // Start with gap then variant
                {
                        // Truth variants
                        Arrays.asList(
                                new TestReferenceBlockConcordanceVariant("<NON_REF>", 1, 10, 99)
                        ),
                        // Eval variants
                        Arrays.asList(
                                new TestReferenceBlockConcordanceVariant("C", 5, 5, 99),
                                new TestReferenceBlockConcordanceVariant("<NON_REF>", 6, 10, 99)
                        ),
                        // Truth block histogram
                        Stream.of(new Object[][] {
                                { "10,99", 1}
                        }).collect(Collectors.toMap(data -> (String) data[0], data -> (Integer) data[1])),
                        // Eval block histogram
                        Stream.of(new Object[][] {
                                { "5,99", 1}
                        }).collect(Collectors.toMap(data -> (String) data[0], data -> (Integer) data[1])),
                        // Confidence concordance
                        Stream.of(new Object[][] {
                                { "99,99", 5}
                        }).collect(Collectors.toMap(data -> (String) data[0], data -> (Integer) data[1])),
                },
                // Start with gap then non_ref
                {
                        // Truth variants
                        Arrays.asList(
                                new TestReferenceBlockConcordanceVariant("<NON_REF>", 1, 10, 99)
                        ),
                        // Eval variants
                        Arrays.asList(
                                new TestReferenceBlockConcordanceVariant("<NON_REF>", 6, 10, 99)
                        ),
                        // Truth block histogram
                        Stream.of(new Object[][] {
                                { "10,99", 1}
                        }).collect(Collectors.toMap(data -> (String) data[0], data -> (Integer) data[1])),
                        // Eval block histogram
                        Stream.of(new Object[][] {
                                { "5,99", 1}
                        }).collect(Collectors.toMap(data -> (String) data[0], data -> (Integer) data[1])),
                        // Confidence concordance
                        Stream.of(new Object[][] {
                                { "99,99", 5}
                        }).collect(Collectors.toMap(data -> (String) data[0], data -> (Integer) data[1])),
                },
                // Start with single block
                {
                        // Truth variants
                        Arrays.asList(
                                new TestReferenceBlockConcordanceVariant("<NON_REF>", 11, 20, 99)
                        ),
                        // Eval variants
                        Arrays.asList(
                                new TestReferenceBlockConcordanceVariant("<NON_REF>", 3, 6, 98),
                                new TestReferenceBlockConcordanceVariant("<NON_REF>", 11, 20, 99)
                        ),
                        // Truth block histogram
                        Stream.of(new Object[][] {
                                { "10,99", 1}
                        }).collect(Collectors.toMap(data -> (String) data[0], data -> (Integer) data[1])),
                        // Eval block histogram
                        Stream.of(new Object[][] {
                                { "4,98", 1},
                                { "10,99", 1}
                        }).collect(Collectors.toMap(data -> (String) data[0], data -> (Integer) data[1])),
                        // Confidence concordance
                        Stream.of(new Object[][] {
                                { "99,99", 10}
                        }).collect(Collectors.toMap(data -> (String) data[0], data -> (Integer) data[1])),
                },
                // No overlap at all
                {
                        // Truth variants
                        Arrays.asList(
                                new TestReferenceBlockConcordanceVariant("<NON_REF>", 11, 20, 99)
                        ),
                        // Eval variants
                        Arrays.asList(
                                new TestReferenceBlockConcordanceVariant("<NON_REF>", 3, 6, 98)
                        ),
                        // Truth block histogram
                        Stream.of(new Object[][] {
                                { "10,99", 1}
                        }).collect(Collectors.toMap(data -> (String) data[0], data -> (Integer) data[1])),
                        // Eval block histogram
                        Stream.of(new Object[][] {
                                { "4,98", 1}
                        }).collect(Collectors.toMap(data -> (String) data[0], data -> (Integer) data[1])),
                        // Confidence concordance
                        Stream.of(new Object[][] {
                        }).collect(Collectors.toMap(data -> (String) data[0], data -> (Integer) data[1])),
                },
                // End with var
                {
                        // Truth variants
                        Arrays.asList(
                                new TestReferenceBlockConcordanceVariant("<NON_REF>", 1, 10, 99)
                        ),
                        // Eval variants
                        Arrays.asList(
                                new TestReferenceBlockConcordanceVariant("<NON_REF>", 1, 9, 99),
                                new TestReferenceBlockConcordanceVariant("C", 10, 10, 99)
                        ),
                        // Truth block histogram
                        Stream.of(new Object[][] {
                                { "10,99", 1}
                        }).collect(Collectors.toMap(data -> (String) data[0], data -> (Integer) data[1])),
                        // Eval block histogram
                        Stream.of(new Object[][] {
                                { "9,99", 1}
                        }).collect(Collectors.toMap(data -> (String) data[0], data -> (Integer) data[1])),
                        // Confidence concordance
                        Stream.of(new Object[][] {
                                { "99,99", 9}
                        }).collect(Collectors.toMap(data -> (String) data[0], data -> (Integer) data[1])),
                },
                // End with var then gap
                {
                        // Truth variants
                        Arrays.asList(
                                new TestReferenceBlockConcordanceVariant("<NON_REF>", 1, 10, 99)
                        ),
                        // Eval variants
                        Arrays.asList(
                                new TestReferenceBlockConcordanceVariant("<NON_REF>", 1, 7, 99),
                                new TestReferenceBlockConcordanceVariant("C", 8, 8, 99)
                        ),
                        // Truth block histogram
                        Stream.of(new Object[][] {
                                { "10,99", 1}
                        }).collect(Collectors.toMap(data -> (String) data[0], data -> (Integer) data[1])),
                        // Eval block histogram
                        Stream.of(new Object[][] {
                                { "7,99", 1}
                        }).collect(Collectors.toMap(data -> (String) data[0], data -> (Integer) data[1])),
                        // Confidence concordance
                        Stream.of(new Object[][] {
                                { "99,99", 7}
                        }).collect(Collectors.toMap(data -> (String) data[0], data -> (Integer) data[1])),
                },
                // End with gap then var
                {
                        // Truth variants
                        Arrays.asList(
                                new TestReferenceBlockConcordanceVariant("<NON_REF>", 1, 10, 99)
                        ),
                        // Eval variants
                        Arrays.asList(
                                new TestReferenceBlockConcordanceVariant("<NON_REF>", 1, 7, 99),
                                new TestReferenceBlockConcordanceVariant("C", 10, 10, 99)
                        ),
                        // Truth block histogram
                        Stream.of(new Object[][] {
                                { "10,99", 1}
                        }).collect(Collectors.toMap(data -> (String) data[0], data -> (Integer) data[1])),
                        // Eval block histogram
                        Stream.of(new Object[][] {
                                { "7,99", 1}
                        }).collect(Collectors.toMap(data -> (String) data[0], data -> (Integer) data[1])),
                        // Confidence concordance
                        Stream.of(new Object[][] {
                                { "99,99", 7}
                        }).collect(Collectors.toMap(data -> (String) data[0], data -> (Integer) data[1])),
                },
        };
    }

    @Test(dataProvider = "provideSyntheticGVCFs")
    public void testSyntheticGVCFs(final List<TestReferenceBlockConcordanceVariant> truthVariants, final List<TestReferenceBlockConcordanceVariant> evalVariants, final Map<String, Integer> expectedTruthBlockHistogram, final Map<String, Integer> expectedEvalBlockHistogram, final Map<String, Integer> expectedConfidenceConcordance) throws Exception {
        Pair<File, File> inputFiles = writeTestGVCFs(truthVariants, evalVariants);

        final Path truthBlockHistogramFile = createTempPath("truth_block_histogram", ".tsv");
        final Path evalBlockHistogramFile = createTempPath("eval_block_histogram", ".tsv");
        final Path confidenceConcordanceHistogramFile = createTempPath("confidence_concordance_histogram", ".tsv");

        final String[] args = {
                "--" + AbstractConcordanceWalker.TRUTH_VARIANTS_LONG_NAME, inputFiles.getLeft().toString(),
                "--" + AbstractConcordanceWalker.EVAL_VARIANTS_LONG_NAME, inputFiles.getRight().toString(),
                "--" + ReferenceBlockConcordance.TRUTH_BLOCK_HISTOGRAM_LONG_NAME, truthBlockHistogramFile.toString(),
                "--" + ReferenceBlockConcordance.EVAL_BLOCK_HISTOGRAM_LONG_NAME, evalBlockHistogramFile.toString(),
                "--" + ReferenceBlockConcordance.CONFIDENCE_CONCORDANCE_HISTOGRAM_LONG_NAME, confidenceConcordanceHistogramFile.toString(),
        };
        runCommandLine(args);

        MetricsFile<?, String> truthBlockMetrics = new MetricsFile<>();
        truthBlockMetrics.read(new FileReader(truthBlockHistogramFile.toFile()));
        if (truthBlockMetrics.getNumHistograms() > 0) {
            Assert.assertEquals(truthBlockMetrics.getNumHistograms(), 1);
            Histogram<String> truthBlockHistogram = truthBlockMetrics.getHistogram();
            // Check value and remove entry from expected histogram...
            truthBlockHistogram.values().forEach(bin -> {
                Assert.assertTrue(expectedTruthBlockHistogram.containsKey(bin.getId()));
                Assert.assertEquals(bin.getValue(), (double) expectedTruthBlockHistogram.get(bin.getId()));
                expectedTruthBlockHistogram.remove(bin.getId());
            });
            // ... and make sure it is empty and all values have been visited
        }
        Assert.assertEquals(expectedTruthBlockHistogram.size(), 0);

        MetricsFile<?, String> evalBlockMetrics = new MetricsFile<>();
        evalBlockMetrics.read(new FileReader(evalBlockHistogramFile.toFile()));
        if (evalBlockMetrics.getNumHistograms() > 0) {
            Assert.assertEquals(evalBlockMetrics.getNumHistograms(), 1);
            Histogram<String> evalBlockHistogram = evalBlockMetrics.getHistogram();
            // Check value and remove entry from expected histogram...
            evalBlockHistogram.values().forEach(bin -> {
                Assert.assertTrue(expectedEvalBlockHistogram.containsKey(bin.getId()));
                Assert.assertEquals(bin.getValue(), (double) expectedEvalBlockHistogram.get(bin.getId()));
                expectedEvalBlockHistogram.remove(bin.getId());
            });
            // ... and make sure it is empty and all values have been visited
        }
        Assert.assertEquals(expectedEvalBlockHistogram.size(), 0);

        MetricsFile<?, String> confidenceConcordanceMetrics = new MetricsFile<>();
        confidenceConcordanceMetrics.read(new FileReader(confidenceConcordanceHistogramFile.toFile()));
        if (confidenceConcordanceMetrics.getNumHistograms() > 0) {
            Assert.assertEquals(confidenceConcordanceMetrics.getNumHistograms(), 1);
            Histogram<String> confidenceConcordanceHistogram = confidenceConcordanceMetrics.getHistogram();
            // Check value and remove entry from expected histogram...
            confidenceConcordanceHistogram.values().forEach(bin -> {
                Assert.assertTrue(expectedConfidenceConcordance.containsKey(bin.getId()));
                Assert.assertEquals(bin.getValue(), (double) expectedConfidenceConcordance.get(bin.getId()));
                expectedConfidenceConcordance.remove(bin.getId());
            });
            // ... and make sure it is empty and all values have been visited
        }
        Assert.assertEquals(expectedConfidenceConcordance.size(), 0);
    }

    @Test
    public void testDoesNotCrashWithNO_VARIATIONAlleles() {
        final File evalVcf = new File(CONCORDANCE_TEST_DIR, "noVariationAlleles.vcf");
        final File truthVcf = new File(CONCORDANCE_TEST_DIR, "noVariationAlleles.vcf");

        final String[] args = {
                "--" + AbstractConcordanceWalker.EVAL_VARIANTS_LONG_NAME, evalVcf.toString(),
                "--" + AbstractConcordanceWalker.TRUTH_VARIANTS_LONG_NAME, truthVcf.toString(),
                "--" + ReferenceBlockConcordance.TRUTH_BLOCK_HISTOGRAM_LONG_NAME, "/dev/null",
                "--" + ReferenceBlockConcordance.EVAL_BLOCK_HISTOGRAM_LONG_NAME, "/dev/null",
                "--" + ReferenceBlockConcordance.CONFIDENCE_CONCORDANCE_HISTOGRAM_LONG_NAME, "/dev/null",
        };

        runCommandLine(args);
    }

    static class TestReferenceBlockConcordanceVariant {
        private final String altAllele;
        private final int start;
        private final int stop;
        private final int confidence;

        public TestReferenceBlockConcordanceVariant(final String altAllele, final int start, final int stop, final int confidence) {
            this.altAllele = altAllele;
            this.start = start;
            this.stop = stop;
            this.confidence = confidence;
        }

        public String getAltAllele() { return altAllele; }
        public int getStart() { return start; }
        public int getStop() { return stop; }
        public int getConfidence() { return confidence; }
    }
}
