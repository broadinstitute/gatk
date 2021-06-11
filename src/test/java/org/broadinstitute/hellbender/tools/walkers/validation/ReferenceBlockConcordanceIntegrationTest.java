package org.broadinstitute.hellbender.tools.walkers.validation;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.Histogram;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.engine.AbstractConcordanceWalker;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.hellbender.utils.variant.writers.GVCFWriter;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.sam.util.Pair;

import java.io.File;
import java.io.FileReader;
import java.nio.file.Path;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

/**
 * Created by Michael Gatzen on 9/9/20.
 */
public class ReferenceBlockConcordanceIntegrationTest extends CommandLineProgramTest{

    private static final String CONCORDANCE_TEST_DIR = toolsTestDir + "concordance/";
    private static final String HAPLOTYPECALLER_TEST_DIR = toolsTestDir + "haplotypecaller/";

    @Test
    public void testIdentical() throws Exception {
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

        final MetricsFile<?, String> truthBlockMetrics = new MetricsFile<>();
        truthBlockMetrics.read(new FileReader(truthBlockHistogramFile.toFile()));
        final MetricsFile<?, String> evalBlockMetrics = new MetricsFile<>();
        evalBlockMetrics.read(new FileReader(evalBlockHistogramFile.toFile()));

        Assert.assertEquals(truthBlockMetrics.getNumHistograms(), 1);
        Assert.assertEquals(evalBlockMetrics.getNumHistograms(), 1);

        final Histogram<String> truthBlockHistogram = truthBlockMetrics.getHistogram();
        final Histogram<String> evalBlockHistogram = evalBlockMetrics.getHistogram();

        // Got this number by counting the <NON_REF> alt alleles in the test GVCF file
        Assert.assertEquals(truthBlockHistogram.getSumOfValues(), 1038);
        Assert.assertEquals(evalBlockHistogram.getSumOfValues(), 1038);

        Assert.assertEquals(truthBlockHistogram, evalBlockHistogram);

        // For confidence concordance, check that there are only values on the diagonal
        final MetricsFile<?, String> confidenceConcordanceMetrics = new MetricsFile<>();
        confidenceConcordanceMetrics.read(new FileReader(confidenceConcordanceHistogramFile.toFile()));
        Assert.assertEquals(confidenceConcordanceMetrics.getNumHistograms(), 1);
        final Histogram<String> confidenceConcordanceHistogram = confidenceConcordanceMetrics.getHistogram();
        confidenceConcordanceHistogram.values().forEach(bin -> {
            final String[] confidenceValues = bin.getId().split(",");
            Assert.assertEquals(confidenceValues[0], confidenceValues[1]);
        });
    }

    @Test
    public void testMultipleContigs() throws Exception {
        final File truthVcf = new File(CONCORDANCE_TEST_DIR, "multiple-contigs-truth.gvcf");
        final File evalVcf = new File(CONCORDANCE_TEST_DIR, "multiple-contigs-eval.gvcf");
        final Path truthBlockHistogramFile = createTempPath("truth_block_histogram", ".tsv");
        final Path evalBlockHistogramFile = createTempPath("eval_block_histogram", ".tsv");
        final Path confidenceConcordanceHistogramFile = createTempPath("confidence_concordance_histogram", ".tsv");

        final Histogram<String> expectedTruthBlockHistogram = new Histogram<>();
        final Histogram<String> expectedEvalBlockHistogram = new Histogram<>();
        final Histogram<String> expectedConfidenceConcordanceHistogram = new Histogram<>();

        // Expected values
        // Truth block histogram
        Stream.of(new Object[][] {
                { "1,99", 1},
                { "1,98", 1},
                { "2,99", 1},
                { "5,99", 3},
        }).collect(Collectors.toMap(data -> (String) data[0], data -> (Integer) data[1]))
          .forEach(expectedTruthBlockHistogram::increment);
        // Eval block histogram
        Stream.of(new Object[][] {
                { "2,99", 2},
                { "3,99", 1},
                { "1,99", 1},
        }).collect(Collectors.toMap(data -> (String) data[0], data -> (Integer) data[1]))
          .forEach(expectedEvalBlockHistogram::increment);
        // Confidence concordance
        Stream.of(new Object[][] {
                { "99,99", 5},
                { "98,99", 1},
        }).collect(Collectors.toMap(data -> (String) data[0], data -> (Integer) data[1]))
          .forEach(expectedConfidenceConcordanceHistogram::increment);

        final String[] args = {
                "--" + AbstractConcordanceWalker.EVAL_VARIANTS_LONG_NAME, evalVcf.toString(),
                "--" + AbstractConcordanceWalker.TRUTH_VARIANTS_LONG_NAME, truthVcf.toString(),
                "--" + ReferenceBlockConcordance.TRUTH_BLOCK_HISTOGRAM_LONG_NAME, truthBlockHistogramFile.toString(),
                "--" + ReferenceBlockConcordance.EVAL_BLOCK_HISTOGRAM_LONG_NAME, evalBlockHistogramFile.toString(),
                "--" + ReferenceBlockConcordance.CONFIDENCE_CONCORDANCE_HISTOGRAM_LONG_NAME, confidenceConcordanceHistogramFile.toString(),
        };
        runCommandLine(args);

        final MetricsFile<?, String> truthBlockMetrics = new MetricsFile<>();
        truthBlockMetrics.read(new FileReader(truthBlockHistogramFile.toFile()));
        Assert.assertEquals(truthBlockMetrics.getNumHistograms(), 1);
        final Histogram<String> truthBlockHistogram = truthBlockMetrics.getHistogram();
        Assert.assertEquals(truthBlockHistogram, expectedTruthBlockHistogram);

        final MetricsFile<?, String> evalBlockMetrics = new MetricsFile<>();
        evalBlockMetrics.read(new FileReader(evalBlockHistogramFile.toFile()));
        Assert.assertEquals(evalBlockMetrics.getNumHistograms(), 1);
        final Histogram<String> evalBlockHistogram = evalBlockMetrics.getHistogram();
        Assert.assertEquals(evalBlockHistogram, expectedEvalBlockHistogram);

        final MetricsFile<?, String> confidenceConcordanceMetrics = new MetricsFile<>();
        confidenceConcordanceMetrics.read(new FileReader(confidenceConcordanceHistogramFile.toFile()));
        Assert.assertEquals(confidenceConcordanceMetrics.getNumHistograms(), 1);
        final Histogram<String> confidenceConcordanceHistogram = confidenceConcordanceMetrics.getHistogram();
        Assert.assertEquals(confidenceConcordanceHistogram, expectedConfidenceConcordanceHistogram);
    }

    private Pair<File, File> writeTestGVCFs(List<VariantContext> truthVariants, List<VariantContext> evalVariants) {
        final File truthFile = createTempFile("truth", ".gvcf");
        final GVCFWriter truthWriter = new GVCFWriter(
                GATKVariantContextUtils.createVCFWriter(truthFile.toPath(), null, false, Options.ALLOW_MISSING_FIELDS_IN_HEADER),
                IntStream.range(1, 100).boxed().collect(Collectors.toList()),
                true
                );
        final VCFHeader header = new VCFHeader(new HashSet<>(), Collections.singletonList("TESTSAMPLE"));
        header.setSequenceDictionary(new SAMSequenceDictionary(Collections.singletonList(new SAMSequenceRecord("test_contig", 1000))));
        truthWriter.writeHeader(header);

        truthVariants.forEach(truthWriter::add);
        truthWriter.close();

        final File evalFile = createTempFile("eval", ".gvcf");
        final GVCFWriter evalWriter = new GVCFWriter(
                GATKVariantContextUtils.createVCFWriter(evalFile.toPath(), null, false, Options.ALLOW_MISSING_FIELDS_IN_HEADER),
                IntStream.range(1, 100).boxed().collect(Collectors.toList()),
                true
        );
        evalWriter.writeHeader(header);

        evalVariants.forEach(evalWriter::add);
        evalWriter.close();

        return new Pair<>(truthFile, evalFile);
    }

    private static VariantContext constructTestVariantContext(final Allele altAllele, final int start, final int stop, final int confidence) {
        final Allele refAllele = Allele.REF_A;
        final List<Allele> alleles = Arrays.asList(refAllele, altAllele);
        final VariantContextBuilder variantContextBuilder = new VariantContextBuilder("TEST", "test_contig", start, stop, alleles);
        final GenotypeBuilder genotypeBuilder = new GenotypeBuilder("TESTSAMPLE", altAllele.isNonRefAllele() ? GATKVariantContextUtils.homozygousAlleleList(refAllele, 2) : alleles);
        if (altAllele.isNonRefAllele()) {
            genotypeBuilder.GQ(confidence);
            genotypeBuilder.PL(new int[] { 0, 0, 0 });
            variantContextBuilder.attribute(VCFConstants.END_KEY, stop);
        }
        final Genotype gt = genotypeBuilder.DP(30).make();
        variantContextBuilder.genotypes(gt);
        return variantContextBuilder.make();
    }

    @DataProvider
    public Object[][] provideSyntheticGVCFs() {
        return new Object[][] {
                // No non_ref blocks
                {
                        // Truth variants
                        Arrays.asList(
                                constructTestVariantContext(Allele.ALT_C, 1, 1, 99)
                        ),
                        // Eval variants
                        Arrays.asList(
                                constructTestVariantContext(Allele.NON_REF_ALLELE, 1, 1, 99)
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
                                constructTestVariantContext(Allele.NON_REF_ALLELE, 1, 1, 99)
                        ),
                        // Eval variants
                        Arrays.asList(
                                constructTestVariantContext(Allele.NON_REF_ALLELE, 1, 1, 99)
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
                                constructTestVariantContext(Allele.NON_REF_ALLELE, 1, 10, 99)
                        ),
                        // Eval variants
                        Arrays.asList(
                                constructTestVariantContext(Allele.ALT_C, 5, 5, 99),
                                constructTestVariantContext(Allele.NON_REF_ALLELE, 6, 10, 99)
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
                                constructTestVariantContext(Allele.NON_REF_ALLELE, 1, 10, 99)
                        ),
                        // Eval variants
                        Arrays.asList(
                                constructTestVariantContext(Allele.NON_REF_ALLELE, 6, 10, 99)
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
                                constructTestVariantContext(Allele.NON_REF_ALLELE, 11, 20, 99)
                        ),
                        // Eval variants
                        Arrays.asList(
                                constructTestVariantContext(Allele.NON_REF_ALLELE, 3, 6, 98),
                                constructTestVariantContext(Allele.NON_REF_ALLELE, 11, 20, 99)
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
                                constructTestVariantContext(Allele.NON_REF_ALLELE, 11, 20, 99)
                        ),
                        // Eval variants
                        Arrays.asList(
                                constructTestVariantContext(Allele.NON_REF_ALLELE, 3, 6, 98)
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
                                constructTestVariantContext(Allele.NON_REF_ALLELE, 1, 10, 99)
                        ),
                        // Eval variants
                        Arrays.asList(
                                constructTestVariantContext(Allele.NON_REF_ALLELE, 1, 9, 99),
                                constructTestVariantContext(Allele.ALT_C, 10, 10, 99)
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
                                constructTestVariantContext(Allele.NON_REF_ALLELE, 1, 10, 99)
                        ),
                        // Eval variants
                        Arrays.asList(
                                constructTestVariantContext(Allele.NON_REF_ALLELE, 1, 7, 99),
                                constructTestVariantContext(Allele.ALT_C, 8, 8, 99)
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
                                constructTestVariantContext(Allele.NON_REF_ALLELE, 1, 10, 99)
                        ),
                        // Eval variants
                        Arrays.asList(
                                constructTestVariantContext(Allele.NON_REF_ALLELE, 1, 7, 99),
                                constructTestVariantContext(Allele.ALT_C, 10, 10, 99)
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
                // Multiple blocks vs one block
                {
                        // Truth variants
                        Arrays.asList(
                                constructTestVariantContext(Allele.NON_REF_ALLELE, 1, 10, 99)
                        ),
                        // Eval variants
                        Arrays.asList(
                                constructTestVariantContext(Allele.NON_REF_ALLELE, 1, 3, 99),
                                constructTestVariantContext(Allele.NON_REF_ALLELE, 4, 6, 98),
                                constructTestVariantContext(Allele.NON_REF_ALLELE, 7, 10, 99)
                        ),
                        // Truth block histogram
                        Stream.of(new Object[][] {
                                { "10,99", 1}
                        }).collect(Collectors.toMap(data -> (String) data[0], data -> (Integer) data[1])),
                        // Eval block histogram
                        Stream.of(new Object[][] {
                                { "3,99", 1},
                                { "3,98", 1},
                                { "4,99", 1}
                        }).collect(Collectors.toMap(data -> (String) data[0], data -> (Integer) data[1])),
                        // Confidence concordance
                        Stream.of(new Object[][] {
                                { "99,99", 7},
                                { "99,98", 3}
                        }).collect(Collectors.toMap(data -> (String) data[0], data -> (Integer) data[1])),
                },
                // Multiple blocks vs one block (ends coinciding)
                {
                        // Truth variants
                        Arrays.asList(
                                constructTestVariantContext(Allele.NON_REF_ALLELE, 1, 3, 99),
                                constructTestVariantContext(Allele.NON_REF_ALLELE, 7, 10, 98)
                        ),
                        // Eval variants
                        Arrays.asList(
                                constructTestVariantContext(Allele.NON_REF_ALLELE, 1, 3, 99),
                                constructTestVariantContext(Allele.NON_REF_ALLELE, 4, 6, 98),
                                constructTestVariantContext(Allele.NON_REF_ALLELE, 7, 10, 99)
                        ),
                        // Truth block histogram
                        Stream.of(new Object[][] {
                                { "3,99", 1},
                                { "4,98", 1}
                        }).collect(Collectors.toMap(data -> (String) data[0], data -> (Integer) data[1])),
                        // Eval block histogram
                        Stream.of(new Object[][] {
                                { "3,99", 1},
                                { "3,98", 1},
                                { "4,99", 1}
                        }).collect(Collectors.toMap(data -> (String) data[0], data -> (Integer) data[1])),
                        // Confidence concordance
                        Stream.of(new Object[][] {
                                { "99,99", 3},
                                { "98,99", 4}
                        }).collect(Collectors.toMap(data -> (String) data[0], data -> (Integer) data[1])),
                },
                // Multiple blocks vs one block (ends overlapping)
                {
                        // Truth variants
                        Arrays.asList(
                                constructTestVariantContext(Allele.NON_REF_ALLELE, 1, 4, 99),
                                constructTestVariantContext(Allele.NON_REF_ALLELE, 7, 10, 98)
                        ),
                        // Eval variants
                        Arrays.asList(
                                constructTestVariantContext(Allele.NON_REF_ALLELE, 1, 3, 99),
                                constructTestVariantContext(Allele.NON_REF_ALLELE, 4, 6, 98),
                                constructTestVariantContext(Allele.NON_REF_ALLELE, 7, 10, 99)
                        ),
                        // Truth block histogram
                        Stream.of(new Object[][] {
                                { "4,99", 1},
                                { "4,98", 1}
                        }).collect(Collectors.toMap(data -> (String) data[0], data -> (Integer) data[1])),
                        // Eval block histogram
                        Stream.of(new Object[][] {
                                { "3,99", 1},
                                { "3,98", 1},
                                { "4,99", 1}
                        }).collect(Collectors.toMap(data -> (String) data[0], data -> (Integer) data[1])),
                        // Confidence concordance
                        Stream.of(new Object[][] {
                                { "99,99", 3},
                                { "99,98", 1},
                                { "98,99", 4}
                        }).collect(Collectors.toMap(data -> (String) data[0], data -> (Integer) data[1])),
                },
                // Multiple blocks vs one block (three 1bp blocks)
                {
                        // Truth variants
                        Arrays.asList(
                                constructTestVariantContext(Allele.NON_REF_ALLELE, 1, 3, 99)
                        ),
                        // Eval variants
                        Arrays.asList(
                                constructTestVariantContext(Allele.NON_REF_ALLELE, 1, 1, 99),
                                constructTestVariantContext(Allele.NON_REF_ALLELE, 2, 2, 98),
                                constructTestVariantContext(Allele.NON_REF_ALLELE, 3, 3, 97)
                        ),
                        // Truth block histogram
                        Stream.of(new Object[][] {
                                { "3,99", 1}
                        }).collect(Collectors.toMap(data -> (String) data[0], data -> (Integer) data[1])),
                        // Eval block histogram
                        Stream.of(new Object[][] {
                                { "1,99", 1},
                                { "1,98", 1},
                                { "1,97", 1}
                        }).collect(Collectors.toMap(data -> (String) data[0], data -> (Integer) data[1])),
                        // Confidence concordance
                        Stream.of(new Object[][] {
                                { "99,99", 1},
                                { "99,98", 1},
                                { "99,97", 1}
                        }).collect(Collectors.toMap(data -> (String) data[0], data -> (Integer) data[1])),
                },
        };
    }

    @Test(dataProvider = "provideSyntheticGVCFs")
    public void testSyntheticGVCFs(final List<VariantContext> truthVariants, final List<VariantContext> evalVariants, final Map<String, Integer> expectedTruthBlocks, final Map<String, Integer> expectedEvalBlocks, final Map<String, Integer> expectedConfidenceConcordance) throws Exception {
        final Pair<File, File> inputFiles = writeTestGVCFs(truthVariants, evalVariants);

        final Path truthBlockHistogramFile = createTempPath("truth_block_histogram", ".tsv");
        final Path evalBlockHistogramFile = createTempPath("eval_block_histogram", ".tsv");
        final Path confidenceConcordanceHistogramFile = createTempPath("confidence_concordance_histogram", ".tsv");

        final Histogram<String> expectedTruthBlockHistogram = new Histogram<>();
        expectedTruthBlocks.forEach(expectedTruthBlockHistogram::increment);

        final Histogram<String> expectedEvalBlockHistogram = new Histogram<>();
        expectedEvalBlocks.forEach(expectedEvalBlockHistogram::increment);

        final Histogram<String> expectedConfidenceConcordanceHistogram = new Histogram<>();
        expectedConfidenceConcordance.forEach(expectedConfidenceConcordanceHistogram::increment);

        final String[] args = {
                "--" + AbstractConcordanceWalker.TRUTH_VARIANTS_LONG_NAME, inputFiles.getLeft().toString(),
                "--" + AbstractConcordanceWalker.EVAL_VARIANTS_LONG_NAME, inputFiles.getRight().toString(),
                "--" + ReferenceBlockConcordance.TRUTH_BLOCK_HISTOGRAM_LONG_NAME, truthBlockHistogramFile.toString(),
                "--" + ReferenceBlockConcordance.EVAL_BLOCK_HISTOGRAM_LONG_NAME, evalBlockHistogramFile.toString(),
                "--" + ReferenceBlockConcordance.CONFIDENCE_CONCORDANCE_HISTOGRAM_LONG_NAME, confidenceConcordanceHistogramFile.toString(),
        };
        runCommandLine(args);

        final MetricsFile<?, String> truthBlockMetrics = new MetricsFile<>();
        truthBlockMetrics.read(new FileReader(truthBlockHistogramFile.toFile()));
        if (expectedTruthBlocks.size() == 0) {
            Assert.assertEquals(truthBlockMetrics.getNumHistograms(), 0);
        } else {
            Assert.assertEquals(truthBlockMetrics.getNumHistograms(), 1);
            final Histogram<String> truthBlockHistogram = truthBlockMetrics.getHistogram();
            Assert.assertEquals(truthBlockHistogram, expectedTruthBlockHistogram);
        }

        final MetricsFile<?, String> evalBlockMetrics = new MetricsFile<>();
        evalBlockMetrics.read(new FileReader(evalBlockHistogramFile.toFile()));
        if (expectedEvalBlocks.size() == 0) {
            Assert.assertEquals(evalBlockMetrics.getNumHistograms(), 0);
        } else {
            Assert.assertEquals(evalBlockMetrics.getNumHistograms(), 1);
            final Histogram<String> evalBlockHistogram = evalBlockMetrics.getHistogram();
            Assert.assertEquals(evalBlockHistogram, expectedEvalBlockHistogram);
        }

        final MetricsFile<?, String> confidenceConcordanceMetrics = new MetricsFile<>();
        confidenceConcordanceMetrics.read(new FileReader(confidenceConcordanceHistogramFile.toFile()));
        if (expectedConfidenceConcordance.size() == 0) {
            Assert.assertEquals(confidenceConcordanceMetrics.getNumHistograms(), 0);
        } else {
            Assert.assertEquals(confidenceConcordanceMetrics.getNumHistograms(), 1);
            final Histogram<String> confidenceConcordanceHistogram = confidenceConcordanceMetrics.getHistogram();
            Assert.assertEquals(confidenceConcordanceHistogram, expectedConfidenceConcordanceHistogram);
        }
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
}
