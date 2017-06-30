package org.broadinstitute.hellbender.tools.exome.eval;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.vcf.VCFFileReader;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.tools.exome.TargetArgumentCollection;
import org.broadinstitute.hellbender.tools.exome.TargetCollection;
import org.broadinstitute.hellbender.utils.hmm.segmentation.HiddenStateSegment;
import org.broadinstitute.hellbender.utils.hmm.segmentation.HiddenStateSegmentRecord;
import org.broadinstitute.hellbender.utils.hmm.segmentation.HiddenStateSegmentRecordReader;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.tools.exome.germlinehmm.CopyNumberTriState;
import org.broadinstitute.hellbender.utils.test.TestResources;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Integration tests for {@link ConvertGSVariantsToSegments} class.
 */
public class ConvertGSVariantsToSegmentsIntegrationTest extends CommandLineProgramTest {

    private static final File LARGE_CNV_TEST_FILE_DIR = new File(TestResources.largeFileTestDir, "cnv");

    private static final File TEST_INPUT_FILE =
            new File("src/test/resources/org/broadinstitute/hellbender/tools/exome/eval", "gs-calls.vcf.gz");

    private static final File TEST_TARGET_FILE =
            new File(LARGE_CNV_TEST_FILE_DIR, "realistic-targets.tab");

    @Override
    public String getTestedClassName() {
        return ConvertGSVariantsToSegments.class.getSimpleName();
    }

    @Test
    public void testRunWithoutFrequencies() throws IOException {
        final File outputFile = createTempFile("output", ".segments");
        Assert.assertTrue(outputFile.delete());
        final List<String> arguments = new ArrayList<>();
        arguments.add("-" + StandardArgumentDefinitions.VARIANT_SHORT_NAME);
        arguments.add(TEST_INPUT_FILE.getAbsolutePath());
        arguments.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        arguments.add(outputFile.getAbsolutePath());
        arguments.add("-" + TargetArgumentCollection.TARGET_FILE_SHORT_NAME);
        arguments.add(TEST_TARGET_FILE.getAbsolutePath());
        runCommandLine(arguments);
        Assert.assertTrue(outputFile.exists());
        checkConsistency(TEST_INPUT_FILE, outputFile);
    }

    private void checkConsistency(final File vcf, final File segments) throws IOException {
        final List<HiddenStateSegmentRecord<CopyNumberTriState, Target>> expectedSegments = composeExpectedSegments(vcf, TargetArgumentCollection.readTargetCollection(TEST_TARGET_FILE));
        final List<HiddenStateSegmentRecord<CopyNumberTriState, Target>> observedSegments =
                new HiddenStateSegmentRecordReader<>(segments, CopyNumberTriState::fromCallString).stream().collect(Collectors.toList());
        Assert.assertEquals(observedSegments.size(), expectedSegments.size());
        for (int i = 0; i < observedSegments.size(); i++) {
            final HiddenStateSegmentRecord<CopyNumberTriState, Target> observedRecord = observedSegments.get(0);
            final HiddenStateSegmentRecord<CopyNumberTriState, Target> expectedRecord = expectedSegments.get(0);
            Assert.assertEquals(observedRecord.getSampleName(), expectedRecord.getSampleName());
            Assert.assertEquals(observedRecord.getSegment().getCall(), expectedRecord.getSegment().getCall());
            Assert.assertEquals(observedRecord.getSegment().getInterval(), expectedRecord.getSegment().getInterval());
            Assert.assertEquals(observedRecord.getSegment().getTargetCount(), expectedRecord.getSegment().getTargetCount());
            Assert.assertEquals(observedRecord.getSegment().getMean(), expectedRecord.getSegment().getMean(), 0.001);
            Assert.assertEquals(observedRecord.getSegment().getEventQuality(), expectedRecord.getSegment().getEventQuality(), 0.001);
            Assert.assertEquals(observedRecord.getSegment().getExactQuality(), expectedRecord.getSegment().getExactQuality(), 0.001, observedRecord.toString() + " vs " + expectedRecord.toString());
            Assert.assertEquals(observedRecord.getSegment().getStdev(), 0.0);
            Assert.assertTrue(Double.isNaN(observedRecord.getSegment().getEndQuality()));
            Assert.assertTrue(Double.isNaN(observedRecord.getSegment().getStartQuality()));
            Assert.assertTrue(Double.isNaN(observedRecord.getSegment().getSomeQuality()));
        }
    }

    private CopyNumberTriState expectedCall(final int cn) {
        if (cn < ConvertGSVariantsToSegments.NEUTRAL_COPY_NUMBER_DEFAULT) {
            return CopyNumberTriState.DELETION;
        } else if (cn == ConvertGSVariantsToSegments.NEUTRAL_COPY_NUMBER_DEFAULT) {
            return CopyNumberTriState.NEUTRAL;
        } else {
            return CopyNumberTriState.DUPLICATION;
        }
    }

    private List<HiddenStateSegmentRecord<CopyNumberTriState, Target>> composeExpectedSegments(final File vcf, final TargetCollection<Target> targets) throws IOException {
        final VCFFileReader reader = new VCFFileReader(vcf, false);
        final List<HiddenStateSegmentRecord<CopyNumberTriState, Target>> result = new ArrayList<>();
        reader.iterator().forEachRemaining(vc -> {
            final int targetCount = targets.indexRange(vc).size();
            for (final Genotype genotype : vc.getGenotypes()) {
                final int cn = Integer.parseInt(genotype.getExtendedAttribute("CN").toString());
                final double[] cnp = Stream.of(genotype.getExtendedAttribute("CNP").toString().replaceAll("\\[\\]", "").split(","))
                        .mapToDouble(Double::parseDouble).toArray();

                final double cnpSum = MathUtils.approximateLog10SumLog10(cnp);
                final CopyNumberTriState call = expectedCall(cn);
                final double exactLog10Prob = expectedExactLog10(call, cnp);
                final HiddenStateSegment<CopyNumberTriState, Target> expectedSegment = new HiddenStateSegment<>(
                        new SimpleInterval(vc), targetCount, Double.parseDouble(genotype.getExtendedAttribute("CNF").toString()),
                        0.000, call, -10.0 * exactLog10Prob, Double.NaN,  Double.NaN, Double.NaN,
                          -10.0 * (cnp[ConvertGSVariantsToSegments.NEUTRAL_COPY_NUMBER_DEFAULT] - cnpSum)
                );
                result.add(new HiddenStateSegmentRecord<>(genotype.getSampleName(), expectedSegment));
            }
        });
        return result;
    }

    private double expectedExactLog10(final CopyNumberTriState call, final double[] cnp) {
        if (call == CopyNumberTriState.DELETION) {
            return MathUtils.log10SumLog10(cnp, ConvertGSVariantsToSegments.NEUTRAL_COPY_NUMBER_DEFAULT, cnp.length) - MathUtils.log10SumLog10(cnp);
        } else if (call == CopyNumberTriState.NEUTRAL) {
            return MathUtils.approximateLog10SumLog10(
                    MathUtils.log10SumLog10(cnp, 0, ConvertGSVariantsToSegments.NEUTRAL_COPY_NUMBER_DEFAULT),
                    cnp.length <= ConvertGSVariantsToSegments.NEUTRAL_COPY_NUMBER_DEFAULT + 1 ? Double.NEGATIVE_INFINITY :
                            MathUtils.log10SumLog10(cnp, ConvertGSVariantsToSegments.NEUTRAL_COPY_NUMBER_DEFAULT + 1, cnp.length)
            ) - MathUtils.log10SumLog10(cnp);
        } else {
            return MathUtils.log10SumLog10(cnp, 0, ConvertGSVariantsToSegments.NEUTRAL_COPY_NUMBER_DEFAULT + 1) - MathUtils.log10SumLog10(cnp);
        }
    }
}
