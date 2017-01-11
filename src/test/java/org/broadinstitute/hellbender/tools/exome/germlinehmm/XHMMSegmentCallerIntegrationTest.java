package org.broadinstitute.hellbender.tools.exome.germlinehmm;

import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.tools.exome.TargetArgumentCollection;
import org.broadinstitute.hellbender.tools.exome.TargetCollection;
import org.broadinstitute.hellbender.tools.exome.germlinehmm.xhmm.XHMMSegmentCaller;
import org.broadinstitute.hellbender.tools.exome.germlinehmm.xhmm.XHMMArgumentCollection;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.hmm.segmentation.HiddenStateSegmentRecord;
import org.broadinstitute.hellbender.utils.hmm.segmentation.HiddenStateSegmentRecordReader;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Integration tests for {@link XHMMSegmentCaller}
 *
 * @author Valentin Ruano-Rugitbio &lt;valentin@broadinstitute.org&gt;
 */
public class XHMMSegmentCallerIntegrationTest extends XHMMSegmentCallerBaseIntegrationTest {

    @Override
    public String getTestedClassName() {
        return XHMMSegmentCaller.class.getSimpleName();
    }

    @Test(dataProvider = "testBadModelArgumentsData", expectedExceptions = IllegalArgumentException.class)
    public void testBadModelArguments(final String argumentShortName, final double value) {
        final List<String> arguments = new ArrayList<>();
        arguments.add("-" + argumentShortName);
        arguments.add(String.valueOf(value));
        final File inputFile = createTempFile("input", ".tab");
        Assert.assertTrue(inputFile.delete());
        arguments.add("-" + StandardArgumentDefinitions.INPUT_SHORT_NAME);
        arguments.add(inputFile.getAbsolutePath());
        final File outputFile = createTempFile("output", ".tab");
        Assert.assertTrue(outputFile.delete());
        arguments.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        arguments.add(outputFile.getAbsolutePath());
        runCommandLine(arguments.toArray(new String[arguments.size()]));
    }

    private List<HiddenStateSegmentRecord<CopyNumberTriState, Target>> readOutputRecords(final File outputFile) {
        try (final HiddenStateSegmentRecordReader<CopyNumberTriState, Target> reader =
                     new HiddenStateSegmentRecordReader<>(outputFile, CopyNumberTriState::fromCallString)) {
            return reader.stream().collect(Collectors.toList());
        } catch (final IOException ex) {
            Assert.fail("problems reading the output file " + outputFile);
            throw new RuntimeException(ex);
        }
    }

    //TODO: this test used to contain a test of concordance with XHMM.  It no longer does that because our model has
    //TODO: diverged from XHMM's.  Eventually the right thing to do is use the simulateChain() method to generate
    //TODO: simulated data for some artificial set of CNV segments and to test concordance with those segments.
    //TODO: however we still use XHMM's emission model, which is both not generative and quite silly.  Once we
    //TODO: have a generative model of coverage we can modify simulateChain() accordingly and then write a concordance
    //TODO: test here.  Until then, we do not have an integration test but we do have our ongoing evaluations, which
    //TODO: show the superiority of our modifications versus the original XHMM model.
    @Test(dataProvider = "simulatedChainData")
    public File testRunCommandLine(final HiddenMarkovModelChain chain) {
        final File inputFile = writeChainInTempFile(chain);
        final File outputFile = createTempFile("output", ".tab");
        runCommandLine(chain, inputFile, outputFile);

        Assert.assertTrue(outputFile.exists());
        final TargetCollection<Target> targets = TargetArgumentCollection.readTargetCollection(REALISTIC_TARGETS_FILE);
        final List<HiddenStateSegmentRecord<CopyNumberTriState, Target>> outputRecords = readOutputRecords(outputFile);
        assertOutputIsInOrder(outputRecords, targets);
        assertOutputHasConsistentNumberOfTargets(outputRecords, targets);
        final Map<String, List<HiddenStateSegmentRecord<CopyNumberTriState, Target>>> outputRecordsBySample =
                splitOutputRecordBySample(outputRecords);
        assertSampleNames(outputRecordsBySample.keySet(), chain);
        for (final List<HiddenStateSegmentRecord<CopyNumberTriState, Target>> sampleRecords : outputRecordsBySample.values()) {
            assertSampleSegmentsCoverAllTargets(sampleRecords, targets);
            assertSampleSegmentsCoordinates(sampleRecords, targets);
        }
        return outputFile;
    }

    public void runCommandLine(final HiddenMarkovModelChain chain, final File inputFile, final File outputFile) {
        // The model arguments:
        final List<String> arguments = new ArrayList<>();
        arguments.add("-" + StandardArgumentDefinitions.INPUT_SHORT_NAME);
        arguments.add(inputFile.getAbsolutePath());
        arguments.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        arguments.add(outputFile.getAbsolutePath());
        loadModelArguments(chain, arguments);
        arguments.add(String.valueOf(XHMMSegmentCaller.ZScoreDimension.NONE.toString()));
        runCommandLine(arguments.toArray(new String[arguments.size()]));
    }


    private void assertSampleSegmentsCoordinates(List<HiddenStateSegmentRecord<CopyNumberTriState, Target>> sampleRecords,
                                                 TargetCollection<Target> targets) {
        for (final HiddenStateSegmentRecord<CopyNumberTriState, Target> record : sampleRecords) {
            final IndexRange range = targets.indexRange(record.getSegment());
            Assert.assertTrue(range.size() > 0);
            Assert.assertEquals(record.getSegment().getContig(),targets.location(range.from).getContig());
            Assert.assertEquals(record.getSegment().getStart(), targets.location(range.from).getStart());
            Assert.assertEquals(record.getSegment().getEnd(), targets.location(range.to - 1).getEnd());
        }
    }

    private void assertSampleSegmentsCoverAllTargets(final List<HiddenStateSegmentRecord<CopyNumberTriState, Target>> sampleRecords,
                                                     final TargetCollection<Target> targets) {
        int next = 0;
        for (final HiddenStateSegmentRecord<CopyNumberTriState, Target> record : sampleRecords) {
            final IndexRange range = targets.indexRange(record.getSegment());
            Assert.assertEquals(range.from, next);
            next = range.to;
        }
    }

    private void assertSampleNames(final Set<String> samples, final HiddenMarkovModelChain chain) {
        final int numberOfSamples = chain.data.size();
        final List<String> sampleNames = IntStream.range(0, chain.data.size()).mapToObj(a -> "SAMPLE_" + a).collect(Collectors.toList());
        Assert.assertEquals(samples.size(), numberOfSamples);
        for (final String sample : sampleNames) {
            Assert.assertTrue(samples.contains(sample));
        }
    }

    private Map<String,List<HiddenStateSegmentRecord<CopyNumberTriState, Target>>> splitOutputRecordBySample(final List<HiddenStateSegmentRecord<CopyNumberTriState, Target>> outputRecords) {
            return outputRecords.stream().collect(Collectors.groupingBy(HiddenStateSegmentRecord::getSampleName));
    }

    private void assertOutputIsInOrder(final List<HiddenStateSegmentRecord<CopyNumberTriState, Target>> outputRecords,
                                       final TargetCollection<Target> targets) {
        for (int i = 1; i < outputRecords.size(); i++) {
            final HiddenStateSegmentRecord<CopyNumberTriState, Target> nextRecord = outputRecords.get(i);
            final HiddenStateSegmentRecord<CopyNumberTriState, Target> previousRecord = outputRecords.get(i - 1);
            final IndexRange nextRange = targets.indexRange(nextRecord.getSegment());
            final IndexRange previousRange = targets.indexRange(previousRecord.getSegment());
            Assert.assertTrue(nextRange.from >= previousRange.from);
        }
    }

    private void assertOutputHasConsistentNumberOfTargets(final List<HiddenStateSegmentRecord<CopyNumberTriState, Target>> outputRecords,
                                                          final TargetCollection<Target> targets) {
        for (final HiddenStateSegmentRecord<CopyNumberTriState, Target> nextRecord : outputRecords) {
            final IndexRange indexRange = targets.indexRange(nextRecord.getSegment());
            Assert.assertEquals(indexRange.to - indexRange.from, nextRecord.getSegment().getTargetCount());
        }
    }

    @DataProvider(name = "testBadModelArgumentsData")
    public Object[][] testBadModelArgumentsData() {
        return new Object[][] {
                {XHMMArgumentCollection.EVENT_START_PROBABILITY_FULL_NAME, -1.0D},
                {XHMMArgumentCollection.MEAN_DELETION_COVERAGE_SHIFT_SHORT_NAME, 1.1D},
                {XHMMArgumentCollection.MEAN_DUPLICATION_COVERAGE_SHIFT_SHORT_NAME, -1.1D},
                {XHMMArgumentCollection.MEAN_EVENT_SIZE_SHORT_NAME, -1.0D},
        };
    }
}
