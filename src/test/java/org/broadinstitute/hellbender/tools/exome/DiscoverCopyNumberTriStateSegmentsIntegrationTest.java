package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.tools.exome.hmm.CopyNumberTriStateHiddenMarkovModelArgumentCollection;
import org.broadinstitute.hellbender.utils.IndexRange;
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
 * Integration tests for {@link DiscoverCopyNumberTriStateSegments}
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class DiscoverCopyNumberTriStateSegmentsIntegrationTest extends CopyNumberTriStateSegmentsCallerIntegrationTest {

    @Override
    public String getTestedClassName() {
        return DiscoverCopyNumberTriStateSegments.class.getSimpleName();
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

    private List<CopyNumberTriStateSegmentRecord> readOutputRecords(final File outputFile) {
        try (final CopyNumberTriStateSegmentRecordReader reader = new CopyNumberTriStateSegmentRecordReader(outputFile)) {
            return reader.stream().collect(Collectors.toList());
        } catch (final IOException ex) {
            Assert.fail("problems reading the output file " + outputFile);
            throw new RuntimeException(ex);
        }
    }

    //TODO: this test used to contain a tet of concordance with XHMM.  It no longer does that because our model has
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
        final List<CopyNumberTriStateSegmentRecord> outputRecords = readOutputRecords(outputFile);
        assertOutputIsInOrder(outputRecords, targets);
        assertOutputHasConsistentNumberOfTargets(outputRecords, targets);
        final Map<String, List<CopyNumberTriStateSegmentRecord>> outputRecordsBySample = splitOutputRecordBySample(outputRecords);
        assertSampleNames(outputRecordsBySample.keySet(), chain);
        for (final List<CopyNumberTriStateSegmentRecord> sampleRecords : outputRecordsBySample.values()) {
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
        arguments.add(String.valueOf(DiscoverCopyNumberTriStateSegments.ZScoreDimension.NONE.toString()));
        runCommandLine(arguments.toArray(new String[arguments.size()]));
    }


    private void assertSampleSegmentsCoordinates(List<CopyNumberTriStateSegmentRecord> sampleRecords, TargetCollection<Target> targets) {
        for (final CopyNumberTriStateSegmentRecord record : sampleRecords) {
            final IndexRange range = targets.indexRange(record.getSegment());
            Assert.assertTrue(range.size() > 0);
            Assert.assertEquals(record.getSegment().getContig(),targets.location(range.from).getContig());
            Assert.assertEquals(record.getSegment().getStart(), targets.location(range.from).getStart());
            Assert.assertEquals(record.getSegment().getEnd(), targets.location(range.to - 1).getEnd());
        }
    }

    private void assertSampleSegmentsCoverAllTargets(final List<CopyNumberTriStateSegmentRecord> sampleRecords, final TargetCollection<Target> targets) {
        int next = 0;
        for (final CopyNumberTriStateSegmentRecord record : sampleRecords) {
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

    private Map<String,List<CopyNumberTriStateSegmentRecord>> splitOutputRecordBySample(final List<CopyNumberTriStateSegmentRecord> outputRecords) {
            return outputRecords.stream().collect(Collectors.groupingBy(CopyNumberTriStateSegmentRecord::getSampleName));
    }

    private void assertOutputIsInOrder(final List<CopyNumberTriStateSegmentRecord> outputRecords, final TargetCollection<Target> targets) {
        for (int i = 1; i < outputRecords.size(); i++) {
            final CopyNumberTriStateSegmentRecord nextRecord = outputRecords.get(i);
            final CopyNumberTriStateSegmentRecord previousRecord = outputRecords.get(i - 1);
            final IndexRange nextRange = targets.indexRange(nextRecord.getSegment());
            final IndexRange previousRange = targets.indexRange(previousRecord.getSegment());
            Assert.assertTrue(nextRange.from >= previousRange.from);
        }
    }

    private void assertOutputHasConsistentNumberOfTargets(final List<CopyNumberTriStateSegmentRecord> outputRecords, final TargetCollection<Target> targets) {
        for (final CopyNumberTriStateSegmentRecord nextRecord : outputRecords) {
            final IndexRange indexRange = targets.indexRange(nextRecord.getSegment());
            Assert.assertEquals(indexRange.to - indexRange.from, nextRecord.getSegment().getTargetCount());
        }
    }

    @DataProvider(name = "testBadModelArgumentsData")
    public Object[][] testBadModelArgumentsData() {
        return new Object[][] {
                {CopyNumberTriStateHiddenMarkovModelArgumentCollection.EVENT_START_PROBABILITY_FULL_NAME, -1.0D},
                {CopyNumberTriStateHiddenMarkovModelArgumentCollection.MEAN_DELETION_COVERAGE_SHIFT_SHORT_NAME, 1.1D},
                {CopyNumberTriStateHiddenMarkovModelArgumentCollection.MEAN_DUPLICATION_COVERAGE_SHIFT_SHORT_NAME, -1.1D},
                {CopyNumberTriStateHiddenMarkovModelArgumentCollection.MEAN_EVENT_SIZE_SHORT_NAME, -1.0D},
        };
    }
}
