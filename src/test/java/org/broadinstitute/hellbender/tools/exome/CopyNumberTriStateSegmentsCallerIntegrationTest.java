package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.exome.germlinehmm.*;
import org.broadinstitute.hellbender.utils.hmm.CopyNumberTriState;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

/**
 * Parent class for integration testers for {@link CopyNumberTriStateSegmentCaller} subclasses.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public abstract class CopyNumberTriStateSegmentsCallerIntegrationTest extends CommandLineProgramTest {

    public static final File REALISTIC_TARGETS_FILE = new File("src/test/resources/org/broadinstitute/hellbender/tools/exome/germlinehmm/realistic-targets.tab");
    public static final TargetCollection<Target> REALISTIC_TARGETS;

    static {
        try (final TargetTableReader reader = new TargetTableReader(REALISTIC_TARGETS_FILE)) {
            REALISTIC_TARGETS = new HashedListTargetCollection<>(reader.stream().collect(Collectors.toList()));
        } catch (final IOException ex) {
            throw new GATKException("could not read the realistic-targets file", ex);
        }
    }
    public static final CopyNumberTriStateHiddenMarkovModel[] TEST_MODELS = new CopyNumberTriStateHiddenMarkovModel[] {
        new CopyNumberTriStateHiddenMarkovModel(1e-4, 70_000, -3, 3)
    };

    public static File writeChainInTempFile(final DiscoverCopyNumberTriStateSegmentsIntegrationTest.HiddenMarkovModelChain chain) {
        final File result = createTempFile("chain-",".tab");
        //final File result = new File("/tmp/input");
        final List<String> sampleNames = IntStream.range(0, chain.data.size()).mapToObj(a -> "SAMPLE_" + a).collect(Collectors.toList());
        final List<String> columnNames = new ArrayList<>(sampleNames.size() + 4);
        columnNames.addAll(Arrays.asList(TargetTableColumn.CONTIG.toString(), TargetTableColumn.START.toString(),
                TargetTableColumn.END.toString(), TargetTableColumn.NAME.toString()));
        columnNames.addAll(sampleNames);
        try (final TableWriter<Integer> writer = new TableWriter<Integer>(result,
                new TableColumnCollection(columnNames)) {
            @Override
            protected void composeLine(final Integer record, final DataLine dataLine) {
                dataLine.append(chain.targets.get(record).getContig())
                        .append(chain.targets.get(record).getStart())
                        .append(chain.targets.get(record).getEnd())
                        .append(chain.targets.get(record).getName());
                for (final List<Double> sampleData : chain.data) {
                    dataLine.append(sampleData.get(record));
                }
            }}) {

            writer.writeAllRecords(IntStream.range(0, chain.targets.size()).boxed().collect(Collectors.toList()));
        } catch (final IOException ex) {
            throw new UncheckedIOException(ex);
        }
        return result;
    }

    @DataProvider(name = "simulatedChainData")
    public static Object[][] simulateChainDataProvider() {
        final Random rdn = new Random(131313);
        final List<Object[]> result = simulateChainData(rdn);
        return result.toArray(new Object[result.size()][]);
    }

    protected static List<Object[]> simulateChainData(Random rdn) {
        return Stream.of(TEST_MODELS)
                .map(m -> new Object[] { simulateChain(m, REALISTIC_TARGETS, 3, rdn) })
                .collect(Collectors.toList());
    }

    private static HiddenMarkovModelChain simulateChain(final CopyNumberTriStateHiddenMarkovModel model, final TargetCollection<Target> targetCollection, final int sampleCount, final Random rdn) {
        final List<Target> targets = targetCollection.targets();
        final List<CopyNumberTriState> truth = new ArrayList<>(targets.size());
        final List<List<Double>> data = new ArrayList<>(sampleCount);
        for (int i = 0; i < sampleCount; i++) {
            final List<Double> sampleData = randomSampleData(model, rdn, targets, truth);
            data.add(sampleData);
        }
        return new HiddenMarkovModelChain(model, targets, truth, data);
    }

    private static List<Double> randomSampleData(CopyNumberTriStateHiddenMarkovModel model, Random rdn, List<Target> targets, List<CopyNumberTriState> truth) {
        truth.addAll(model.generateHiddenStateChain(targets));
        return truth.stream().map(state -> model.randomDatum(state, rdn)).collect(Collectors.toList());
    }


    // This is a meta-test just to double check that simulateChain does provide a result.
    // Otherwise this would result in a silent skip of the main test of this class.
    @Test
    public void testSimulateChainDataProvider() {
        final Random rdn = new Random(131313);
        simulateChainData(rdn);
    }

    protected void loadModelArguments(HiddenMarkovModelChain chain, List<String> arguments) {
        arguments.add("-" + CopyNumberTriStateHiddenMarkovModelArgumentCollection.MEAN_DELETION_COVERAGE_SHIFT_SHORT_NAME);
        arguments.add(String.valueOf(chain.model.getDeletionMean()));
        arguments.add("-" + CopyNumberTriStateHiddenMarkovModelArgumentCollection.MEAN_DUPLICATION_COVERAGE_SHIFT_SHORT_NAME);
        arguments.add(String.valueOf(chain.model.getDuplicationMean()));
        arguments.add("-" + CopyNumberTriStateHiddenMarkovModelArgumentCollection.EVENT_START_PROBABILITY_SHORT_NAME);
        arguments.add(String.valueOf(chain.model.getEventStartProbability()));
        arguments.add("-" + CopyNumberTriStateHiddenMarkovModelArgumentCollection.MEAN_EVENT_SIZE_SHORT_NAME);
        arguments.add(String.valueOf(chain.model.getMeanEventSize()));
        arguments.add("-" + DiscoverCopyNumberTriStateSegments.ZSCORE_DIMENSION_SHORT_NAME);
    }

    public static class HiddenMarkovModelChain {

        public final List<Target> targets;
        public final List<CopyNumberTriState> truth;
        public final List<List<Double>> data;
        public final CopyNumberTriStateHiddenMarkovModel model;

        protected HiddenMarkovModelChain(final CopyNumberTriStateHiddenMarkovModel model, final List<Target> targets, final List<CopyNumberTriState> truth, final List<List<Double>> data) {
            this.targets = targets;
            this.truth = truth;
            this.data = data;
            this.model = model;
        }
    }
}
