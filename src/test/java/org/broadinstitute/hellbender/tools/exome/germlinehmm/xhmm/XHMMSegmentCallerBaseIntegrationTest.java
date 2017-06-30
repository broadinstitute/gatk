package org.broadinstitute.hellbender.tools.exome.germlinehmm.xhmm;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.exome.*;
import org.broadinstitute.hellbender.tools.exome.germlinehmm.CopyNumberTriState;
import org.broadinstitute.hellbender.utils.test.TestResources;
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
 * Parent class for integration testers for {@link XHMMSegmentCallerBase} subclasses.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public abstract class XHMMSegmentCallerBaseIntegrationTest extends CommandLineProgramTest {

    private static final File LARGE_CNV_TEST_FILE_DIR = new File(TestResources.largeFileTestDir, "cnv");

    public static final File TRUNCATED_REALISTIC_TARGETS_FILE = new File(LARGE_CNV_TEST_FILE_DIR, "truncated-realistic-targets.tab");
    public static final TargetCollection<Target> REALISTIC_TARGETS;

    static {
        try (final TargetTableReader reader = new TargetTableReader(TRUNCATED_REALISTIC_TARGETS_FILE)) {
            REALISTIC_TARGETS = new HashedListTargetCollection<>(reader.stream().collect(Collectors.toList()));
        } catch (final IOException ex) {
            throw new GATKException("could not read the realistic-targets file", ex);
        }
    }
    public static final XHMMModel[] TEST_MODELS = new XHMMModel[] {
        new XHMMModel(1e-4, 70_000, -3, 3)
    };

    public static File writeChainInTempFile(final XHMMData chain) {
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

    private static XHMMData simulateChain(final XHMMModel model, final TargetCollection<Target> targetCollection, final int sampleCount, final Random rdn) {
        final List<Target> targets = targetCollection.targets();
        final List<CopyNumberTriState> truth = new ArrayList<>(targets.size());
        final List<List<Double>> data = new ArrayList<>(sampleCount);
        for (int i = 0; i < sampleCount; i++) {
            final List<Double> sampleData = randomSampleData(model, rdn, targets, truth);
            data.add(sampleData);
        }
        return new XHMMData(model, targets, truth, data);
    }

    private static List<Double> randomSampleData(XHMMModel model, Random rdn, List<Target> targets, List<CopyNumberTriState> truth) {
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

    protected void loadModelArguments(XHMMData chain, List<String> arguments) {
        arguments.add("-" + XHMMArgumentCollection.MEAN_DELETION_COVERAGE_SHIFT_SHORT_NAME);
        arguments.add(String.valueOf(chain.model.getDeletionMean()));
        arguments.add("-" + XHMMArgumentCollection.MEAN_DUPLICATION_COVERAGE_SHIFT_SHORT_NAME);
        arguments.add(String.valueOf(chain.model.getDuplicationMean()));
        arguments.add("-" + XHMMArgumentCollection.EVENT_START_PROBABILITY_SHORT_NAME);
        arguments.add(String.valueOf(chain.model.getEventStartProbability()));
        arguments.add("-" + XHMMArgumentCollection.MEAN_EVENT_SIZE_SHORT_NAME);
        arguments.add(String.valueOf(chain.model.getMeanEventSize()));
        arguments.add("-" + XHMMSegmentCaller.ZSCORE_DIMENSION_SHORT_NAME);
    }

    public static class XHMMData {

        public final List<Target> targets;
        public final List<CopyNumberTriState> truth;
        public final List<List<Double>> data;
        public final XHMMModel model;

        protected XHMMData(final XHMMModel model, final List<Target> targets, final List<CopyNumberTriState> truth, final List<List<Double>> data) {
            this.targets = targets;
            this.truth = truth;
            this.data = data;
            this.model = model;
        }
    }
}
