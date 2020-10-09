package org.broadinstitute.hellbender.tools.dragstr;

import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.dragstr.DragstrParamUtils;
import org.broadinstitute.hellbender.utils.dragstr.DragstrParams;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;

/**
 * Tests the {@link DragstrParametersEstimator} code.
 */
public class DragstrParametersEstimatorTest {

    private static final Path DRAGEN_USED_CASES = new File("src/test/resources/large/org/broadinstitute/hellbender/tools/dragstr/dragen-param-estimator-input.tab.gz").toPath();
    private static final Path DRAGEN_ESPECTED_RESULT = new File("src/test/resources/large/org/broadinstitute/hellbender/tools/dragstr/dragen-param-estimator-expected-output.txt").toPath();
    private static final DragstrHyperParameters HYPER_PARAMETERS = new DragstrHyperParameters();

    @Test
    public void testMatchWithDragenValues()  {
        final StratifiedDragstrLocusCases cases = composeInputCases();
        final DragstrParametersEstimator subject = new DragstrParametersEstimator(HYPER_PARAMETERS);
        final DragstrParams actual = Utils.runInParallel(0, () -> subject.estimate(cases));
        final DragstrParams expected = DragstrParamUtils.parse(new GATKPath(DRAGEN_ESPECTED_RESULT.toString()));
        Assert.assertEquals(actual.maximumPeriod(), expected.maximumPeriod());
        Assert.assertEquals(actual.maximumRepeats(), expected.maximumRepeats());
        for (int i = 1; i <= actual.maximumPeriod(); i++) {
            for (int j = 1; j <= actual.maximumRepeats(); j++) {
                Assert.assertEquals(roundOf(actual.gop(i,j)), roundOf(expected.gop(i,j)));
                Assert.assertEquals(roundOf(actual.gcp(i,j)), roundOf(expected.gcp(i,j)));
                Assert.assertEquals(roundOf(actual.api(i,j)), roundOf(expected.api(i,j)));
            }
        }
    }

    private double roundOf(final double dbl) {
        return Math.round(100 * dbl) / 100.0;
    }

    private StratifiedDragstrLocusCases composeInputCases() {
        try (final DRAGENEstimationDataReader reader = new DRAGENEstimationDataReader(DRAGEN_USED_CASES)) {
            return reader.stream()
                    .map(DRAGENEstimationDataRecord::toCase)
                    .reduce(StratifiedDragstrLocusCases.make(HYPER_PARAMETERS.maxPeriod, HYPER_PARAMETERS.maxRepeatLength),
                            (accu, caze) -> { accu.add(caze); return accu; },
                            StratifiedDragstrLocusCases::merge);
        } catch (final IOException ex) {
            Assert.fail("can't read input data test file: " + DRAGEN_USED_CASES, ex);
            throw new RuntimeException("this should not be reached");
        }
    }

    static class DRAGENEstimationDataReader extends TableReader<DRAGENEstimationDataRecord> {
        DRAGENEstimationDataReader(final Path path) throws IOException {
            super(path);
        }

        @Override
        protected DRAGENEstimationDataRecord createRecord(DataLine dataLine) {
            return new DRAGENEstimationDataRecord(
                    dataLine.getInt(), // 0
                    dataLine.getLong(), // 1
                    dataLine.getInt(), // 2
                    dataLine.getInt(),  // 3
                    dataLine.seek(6).getInt(), // 6 because we skip operiod (4) and orepeat (5) columns.
                    dataLine.getInt()); // 7
        }
    }

    static class DRAGENEstimationDataRecord {
        private final int chrIdx;
        private final long pos;
        private final int period;
        private final int repeats;
        private final int n;
        private final int k;

        private DRAGENEstimationDataRecord(int chrIdx, long pos, int period, int repeats, int n, int k) {
            this.chrIdx = chrIdx;
            this.pos = pos;
            this.period = period;
            this.repeats = repeats;
            this.n = n;
            this.k = k;
        }

        private DragstrLocusCase toCase() {
            // 0 mask is irrelevant as these don't affect the test.
            final DragstrLocus locus = DragstrLocus.make(chrIdx, pos, (byte) period, (short) (period * repeats), 0);
            // 60,0 minMQ,nSup are not relevant as these do not affect the test.
            return DragstrLocusCase.create(locus, n, k, 60, 0);
        }
    }
}
