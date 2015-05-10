package org.broadinstitute.hellbender.tools.picard.interval;

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.*;

import static htsjdk.samtools.util.IntervalList.fromFile;
import static java.util.Arrays.asList;
import static org.broadinstitute.hellbender.tools.picard.interval.IntervalListScatterer.Mode;
import static org.broadinstitute.hellbender.tools.picard.interval.IntervalListScatterer.Mode.BALANCING_WITHOUT_INTERVAL_SUBDIVISION;
import static org.broadinstitute.hellbender.tools.picard.interval.IntervalListScatterer.Mode.INTERVAL_SUBDIVISION;
import static org.testng.Assert.assertEquals;

/**
 * Very basic test for scatter functionality in IntervalListTools
 */
public final class IntervalListScattererTest extends CommandLineProgramTest {
    private static final IntervalList LIST_TO_SCATTER;

    static {
        LIST_TO_SCATTER = fromFile(new File(getTestDataDir(), "picard/interval/scatterable.interval_list"));
        assertEquals(LIST_TO_SCATTER.getUniqueBaseCount(), 200, "Wrong unique base count");
    }

    private static class Testcase {
        final IntervalList source;
        final List<IntervalList> expectedScatter;
        final int scatterWidth;
        final Mode mode;

        @Override
        public String toString() {
            return "Testcase{" +
                    "scatterWidth=" + scatterWidth +
                    ", mode=" + mode +
                    '}';
        }

        private Testcase(final IntervalList source, final int scatterWidth, final Mode mode, final List<IntervalList> expectedScatter) {
            this.source = source;
            this.expectedScatter = expectedScatter;
            this.scatterWidth = scatterWidth;
            this.mode = mode;
        }
    }

    private static final List<Testcase> testcases = new ArrayList<>();

    static {
        testcases.add(new Testcase(
                LIST_TO_SCATTER, 2, INTERVAL_SUBDIVISION,
                asList(
                        composeIntervalList(LIST_TO_SCATTER, "1",
                                30000, 30098,
                                30100, 30100
                        ),
                        composeIntervalList(LIST_TO_SCATTER, "1",
                                30101, 30150,
                                30200, 30249
                        )
                )
        ));

        testcases.add(new Testcase(
                LIST_TO_SCATTER, 4, INTERVAL_SUBDIVISION,
                asList(
                        composeIntervalList(LIST_TO_SCATTER, "1",
                                30000, 30049
                        ),
                        composeIntervalList(LIST_TO_SCATTER, "1",
                                30050, 30098,
                                30100, 30100
                        ),
                        composeIntervalList(LIST_TO_SCATTER, "1",
                                30101, 30150
                        ),
                        composeIntervalList(LIST_TO_SCATTER, "1",
                                30200, 30249
                        )
                )
        ));

        testcases.add(new Testcase(
                LIST_TO_SCATTER, 5, INTERVAL_SUBDIVISION,
                asList(
                        composeIntervalList(LIST_TO_SCATTER, "1",
                                30000, 30039
                        ),
                        composeIntervalList(LIST_TO_SCATTER, "1",
                                30040, 30079
                        ),
                        composeIntervalList(LIST_TO_SCATTER, "1",
                                30080, 30098,
                                30100, 30120
                        ),
                        composeIntervalList(LIST_TO_SCATTER, "1",
                                30121, 30150,
                                30200, 30209
                        ),
                        composeIntervalList(LIST_TO_SCATTER, "1",
                                30210, 30249
                        )
                )
        ));

        testcases.add(new Testcase(
                LIST_TO_SCATTER, 6, BALANCING_WITHOUT_INTERVAL_SUBDIVISION,
                asList(
                        composeIntervalList(LIST_TO_SCATTER, "1",
                                30000, 30098
                        ),
                        composeIntervalList(LIST_TO_SCATTER, "1",
                                30100, 30150
                        ),
                        composeIntervalList(LIST_TO_SCATTER, "1",
                                30200, 30249
                        )
                )
        ));

        testcases.add(new Testcase(
                LIST_TO_SCATTER, 2, BALANCING_WITHOUT_INTERVAL_SUBDIVISION,
                asList(
                        composeIntervalList(LIST_TO_SCATTER, "1",
                                30000, 30098
                        ),
                        composeIntervalList(LIST_TO_SCATTER, "1",
                                30100, 30150,
                                30200, 30249
                        )
                )
        ));

        testcases.add(new Testcase(
                LIST_TO_SCATTER, 1, BALANCING_WITHOUT_INTERVAL_SUBDIVISION,
                asList(
                        composeIntervalList(LIST_TO_SCATTER, "1",
                                30000, 30098,
                                30100, 30150,
                                30200, 30249
                        )
                )
        ));
    }

    @DataProvider
    public Object[][] testScatterTestcases() {
        final Object[][] objects = new Object[testcases.size()][];
        for (int i = 0; i < objects.length; i++) {
            objects[i] = new Object[]{testcases.get(i)};
        }
        return objects;
    }

    @Test(dataProvider = "testScatterTestcases")
    public void testScatter(final Testcase tc) {
        final IntervalListScatterer scatterer = new IntervalListScatterer(tc.mode);
        final List<IntervalList> scatter = scatterer.scatter(tc.source, tc.scatterWidth);
        assertEquals(scatter, tc.expectedScatter);
    }

    private static IntervalList composeIntervalList(final IntervalList source, final String chromosome, final int... segmentsByPair) {
        final IntervalList intervals = new IntervalList(source.getHeader());
        for (int i = 0; i < segmentsByPair.length; i += 2) {
            final Interval parentInterval = lookupIntervalContainingLocus(source, chromosome, segmentsByPair[i]);
            intervals.add(new Interval("1", segmentsByPair[i], segmentsByPair[i + 1], parentInterval.isNegativeStrand(), parentInterval.getName()));
        }
        return intervals;
    }

    private static Interval lookupIntervalContainingLocus(final IntervalList source, final String chromosome, final int startIndex) {
        for (final Interval interval : source) {
            if (interval.getContig().equals(chromosome) && startIndex >= interval.getStart() && startIndex <= interval.getEnd()) {
                return interval;
            }
        }
        throw new NoSuchElementException();
    }
}
