package org.broadinstitute.hellbender.tools.exome;

import org.apache.commons.collections4.list.SetUniqueList;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.stream.Collectors;

/**
 * Units tests for {@link ReadCountCollection}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class ReadCountCollectionUnitTest {

    private static final int[] CORRECT_COLUMN_COUNTS = {1, 2, 3, 6, 12};

    private static final int[] CORRECT_TARGET_COUNTS = {1, 2, 3, 6, 12, 101, 1001, 10001};

    @Test(dataProvider="correctInstantiationData")
    public void testCorrectInstantiation(final ReadCountCollectionInfo info) {
        final ReadCountCollection readCountCollection = info.newInstance();
        Assert.assertNotNull(readCountCollection);
    }

    @Test(dataProvider="correctInstantiationData",dependsOnMethods = "testCorrectInstantiation")
    public void testCounts(final ReadCountCollectionInfo info) {
        final ReadCountCollection readCountCollection = info.newInstance();
        final RealMatrix countMatrix = readCountCollection.counts();
        Assert.assertNotNull(countMatrix);
        Assert.assertEquals(countMatrix.getRowDimension(), info.targetCount);
        Assert.assertEquals(countMatrix.getColumnDimension(), info.columnCount);
        for (int i = 0; i < countMatrix.getRowDimension(); i++) {
            for (int j = 0; j < countMatrix.getColumnDimension(); j++) {
                Assert.assertEquals(countMatrix.getEntry(i,j),info.counts[i][j]);
            }
        }
    }

    @Test(dataProvider="correctInstantiationData", dependsOnMethods = "testCorrectInstantiation")
    public void testTargetCount(final ReadCountCollectionInfo info) {
        final ReadCountCollection readCountCollection = info.newInstance();
        Assert.assertEquals(readCountCollection.targets().size(), info.targetCount);
    }

    @Test(dataProvider="correctInstantiationData",dependsOnMethods = "testCorrectInstantiation")
    public void testTargetNames(final ReadCountCollectionInfo info) {
        final ReadCountCollection readCountCollection = info.newInstance();
        final List<String> targetNames = readCountCollection.targets().stream().map(Target::getName).collect(Collectors.toList());
        Assert.assertNotNull(targetNames);
        Assert.assertEquals(new ArrayList<>(targetNames), new ArrayList<>(info.targetNames));
    }

    @Test(dataProvider="correctInstantiationData",dependsOnMethods = "testCorrectInstantiation")
    public void testIntervals(final ReadCountCollectionInfo info) {
        final ReadCountCollection readCountCollection = info.newInstance();
        if (info.intervals != null) {
            final List<SimpleInterval> intervals = readCountCollection.targets().stream().map(Target::getInterval).collect(Collectors.toList());
            Assert.assertEquals(new ArrayList<>(intervals),new ArrayList<>(info.intervals));
        } else {
            Assert.assertFalse(readCountCollection.targets().stream().anyMatch(t -> t.getInterval() != null));
        }
    }

    @Test(dataProvider="wrongInstantiationData",
            expectedExceptions = IllegalArgumentException.class)
    public void testWrongInstantiation(final ReadCountCollectionInfo info, final String caseName) {
        try {
            info.newInstance();
        } catch (final RuntimeException ex) {
            org.testng.log4testng.Logger.getLogger(getClass()).info("testWrongInstantiation exception message: " + ex.getMessage());
            throw ex;
        }
        Assert.fail("Exception not thrown: case " + caseName);
    }

    private static class ReadCountCollectionInfo {
        private final int columnCount;
        private int targetCount;
        private List<String> columnNames;
        private List<String> targetNames;
        private List<SimpleInterval> intervals;
        private double[][] counts;

        private ReadCountCollectionInfo(final int columnCount,
                                        final int targetCount,
                                        final List<String> columnNames,
                                        final List<String> targetNames,
                                        final List<SimpleInterval> intervals,
                                        final double[][] counts) {
            this.columnCount = columnCount;
            this.targetCount = targetCount;
            this.columnNames = columnNames;
            this.targetNames = targetNames;
            this.intervals = intervals;
            this.counts = counts;
        }

        private ReadCountCollection newInstance() {
            final List<Target> targets = new ArrayList<>();
            for (int i = 0; i < targetNames.size(); i++) {
                final String targetName = targetNames.get(i);
                targets.add(targetName == null ? null : new Target(targetName,intervals == null ? null : intervals.get(i)));
            }

            return new ReadCountCollection(
                    SetUniqueList.setUniqueList(targets),
                    SetUniqueList.setUniqueList(columnNames),
                    new Array2DRowRealMatrix(counts));
        }

        public ReadCountCollectionInfo makeATargetNameNull(final int index) {
            targetNames.set(index,null);
            return this;
        }

        public ReadCountCollectionInfo makeAColumnNameNull(final int index) {
            columnNames.set(index,null);
            return this;
        }

        public ReadCountCollectionInfo changeTargetCount(final int delta) {
           this.targetCount -= delta;
            if (delta > 0) {
                for (int i = 0; i < delta; i++) {
                    this.targetNames.add(targetNames.get(i) + "_extra");
                    if (intervals != null) {
                        final SimpleInterval interval = intervals.get(i);
                        this.intervals.add(intervals.get(i) == null ? null : new SimpleInterval("other_chr",interval.getStart(),interval.getEnd()));
                    }
                }
            } else if (delta < 0) {
                if (intervals != null) {
                    intervals = intervals.subList(0,intervals.size() + delta);
                }
                targetNames = targetNames.subList(0,targetNames.size() + delta);
            }
            this.targetCount += delta;
            return this;
        }

        public ReadCountCollectionInfo changeColumnCount(final int delta) {
            if (delta > 0) {
                for (int i = 0; i < delta; i++)
                columnNames.add(columnNames.get(i) + "_extra");
            } else if (delta < 0) {
                columnNames = columnNames.subList(0,columnNames.size() + delta);
            }
            return this;
        }
    }

    private ReadCountCollectionInfo newInstanceInfo(final int columnCount, final int targetCount, final boolean withIntervals) {
        final Random rdn = new Random((11 + columnCount) * (13 * targetCount) * (withIntervals ? 37 : 9) );
        final List<String> columnNames = new ArrayList<>(columnCount);
        for (int i = 0; i < columnCount; i++) {
            columnNames.add("col_" + Math.abs(rdn.nextInt()) + "_" + i );
        }
        final List<String> targetNames;
        targetNames = new ArrayList<>(targetCount);
        for (int j = 0; j < targetCount; j++) {
            targetNames.add("tgt_" + Math.abs(rdn.nextInt()) + "_" + j);
        }
        final double[][] counts = new double[targetCount][columnCount];
        for (int i = 0; i < counts.length; i++) {
            for (int j = 0; j < counts[i].length; j++) {
                counts[i][j] = rdn.nextGaussian();
            }
        }
        final List<SimpleInterval> intervals;
        if (withIntervals) {
            intervals = new ArrayList<>(targetCount);
            int start = 1;
            for (int j = 0; j < targetCount; j++) {
                final int newStart = start + rdn.nextInt(100);
                final int newEnd = newStart + 1 + rdn.nextInt(200);
                intervals.add(new SimpleInterval("seq1",newStart,newEnd));
                start = newEnd + 1;
            }
        } else {
            intervals = null;
        }
        return new ReadCountCollectionInfo(columnCount,targetCount,columnNames,targetNames,intervals,counts);
    }

    @DataProvider(name="wrongInstantiationData")
    public Object[][] wrongInstantiationData() {
        final List<Object[]> result = new ArrayList<>();
        result.addAll(Arrays.asList(new Object[]{newInstanceInfo(0, 2, true), "0"},
                new Object[]{newInstanceInfo(10, 100, true).makeATargetNameNull(78), "12"},
                new Object[]{newInstanceInfo(10, 100, true).changeTargetCount(-1), "12"},
                new Object[]{newInstanceInfo(10, 100, true).changeColumnCount(1),"12"},
                new Object[]{newInstanceInfo(10, 100, true).makeAColumnNameNull(5), "13"}));
        return result.toArray(new Object[result.size()][]);
    }

    @DataProvider(name="correctInstantiationData")
    public Object[][] correctInstantiationData() {
        final Object[][] result = new Object[CORRECT_COLUMN_COUNTS.length * CORRECT_TARGET_COUNTS.length * 2][];
        for (int i = 0; i < CORRECT_COLUMN_COUNTS.length; i++) {
            for (int j = 0; j < CORRECT_TARGET_COUNTS.length; j++) {
                final int baseIndex = 2 * (i * (CORRECT_TARGET_COUNTS.length) + j);
                result[baseIndex] = new Object[] { newInstanceInfo(CORRECT_COLUMN_COUNTS[i], CORRECT_TARGET_COUNTS[j], true) };
                result[baseIndex + 1] = new Object[] { newInstanceInfo (CORRECT_COLUMN_COUNTS[i], CORRECT_TARGET_COUNTS[j], false)};
            }
        }
        return result;
    }



}
