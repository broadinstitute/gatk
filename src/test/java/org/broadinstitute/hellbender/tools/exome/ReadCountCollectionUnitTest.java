package org.broadinstitute.hellbender.tools.exome;

import htsjdk.samtools.util.TestUtil;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import org.testng.internal.junit.ArrayAsserts;

import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Units tests for {@link ReadCountCollection}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class ReadCountCollectionUnitTest extends GATKBaseTest {

    private static final int[] CORRECT_COLUMN_COUNTS = {1, 10};

    private static final int[] CORRECT_TARGET_COUNTS = {1, 10, 100};

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
            logger.info("testWrongInstantiation exception message: " + ex.getMessage());
            throw ex;
        }
        Assert.fail("Exception not thrown: case " + caseName);
    }

    @Test(dataProvider="targetArrangeData")
    public void testArrangeTargets(final ReadCountCollectionInfo info, final List<String> newOrder) {
        final List<Target> targetsNewOrder = newOrder.stream()
                .map(Target::new).collect(Collectors.toList());
        final ReadCountCollection subject = info.newInstance();
        final ReadCountCollection result = subject.arrangeTargets(targetsNewOrder);
        final List<Target> afterTargets = result.targets();
        Assert.assertEquals(afterTargets.size(),  targetsNewOrder.size());
        Assert.assertFalse(afterTargets.stream().anyMatch(t -> !targetsNewOrder.contains(t)));
        final RealMatrix beforeCounts = subject.counts();
        final RealMatrix afterCounts = result.counts();
        Assert.assertEquals(beforeCounts.getColumnDimension(), afterCounts.getColumnDimension());
        Assert.assertEquals(afterCounts.getRowDimension(), targetsNewOrder.size());
        final int[] beforeIndexes = new int[targetsNewOrder.size()];
        final int[] afterIndexes = new int[targetsNewOrder.size()];
        int nextIdx = 0;
        for (final Target target : targetsNewOrder) {
            final int beforeIndex = subject.targets().indexOf(target);
            final int afterIndex = result.targets().indexOf(target);
            beforeIndexes[nextIdx] = beforeIndex;
            afterIndexes[nextIdx++] = afterIndex;
        }
        // check that the counts are exactly the same.
        for (int i = 0; i < beforeIndexes.length; i++) {
            final double[] before = beforeCounts.getRow(beforeIndexes[i]);
            final double[] after = afterCounts.getRow(afterIndexes[i]);
            Assert.assertEquals(before, after);
        }
    }


    @Test(dataProvider="targetSubsetData")
    public void testSubsetTargets(final ReadCountCollectionInfo info, final Set<String> targetNamesToKeep) {
        final Set<Target> targetsToKeep = targetNamesToKeep.stream()
                .map(Target::new).collect(Collectors.toSet());
        final ReadCountCollection subject = info.newInstance();
        final ReadCountCollection result = subject.subsetTargets(targetsToKeep);
        final List<Target> afterTargets = result.targets();
        Assert.assertEquals(afterTargets.size(),  targetsToKeep.size());
        Assert.assertFalse(afterTargets.stream().anyMatch(t -> !targetsToKeep.contains(t)));
        final RealMatrix beforeCounts = subject.counts();
        final RealMatrix afterCounts = result.counts();
        Assert.assertEquals(beforeCounts.getColumnDimension(), afterCounts.getColumnDimension());
        Assert.assertEquals(afterCounts.getRowDimension(), targetNamesToKeep.size());
        final int[] beforeIndexes = new int[targetsToKeep.size()];
        final int[] afterIndexes = new int[targetsToKeep.size()];
        int nextIdx = 0;
        for (final Target target : targetsToKeep) {
            final int beforeIndex = subject.targets().indexOf(target);
            final int afterIndex = result.targets().indexOf(target);
            beforeIndexes[nextIdx] = beforeIndex;
            afterIndexes[nextIdx++] = afterIndex;
        }
        // check that the order of targets in the output is kept to the original order.
        for (int i = 0; i < beforeIndexes.length; i++) {
            for (int j = 0; j < beforeIndexes.length; j++) {
                Assert.assertEquals(Integer.compare(beforeIndexes[i], beforeIndexes[j]),
                        Integer.compare(afterIndexes[i], afterIndexes[j]));
            }
        }
        // check that the counts are exactly the same.
        for (int i = 0; i < beforeIndexes.length; i++) {
            final double[] before = beforeCounts.getRow(beforeIndexes[i]);
            final double[] after = afterCounts.getRow(afterIndexes[i]);
            Assert.assertEquals(before, after);
        }
    }

    @Test(dataProvider="wrongTargetSubsetData", expectedExceptions = IllegalArgumentException.class)
    public void testTargetSubsetColumns(final ReadCountCollectionInfo info, final Set<String> targetNamesToKeep) {
        if (targetNamesToKeep == null) {
            info.newInstance().subsetTargets(null);
        } else {
            final Set<Target> targets = targetNamesToKeep.stream().map(Target::new).collect(Collectors.toSet());
            info.newInstance().subsetTargets(targets);
        }
    }

    @Test(dataProvider="wrongTargetArrangeData", expectedExceptions = IllegalArgumentException.class)
    public void testTargetSubsetColumns(final ReadCountCollectionInfo info, final List<String> newOrder) {
        if (newOrder == null) {
            info.newInstance().subsetTargets(null);
        } else {
            final List<Target> targets = newOrder.stream().map(Target::new).collect(Collectors.toList());
            info.newInstance().arrangeTargets(targets);
        }
    }

    @Test(dataProvider="wrongColumnSubsetData", expectedExceptions = IllegalArgumentException.class)
    public void testWrongSubsetColumns(final ReadCountCollectionInfo info, final Set<String> columnsToKeep) {
        info.newInstance().subsetColumns(columnsToKeep);
    }

    @Test(dataProvider="columnSubsetData")
    public void testSubsetColumns(final ReadCountCollectionInfo info, final Set<String> columnsToKeep) {
        final ReadCountCollection subject = info.newInstance();
        final ReadCountCollection result = subject.subsetColumns(columnsToKeep);
        final List<String> afterColumns = result.columnNames();
        Assert.assertEquals(afterColumns.size(),  columnsToKeep.size());
        Assert.assertFalse(afterColumns.stream().anyMatch(t -> !columnsToKeep.contains(t)));
        final RealMatrix beforeCounts = subject.counts();
        final RealMatrix afterCounts = result.counts();
        Assert.assertEquals(afterCounts.getColumnDimension(), columnsToKeep.size());
        Assert.assertEquals(afterCounts.getRowDimension(), beforeCounts.getRowDimension());
        final int[] beforeIndexes = new int[columnsToKeep.size()];
        final int[] afterIndexes = new int[columnsToKeep.size()];
        int nextIdx = 0;
        for (final String column : columnsToKeep) {
            final int beforeIndex = subject.columnNames().indexOf(column);
            final int afterIndex = result.columnNames().indexOf(column);
            beforeIndexes[nextIdx] = beforeIndex;
            afterIndexes[nextIdx++] = afterIndex;
        }
        // check that the order of targets in the output is kept to the original order.
        for (int i = 0; i < beforeIndexes.length; i++) {
            for (int j = 0; j < beforeIndexes.length; j++) {
                Assert.assertEquals(Integer.compare(beforeIndexes[i], beforeIndexes[j]),
                        Integer.compare(afterIndexes[i], afterIndexes[j]));
            }
        }
        // check that the counts are exactly the same.
        for (int i = 0; i < beforeIndexes.length; i++) {
            final double[] before = beforeCounts.getColumn(beforeIndexes[i]);
            final double[] after = afterCounts.getColumn(afterIndexes[i]);
            Assert.assertEquals(before, after);
        }
    }

    @Test(dataProvider = "getColumnOnSpecifiedTargetsTestData")
    public void testGetColumnOnSpecifiedTargets(final ReadCountCollection rcc,
                                                final int columnIndex,
                                                final List<Target> targets, final double[] expected) {
        final double[] result = rcc.getColumnOnSpecifiedTargets(columnIndex, targets, true);
        ArrayAsserts.assertArrayEquals(expected, result, 1e-12);
    }

    @Test(dataProvider = "getRowByIndexTestData")
    public void testGetRowByIndex(final ReadCountCollection rcc, final int rowIndex, final double[] expected) {
        final double[] result = rcc.getRow(rowIndex);
        ArrayAsserts.assertArrayEquals(expected, result, 1e-12);
    }

    @Test(dataProvider = "getRowByTargetTestData")
    public void testGetRowByTarget(final ReadCountCollection rcc, final Target target, final double[] expected) {
        final double[] result = rcc.getRow(target);
        ArrayAsserts.assertArrayEquals(expected, result, 1e-12);
    }

    @Test(dataProvider = "getRowIndexOutOfRangeTestData", expectedExceptions = IllegalArgumentException.class)
    public void testGetRowIndexOutOfRange(final ReadCountCollection rcc, final int rowIndex) {
        rcc.getRow(rowIndex);
    }

    @Test(dataProvider = "getRowNonexistentTargetTestData", expectedExceptions = IllegalArgumentException.class)
    public void testGetRowNonexistentTarget(final ReadCountCollection rcc, final Target nonexistentTarget) {
        rcc.getRow(nonexistentTarget);
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

            return new ReadCountCollection(targets, columnNames, new Array2DRowRealMatrix(counts));
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
            targetCount -= delta;
            if (delta > 0) {
                for (int i = 0; i < delta; i++) {
                    targetNames.add(targetNames.get(i) + "_extra");
                    if (intervals != null) {
                        final SimpleInterval interval = intervals.get(i);
                        intervals.add(intervals.get(i) == null ? null : new SimpleInterval("other_chr",interval.getStart(),interval.getEnd()));
                    }
                }
            } else if (delta < 0) {
                if (intervals != null) {
                    intervals = intervals.subList(0,intervals.size() + delta);
                }
                targetNames = targetNames.subList(0,targetNames.size() + delta);
            }
            targetCount += delta;
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

    @DataProvider(name = "getColumnOnSpecifiedTargetsTestData")
    public Object[][] getColumnOnSpecifiedTargetsTestData() {
        final ReadCountCollectionInfo info = newInstanceInfo(3, 10, true);
        final ReadCountCollection rcc = info.newInstance();
        final List<Target> allTargets = rcc.targets();
        final List<Target> newTargets_1 = Arrays.stream(new Target[] {
                allTargets.get(0), allTargets.get(2), allTargets.get(4), allTargets.get(6)})
                .collect(Collectors.toList());
        final List<Target> newTargets_2 = Arrays.stream(new Target[] {
                allTargets.get(3), allTargets.get(1), allTargets.get(5)})
                .collect(Collectors.toList());
        final List<Target> newTargets_3 = Arrays.stream(new Target[] {
                allTargets.get(2), allTargets.get(2), allTargets.get(8)})
                .collect(Collectors.toList());
        final double[] col1 = rcc.getColumn(1);
        final double[] col2 = rcc.getColumn(2);
        return new Object[][] {
                {rcc, 1, allTargets, col1},
                {rcc, 1, newTargets_1, new double[] {col1[0], col1[2], col1[4], col1[6]}},
                {rcc, 1, newTargets_2, new double[] {col1[3], col1[1], col1[5]}},
                {rcc, 1, newTargets_3, new double[] {col1[2], col1[2], col1[8]}},
                {rcc, 2, allTargets, col2},
                {rcc, 2, newTargets_1, new double[] {col2[0], col2[2], col2[4], col2[6]}},
                {rcc, 2, newTargets_2, new double[] {col2[3], col2[1], col2[5]}},
                {rcc, 2, newTargets_3, new double[] {col2[2], col2[2], col2[8]}}
        };
    }

    @DataProvider(name = "getRowByIndexTestData")
    public Object[][] getRowByIndexTestData() {
        final ReadCountCollectionInfo info = newInstanceInfo(3, 10, true);
        final ReadCountCollection rcc = info.newInstance();
        return new Object[][] {
                {rcc, 0, rcc.counts().getRow(0)},
                {rcc, 1, rcc.counts().getRow(1)},
                {rcc, 2, rcc.counts().getRow(2)}
        };
    }

    @DataProvider(name = "getRowByTargetTestData")
    public Object[][] getRowByTargetTestData() {
        final ReadCountCollectionInfo info = newInstanceInfo(3, 10, true);
        final ReadCountCollection rcc = info.newInstance();
        return new Object[][] {
                {rcc, rcc.targets().get(0), rcc.counts().getRow(0)},
                {rcc, rcc.targets().get(1), rcc.counts().getRow(1)},
                {rcc, rcc.targets().get(2), rcc.counts().getRow(2)}
        };
    }

    @DataProvider(name = "getRowIndexOutOfRangeTestData")
    public Object[][] getRowIndexOutOfRangeTestData() {
        final ReadCountCollectionInfo info = newInstanceInfo(3, 10, true);
        final ReadCountCollection rcc = info.newInstance();
        return new Object[][] {{rcc, -5}, {rcc, -1}, {rcc, 10}, {rcc, 15}};
    }

    @DataProvider(name = "getRowNonexistentTargetTestData")
    public Object[][] getRowNonexistentTargetTestData() {
        final ReadCountCollectionInfo info = newInstanceInfo(3, 10, true);
        final ReadCountCollection rcc = info.newInstance();
        return new Object[][] {{rcc, new Target("MY_FABULOUS_NONEXISTENT_TARGET")}};
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

    @DataProvider(name="targetArrangeData")
    public Object[][] targetArrangeData() {
        final Object[][] correctInstantiationData = correctInstantiationData();
        final Random random = new Random(13);
        final List<Object[]> result = new ArrayList<>(correctInstantiationData.length * 4);
        for (final Object[] params : correctInstantiationData) {
            final ReadCountCollectionInfo info = (ReadCountCollectionInfo) params[0];
            final List<String> targetNames = info.targetNames;
            // no actual sub-setting.
            if (info.targetNames.size() > 0) {
                final List<String> shuffledTargetNames = new ArrayList<>(targetNames);
                Collections.shuffle(shuffledTargetNames, random);
                        result.add(new Object[]{info, shuffledTargetNames });
                result.add(new Object[]{info, Collections.singletonList(targetNames.get(0))});
                result.add(new Object[]{info, Collections.singletonList(targetNames.get(0))});
                if (info.targetNames.size() >= 2) {
                    result.add(new Object[]{info, shuffledTargetNames.subList(0, targetNames.size() / 2)});
                }
            }
        }
        return result.toArray(new Object[result.size()][]);
    }

    @DataProvider(name="targetSubsetData")
    public Object[][] targetSubsetData() {
        final Object[][] correctInstantiationData = correctInstantiationData();
        final List<Object[]> result = new ArrayList<>(correctInstantiationData.length * 4);
        for (final Object[] params : correctInstantiationData) {
            final ReadCountCollectionInfo info = (ReadCountCollectionInfo) params[0];
            final List<String> targetNames = info.targetNames;
            // no actual sub-setting.
            result.add(new Object[] { info, new HashSet<>( targetNames ) });
            if (info.targetNames.size() > 0) {
                result.add(new Object[]{info, Collections.singleton(targetNames.get(0))});
                if (info.targetNames.size() > 1) {
                    result.add(new Object[] { info, Collections.singleton(targetNames.get(targetNames.size() - 1))});
                    if (info.targetNames.size() > 2) {
                        result.add(new Object[] { info, Collections.singleton(targetNames.get(targetNames.size() / 2))});
                        result.add(new Object[] { info,
                                new HashSet<>(Arrays.asList(targetNames.get(0), targetNames.get(targetNames.size() - 1))) });
                    }
                }
            }
        }
        return result.toArray(new Object[result.size()][]);
    }

    @DataProvider(name="columnSubsetData")
    public Object[][] columnSubsetData() {
        final Object[][] correctInstantiationData = correctInstantiationData();
        final List<Object[]> result = new ArrayList<>(correctInstantiationData.length * 4);
        for (final Object[] params : correctInstantiationData) {
            final ReadCountCollectionInfo info = (ReadCountCollectionInfo) params[0];
            final List<String> columnNames = info.columnNames;
            // no actual sub-setting.
            result.add(new Object[] { info, new HashSet<>( columnNames ) });
            if (info.columnNames.size() > 0) {
                result.add(new Object[]{info, Collections.singleton(columnNames.get(0))});
                if (info.columnNames.size() > 1) {
                    result.add(new Object[] { info, Collections.singleton(columnNames.get(columnNames.size() - 1))});
                    if (info.columnNames.size() > 2) {
                        result.add(new Object[] { info, Collections.singleton(columnNames.get(columnNames.size() / 2))});
                        result.add(new Object[] { info,
                                new HashSet<>(Arrays.asList(columnNames.get(0), columnNames.get(columnNames.size() - 1))) });
                    }
                }
            }
        }
        return result.toArray(new Object[result.size()][]);
    }

    @DataProvider(name="wrongColumnSubsetData")
    public Object[][] wrongColumnSubsetData() {
        final Object[][] correctInstantiationData = correctInstantiationData();
        final List<Object[]> result = new ArrayList<>(correctInstantiationData.length * 4);
        for (final Object[] params : correctInstantiationData) {
            final ReadCountCollectionInfo info = (ReadCountCollectionInfo) params[0];
            final List<String> columnNames = info.columnNames;
            // try to remove all columns:
            result.add(new Object[] { info, new HashSet<>( ) });
            result.add(new Object[] { info, Collections.singleton("NOT_A_COLUMN")});
            result.add(new Object[] { info, null });
            result.add(new Object[] { info, Collections.singleton(null)});
        }
        return result.toArray(new Object[result.size()][]);
    }

    @DataProvider(name="wrongTargetSubsetData")
    public Object[][] wrongTargetSubsetData() {
        final Object[][] correctInstantiationData = correctInstantiationData();
        final List<Object[]> result = new ArrayList<>(correctInstantiationData.length * 4);
        for (final Object[] params : correctInstantiationData) {
            final ReadCountCollectionInfo info = (ReadCountCollectionInfo) params[0];
            // try to remove all columns:
            result.add(new Object[] { info, new HashSet<>( ) });
            result.add(new Object[] { info, Collections.singleton("NOT_A_TARGET")});
            result.add(new Object[] { info, null });
            result.add(new Object[] { info, Collections.singleton(null)});
        }
        return result.toArray(new Object[result.size()][]);
    }

    @DataProvider(name="wrongTargetArrangeData")
    public Object[][] wrongTargetArrangeData() {
        final Object[][] correctInstantiationData = correctInstantiationData();
        final List<Object[]> result = new ArrayList<>(correctInstantiationData.length * 4);
        for (final Object[] params : correctInstantiationData) {
            final ReadCountCollectionInfo info = (ReadCountCollectionInfo) params[0];
            // try to remove all columns:
            result.add(new Object[] { info, new ArrayList<>( ) });
            result.add(new Object[] { info, Collections.singletonList("NOT_A_TARGET")});
            result.add(new Object[] { info, null });
            result.add(new Object[] { info, Collections.singletonList(null)});
        }
        return result.toArray(new Object[result.size()][]);
    }

    @Test(dataProvider = "correctInstantiationData")
    public void testSerializeDeserialize(final ReadCountCollectionInfo info) {
        final ReadCountCollection newInstance = info.newInstance();
        try {
            TestUtil.serializeAndDeserialize(newInstance);
        } catch (final IOException | ClassNotFoundException e) {
            Assert.fail("Exception thrown", e);
        }
    }

}
