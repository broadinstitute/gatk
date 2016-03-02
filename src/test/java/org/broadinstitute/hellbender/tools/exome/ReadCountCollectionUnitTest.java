package org.broadinstitute.hellbender.tools.exome;

import htsjdk.samtools.util.TestUtil;
import org.apache.commons.collections4.list.SetUniqueList;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.stat.descriptive.moment.Variance;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Units tests for {@link ReadCountCollection}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class ReadCountCollectionUnitTest extends BaseTest {

    private static final int[] CORRECT_COLUMN_COUNTS = {1, 2, 3, 6, 12};

    private static final int[] CORRECT_TARGET_COUNTS = {1, 2, 3, 6, 12, 101, 1001, 10001};

    @Test(dataProvider="correctInstantiationData")
    public void testCorrectInstantiation(final ReadCountCollectionInfo info) {
        final ReadCountCollection readCountCollection = info.newInstance();
        Assert.assertNotNull(readCountCollection);
    }

    @Test(dataProvider ="correctInstantiationData")
    public void testNormalizeBySampleAverageNoWeighting(final ReadCountCollectionInfo info) {
        final boolean weightByTargetSize = false;
        final ReadCountCollection counts = info.newInstance();
        final ReadCountCollection normalizedCounts = counts.normalizeByColumnAverages(weightByTargetSize);
        assertNormalizedBySampleAverage(counts, normalizedCounts, weightByTargetSize);
    }

    @Test(dataProvider ="correctInstantiationData")
    public void testNormalizeBySampleAverageWeighting(final ReadCountCollectionInfo info) {
        if (info.intervals == null) { return; }
        final boolean weightByTargetSize = true;
        final ReadCountCollection counts = info.newInstance();
        final ReadCountCollection normalizedCounts = counts.normalizeByColumnAverages(weightByTargetSize);
        assertNormalizedBySampleAverage(counts, normalizedCounts, weightByTargetSize);
    }

    @Test(dataProvider ="correctInstantiationData")
    public void testConvertToZScores(final ReadCountCollectionInfo info) {
        final ReadCountCollection counts = info.newInstance();
        final ReadCountCollection zScoreCounts = counts.zScoreCounts();
        assertZScores(counts, zScoreCounts);
    }

    @Test(dataProvider ="correctInstantiationData", expectedExceptions = IllegalArgumentException.class)
    public void testNormalizeBySampleAverageWithweightingNoTargetIntervals(final ReadCountCollectionInfo info) {
        if (info.intervals != null) {
            throw new IllegalArgumentException("this instance given by the DataProvider has intervals and is thus irrelevant");
        }
        final boolean weightByTargetSize = true;
        final ReadCountCollection counts = info.newInstance();
        final ReadCountCollection normalizedCounts = counts.normalizeByColumnAverages(weightByTargetSize);
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

    private void assertNormalizedBySampleAverage(final ReadCountCollection counts,
                                                 final ReadCountCollection normalizedCounts, final boolean weightByTargetSize) {
        Assert.assertEquals(counts.targets(), normalizedCounts.targets());
        Assert.assertEquals(counts.columnNames(), normalizedCounts.columnNames());

        final RealVector weights = new ArrayRealVector(counts.targets().stream()
                .mapToDouble(t -> weightByTargetSize ? t.length() : 1).toArray());

        final RealMatrix countsMatrix = counts.counts();
        final RealMatrix normalizedCountsMatrix = normalizedCounts.counts();

        Assert.assertEquals(countsMatrix.getRowDimension(), normalizedCountsMatrix.getRowDimension());
        Assert.assertEquals(countsMatrix.getColumnDimension(), normalizedCountsMatrix.getColumnDimension());

        for (int col = 0; col < countsMatrix.getColumnDimension(); col++) {
            final RealVector countsColumn = countsMatrix.getColumnVector(col);
            final RealVector normalizedCountsColumn = normalizedCountsMatrix.getColumnVector(col);

            //check that the two columns are parallel
            Assert.assertEquals(Math.abs(countsColumn.cosine(normalizedCountsColumn)), 1, 1e-8);

            final double inputWeightedAverage = countsColumn.dotProduct(weights) / weights.getL1Norm();
            Assert.assertEquals(normalizedCountsColumn.mapMultiply(inputWeightedAverage).getL1Distance(countsColumn), 0, 1e-8);
        }
    }

    private void assertZScores(final ReadCountCollection counts, final ReadCountCollection zScoreCounts) {
        Assert.assertEquals(counts.targets(), zScoreCounts.targets());
        Assert.assertEquals(counts.columnNames(), zScoreCounts.columnNames());

        final RealMatrix countsMatrix = counts.counts();
        final RealMatrix zScoreCountsMatrix = zScoreCounts.counts();

        for (int row = 0; row < countsMatrix.getRowDimension(); row++) {
            final double mean = GATKProtectedMathUtils.mean(countsMatrix.getRow(row));
            final double std = Math.sqrt(new Variance().evaluate(countsMatrix.getRow(row)));
            for (int col = 0; col < countsMatrix.getColumnDimension(); col++) {
                if (countsMatrix.getColumnDimension() > 1) {
                    Assert.assertEquals(zScoreCountsMatrix.getEntry(row, col), (countsMatrix.getEntry(row, col) - mean) / std, 1e-8);
                } else {
                    Assert.assertTrue(Double.isNaN(zScoreCountsMatrix.getEntry(row, col)));
                }
            }
        }

    }
}
