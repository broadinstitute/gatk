package org.broadinstitute.hellbender.utils.recalibration;

import org.broadinstitute.hellbender.tools.dataflow.transforms.bqsr.RecalibrationTablesMerger;
import org.broadinstitute.hellbender.utils.collections.NestedIntegerArray;
import org.broadinstitute.hellbender.utils.recalibration.covariates.StandardCovariateList;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Tests RecalibrationTablesMerger. This is based on the RecalibrationTablesUnitTest, since the merging behavior should
 * be the same if everything's working as expected.
 */
public final class RecalibrationTablesMergerUnitTest extends BaseTest {
    private RecalibrationTables tables;
    private StandardCovariateList covariates;
    private int numReadGroups = 6;
    final byte qualByte = 1;
    final List<Integer> combineStates = Arrays.asList(0, 1, 2);

    @BeforeMethod
    private void makeTables() {
        final List<String> readGroups= IntStream.range(1, numReadGroups).mapToObj(i -> "readgroup"+i).collect(Collectors.toList());
        covariates = new StandardCovariateList(new RecalibrationArgumentCollection(), readGroups);
        tables = new RecalibrationTables(covariates, numReadGroups);
        fillTable(tables);
    }

    private void fillTable(final RecalibrationTables tables) {
        for ( int iterations = 0; iterations < 10; iterations++ ) {
            for ( final EventType et : EventType.values() ) {
                for ( final int rg : combineStates) {
                    final double error = rg % 2 == 0 ? 1 : 0;
                    RecalUtils.incrementDatumOrPutIfNecessary(tables.getReadGroupTable(), qualByte, error, rg, et.ordinal());
                    for ( final int qual : combineStates) {
                        RecalUtils.incrementDatumOrPutIfNecessary(tables.getQualityScoreTable(), qualByte, error, rg, qual, et.ordinal());
                        for ( final int cycle : combineStates)
                            RecalUtils.incrementDatumOrPutIfNecessary(tables.getTable(2), qualByte, error, rg, qual, cycle, et.ordinal());
                        for ( final int context : combineStates)
                            RecalUtils.incrementDatumOrPutIfNecessary(tables.getTable(3), qualByte, error, rg, qual, context, et.ordinal());
                    }
                }
            }
        }
    }

    @Test
    public void testCombine1() {
        final RecalibrationTables toMerge = new RecalibrationTables(covariates, numReadGroups);
        fillTable(toMerge);

        RecalibrationTablesMerger merger = new RecalibrationTablesMerger();
        RecalibrationTables acc = merger.createAccumulator();
        acc = merger.addInput(acc, toMerge);
        acc = merger.addInput(acc, tables);
        RecalibrationTables merged = merger.extractOutput(acc);

        for ( int i = 0; i < tables.numTables(); i++ ) {
            NestedIntegerArray<RecalDatum> table = tables.getTable(i);
            NestedIntegerArray<RecalDatum> mergedTable = merged.getTable(i);

            Assert.assertEquals(table.getAllLeaves().size(), mergedTable.getAllLeaves().size());
            for ( final NestedIntegerArray.Leaf<RecalDatum> leaf : table.getAllLeaves() ) {
                final RecalDatum mergedValue = mergedTable.get(leaf.keys);
                Assert.assertNotNull(mergedValue);
                Assert.assertEquals(mergedValue.getNumObservations(), leaf.value.getNumObservations() * 2);
                Assert.assertEquals(mergedValue.getNumMismatches(), leaf.value.getNumMismatches() * 2);
            }
        }
    }

    @Test
    public void testCombineWithMergeAccumulators() {
        final RecalibrationTables toMerge = new RecalibrationTables(covariates, numReadGroups);
        fillTable(toMerge);
        final RecalibrationTables toMerge2 = new RecalibrationTables(covariates, numReadGroups);
        fillTable(toMerge2);

        RecalibrationTablesMerger merger = new RecalibrationTablesMerger();
        RecalibrationTables acc = merger.createAccumulator();
        acc = merger.addInput(acc, toMerge);
        acc = merger.addInput(acc, tables);
        RecalibrationTables acc2 = merger.createAccumulator();
        acc2 = merger.addInput(acc2, toMerge2);

        RecalibrationTables mergedAccumulator = merger.mergeAccumulators(Arrays.asList(new RecalibrationTables[]{acc, acc2}));
        RecalibrationTables merged = merger.extractOutput(mergedAccumulator);

        for ( int i = 0; i < this.tables.numTables(); i++ ) {
            NestedIntegerArray<RecalDatum> table = this.tables.getTable(i);
            NestedIntegerArray<RecalDatum> mergedTable = merged.getTable(i);

            Assert.assertEquals(table.getAllLeaves().size(), mergedTable.getAllLeaves().size());
            for ( final NestedIntegerArray.Leaf<RecalDatum> leaf : table.getAllLeaves() ) {
                final RecalDatum mergedValue = mergedTable.get(leaf.keys);
                Assert.assertNotNull(mergedValue);
                Assert.assertEquals(mergedValue.getNumObservations(), leaf.value.getNumObservations() * 3);
                Assert.assertEquals(mergedValue.getNumMismatches(), leaf.value.getNumMismatches() * 3);
            }
        }
    }

    @Test
    public void testCombineEmptyOther() {
        final RecalibrationTables toMerge = new RecalibrationTables(covariates, numReadGroups);

        RecalibrationTablesMerger merger = new RecalibrationTablesMerger();
        RecalibrationTables acc = merger.createAccumulator();
        acc = merger.addInput(acc, toMerge);
        acc = merger.addInput(acc, tables);
        RecalibrationTables merged = merger.extractOutput(acc);

        for ( int i = 0; i < tables.numTables(); i++ ) {
            NestedIntegerArray<RecalDatum> table = tables.getTable(i);
            NestedIntegerArray<RecalDatum> mergedTable = merged.getTable(i);

            Assert.assertEquals(table.getAllLeaves().size(), mergedTable.getAllLeaves().size());
            for ( final NestedIntegerArray.Leaf<RecalDatum> leaf : table.getAllLeaves() ) {
                final RecalDatum mergedValue = mergedTable.get(leaf.keys);
                Assert.assertNotNull(mergedValue);
                Assert.assertEquals(mergedValue.getNumObservations(), leaf.value.getNumObservations());
                Assert.assertEquals(mergedValue.getNumMismatches(), leaf.value.getNumMismatches());
            }
        }
    }

    @Test
    public void testCombinePartial() {
        final RecalibrationTables toMerge = new RecalibrationTables(covariates, numReadGroups);
        for ( final int rg : combineStates) {
            RecalUtils.incrementDatumOrPutIfNecessary(toMerge.getTable(3), qualByte, 1, rg, 0, 0, 0);
        }

        RecalibrationTablesMerger merger = new RecalibrationTablesMerger();
        RecalibrationTables acc = merger.createAccumulator();
        acc = merger.addInput(acc, toMerge);
        acc = merger.addInput(acc, tables);
        RecalibrationTables merged = merger.extractOutput(acc);

        for ( int i = 0; i < tables.numTables(); i++ ) {
            NestedIntegerArray<RecalDatum> table = tables.getTable(i);
            NestedIntegerArray<RecalDatum> mergedTable = merged.getTable(i);

            Assert.assertEquals(table.getAllLeaves().size(), mergedTable.getAllLeaves().size());
            for ( final NestedIntegerArray.Leaf<RecalDatum> leaf : table.getAllLeaves() ) {
                final RecalDatum mergedValue = mergedTable.get(leaf.keys);
                Assert.assertNotNull(mergedValue);

                final int delta = i == 3 && leaf.keys[1] == 0 && leaf.keys[2] == 0 && leaf.keys[3] == 0 ? 1 : 0;
                Assert.assertEquals(mergedValue.getNumObservations(), leaf.value.getNumObservations() + delta);
                Assert.assertEquals(mergedValue.getNumMismatches(), leaf.value.getNumMismatches() + delta);
            }
        }
    }
}
