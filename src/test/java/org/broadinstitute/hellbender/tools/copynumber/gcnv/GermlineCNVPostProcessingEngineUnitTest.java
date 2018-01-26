package org.broadinstitute.hellbender.tools.copynumber.gcnv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CopyNumberPosteriorDistribution;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.LocatableCopyNumberPosteriorDistribution;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Unit tests for {@link GermlineCNVPostProcessingEngine}
 */
public class GermlineCNVPostProcessingEngineUnitTest extends CommandLineProgramTest {

    //objects common to more than a single test in this class
    private final static SimpleInterval testInterval = new SimpleInterval("1", 1, 10000);
    private final static List<IntegerCopyNumberState> testCopyNumberStateList = new ArrayList<>();
    private final static IntegerCopyNumberStateCollection testCopyNumberStateCollection;

    static {
        testCopyNumberStateList.add(new IntegerCopyNumberState(0));
        testCopyNumberStateList.add(new IntegerCopyNumberState(1));
        testCopyNumberStateList.add(new IntegerCopyNumberState(2));
        testCopyNumberStateList.add(new IntegerCopyNumberState(3));
        testCopyNumberStateList.add(new IntegerCopyNumberState(4));
        testCopyNumberStateList.add(new IntegerCopyNumberState(5));
        testCopyNumberStateCollection =
                new IntegerCopyNumberStateCollection(testCopyNumberStateList.stream().map(s -> s.toString()).collect(Collectors.toList()));
    }

    @DataProvider(name = "examplePosteriorRecords")
    public Object[][] testData() {
        final Map<IntegerCopyNumberState, Double> copyNumberPosteriorDistribution = new HashMap<>();
        copyNumberPosteriorDistribution.put(testCopyNumberStateList.get(0), -22.133276724683618);
        copyNumberPosteriorDistribution.put(testCopyNumberStateList.get(1), -4.3825075766712871);
        copyNumberPosteriorDistribution.put(testCopyNumberStateList.get(2), -0.012754265709645551);
        copyNumberPosteriorDistribution.put(testCopyNumberStateList.get(3), -8.6265377789688955);
        copyNumberPosteriorDistribution.put(testCopyNumberStateList.get(4), -21.918647174298602);
        copyNumberPosteriorDistribution.put(testCopyNumberStateList.get(5), -22.133276474165779);
        final CopyNumberPosteriorDistribution posteriorRecord = new CopyNumberPosteriorDistribution(copyNumberPosteriorDistribution);
        final LocatableCopyNumberPosteriorDistribution locatablePosteriorRecord = new LocatableCopyNumberPosteriorDistribution(testInterval, posteriorRecord);
        final List<Integer> PLs = new ArrayList<>(
                Arrays.asList(new Integer[] {96, 18, 0, 37, 95, 96}));
        final int expectedMAPCopyNumber = 2;
        final int expectedGQ = 18;

        return new Object[][] {
            {locatablePosteriorRecord, testCopyNumberStateCollection, PLs, expectedMAPCopyNumber, expectedGQ}
        };
    }

    //Test that the exception is thrown if normalized posterior probabilities do not add up to 1
    //TODO move this to the CopyNumberPosteriorRecordUnitTest class
    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testPosteriorValidation() {
        final Map<IntegerCopyNumberState, Double> copyNumberPosteriorDistribution = new HashMap<>();
        copyNumberPosteriorDistribution.put(testCopyNumberStateList.get(0), -22.133276724683618);
        copyNumberPosteriorDistribution.put(testCopyNumberStateList.get(1), -4.3825075766712871);
        copyNumberPosteriorDistribution.put(testCopyNumberStateList.get(2), -0.112754265709645551);
        copyNumberPosteriorDistribution.put(testCopyNumberStateList.get(3), -8.6265377789688955);
        copyNumberPosteriorDistribution.put(testCopyNumberStateList.get(4), -21.918647174298602);
        copyNumberPosteriorDistribution.put(testCopyNumberStateList.get(5), -22.133276474165779);
        final CopyNumberPosteriorDistribution posteriorRecord = new CopyNumberPosteriorDistribution(copyNumberPosteriorDistribution);
        final LocatableCopyNumberPosteriorDistribution locatablePosteriorRecord = new LocatableCopyNumberPosteriorDistribution(testInterval, posteriorRecord);
        test(locatablePosteriorRecord, testCopyNumberStateCollection, null, 0, 0);
    }

    //Test that exception is thrown when a copy number state collection has less than 3 states
    @Test(expectedExceptions = IllegalStateException.class)
    public void testNumberOfCopyStatesValidation() {
        final File outputFile = createTempFile("test", ".tsv");
        final String testSampleName = "sampleName";
        final SAMSequenceDictionary testSAMSequenceDictionary = new SAMSequenceDictionary();
        final VariantContextWriter outputWriter = GATKVariantContextUtils.createVCFWriter(outputFile, null, false);
        final List<IntegerCopyNumberState> shortCopyNumberStateList = new ArrayList<>(
                Arrays.asList(new IntegerCopyNumberState[] {new IntegerCopyNumberState(0), new IntegerCopyNumberState(1)})
        );
        final IntegerCopyNumberStateCollection shortCopyNumberStateCollection =
                new IntegerCopyNumberStateCollection(shortCopyNumberStateList.stream().map(s -> s.toString())
                        .collect(Collectors.toList()));
        new GermlineCNVPostProcessingEngine(outputWriter, shortCopyNumberStateCollection, testSampleName, testSAMSequenceDictionary);
    }

    @Test(dataProvider = "examplePosteriorRecords")
    public void test(final LocatableCopyNumberPosteriorDistribution posteriorLocatableRecord,
                     final IntegerCopyNumberStateCollection copyNumberStateCollection,
                     final List<Integer> expectedPLs,
                     final int expectedMAPCopyNumber,
                     final int expectedGQ) {
        final List<Integer> actualPLs = GermlineCNVPostProcessingEngine.getCopyNumberPLVector(posteriorLocatableRecord, copyNumberStateCollection);
        final int actualMAPCopyNumber = GermlineCNVPostProcessingEngine.calculateMAPCopyNumberState(posteriorLocatableRecord, copyNumberStateCollection);
        final int actualGQ = GermlineCNVPostProcessingEngine.calculateGenotypeQuality(posteriorLocatableRecord, copyNumberStateCollection);
        Assert.assertEquals(actualMAPCopyNumber, expectedMAPCopyNumber);
        IntStream.range(0, expectedPLs.size())
                .forEach(i -> Assert.assertEquals(actualPLs.get(i).intValue(), expectedPLs.get(i).intValue()));
        Assert.assertEquals(actualGQ, expectedGQ);
    }

}