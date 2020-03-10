package org.broadinstitute.hellbender.utils.recalibration;

import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.engine.ReferenceMemorySource;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public final class BaseRecalibrationEngineUnitTest {

    @Test
    public void basicDBQSRFractionalErrorTestEnd() {
        byte[] baq = "@@@@@@@FGH".getBytes();
        int[] errors = new int[baq.length];
        Arrays.fill(errors, 0);
        errors[7] = 1;
        double[] answer = new double[baq.length];
        Arrays.fill(answer, 0.0);
        answer[6] = answer[7] = answer[8] = answer[9] = 1.0 / 4.0;
        double[] result = BaseRecalibrationEngine.calculateFractionalErrorArray(errors, baq);
        for( int iii = 0; iii < answer.length; iii++) {
            Assert.assertEquals(result[iii], answer[iii], 1E-6);
        }
    }

    @Test
    public void basicDBQSRFractionalErrorTestStart() {
        byte[] baq = "FFF@@@@@@@".getBytes();
        int[] errors = new int[baq.length];
        Arrays.fill(errors, 0);
        errors[2] = 1;
        double[] answer = new double[baq.length];
        Arrays.fill(answer, 0.0);
        answer[0] = answer[1] = answer[2] = answer[3] = 1.0 / 4.0;
        double[] result = BaseRecalibrationEngine.calculateFractionalErrorArray(errors, baq);
        for( int iii = 0; iii < answer.length; iii++) {
            Assert.assertEquals(result[iii], answer[iii], 1E-6);
        }
    }

    @Test
    public void basicDBQSRFractionalErrorTestNoBAQ() {
        byte[] baq = "@@@@@@@@@@".getBytes();
        int[] errors = new int[baq.length];
        Arrays.fill(errors, 0);
        errors[7] = 1;
        double[] answer = new double[baq.length];
        Arrays.fill(answer, 0.0);
        answer[7] = 1.0;
        double[] result = BaseRecalibrationEngine.calculateFractionalErrorArray(errors, baq);
        for( int iii = 0; iii < answer.length; iii++) {
            Assert.assertEquals(result[iii], answer[iii], 1E-6);
        }
    }

    @Test
    public void basicDBQSRFractionalErrorTestBAQOffset() {
        byte[] baq = "@FGH@@@@@@".getBytes();
        int[] errors = new int[baq.length];
        Arrays.fill(errors, 0);
        errors[7] = 1;
        double[] answer = new double[baq.length];
        Arrays.fill(answer, 0.0);
        answer[7] = 1.0;
        double[] result = BaseRecalibrationEngine.calculateFractionalErrorArray(errors, baq);
        for( int iii = 0; iii < answer.length; iii++) {
            Assert.assertEquals(result[iii], answer[iii], 1E-6);
        }
    }

    @Test
    public void basicDBQSRFractionalErrorTestMiddle() {
        byte[] baq = "@@@FGH@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@".getBytes();
        int[] errors = new int[baq.length];
        Arrays.fill(errors, 0);
        errors[4] = 1;
        double[] answer = new double[baq.length];
        Arrays.fill(answer, 0.0);
        answer[2] = answer[3] = answer[4] = answer[5] = answer[6] = 1.0 / 5.0;
        double[] result = BaseRecalibrationEngine.calculateFractionalErrorArray(errors, baq);
        for( int iii = 0; iii < answer.length; iii++) {
            Assert.assertEquals(result[iii], answer[iii], 1E-6);
        }
    }

    @DataProvider(name = "CalculateIsIndelData")
    public Object[][] makeCalculateIsIndelData() {
        List<Object[]> tests = new ArrayList<>();

        // this functionality can be adapted to provide input data for whatever you might want in your data
        for ( final EventType model : Arrays.asList(EventType.BASE_DELETION, EventType.BASE_INSERTION) ) {
            for ( final boolean neg : Arrays.asList(true, false) ) {
                for ( final int readLen : Arrays.asList(1, 3, 10)) {
                    tests.add(new Object[]{readLen + "M", neg, model, new int[readLen]});
                }
            }
        }

        tests.add(new Object[]{"1D3M",   false, EventType.BASE_DELETION, new int[]{0,0,0}});
        tests.add(new Object[]{"1M1D2M", false, EventType.BASE_DELETION, new int[]{1,0,0}});
        tests.add(new Object[]{"2M1D1M", false, EventType.BASE_DELETION, new int[]{0,1,0}});
        tests.add(new Object[]{"3M1D",   false, EventType.BASE_DELETION, new int[]{0,0,1}});

        tests.add(new Object[]{"1D3M",   true, EventType.BASE_DELETION, new int[]{1,0,0}});
        tests.add(new Object[]{"1M1D2M", true, EventType.BASE_DELETION, new int[]{0,1,0}});
        tests.add(new Object[]{"2M1D1M", true, EventType.BASE_DELETION, new int[]{0,0,1}});
        tests.add(new Object[]{"3M1D",   true, EventType.BASE_DELETION, new int[]{0,0,0}});

        tests.add(new Object[]{"4M1I",   false, EventType.BASE_INSERTION, new int[]{0,0,0,1,0}});
        tests.add(new Object[]{"3M1I1M", false, EventType.BASE_INSERTION, new int[]{0,0,1,0,0}});
        tests.add(new Object[]{"2M1I2M", false, EventType.BASE_INSERTION, new int[]{0,1,0,0,0}});
        tests.add(new Object[]{"1M1I3M", false, EventType.BASE_INSERTION, new int[]{1,0,0,0,0}});
        tests.add(new Object[]{"1I4M",   false, EventType.BASE_INSERTION, new int[]{0,0,0,0,0}});

        tests.add(new Object[]{"4M1I",   true, EventType.BASE_INSERTION, new int[]{0,0,0,0,0}});
        tests.add(new Object[]{"3M1I1M", true, EventType.BASE_INSERTION, new int[]{0,0,0,0,1}});
        tests.add(new Object[]{"2M1I2M", true, EventType.BASE_INSERTION, new int[]{0,0,0,1,0}});
        tests.add(new Object[]{"1M1I3M", true, EventType.BASE_INSERTION, new int[]{0,0,1,0,0}});
        tests.add(new Object[]{"1I4M",   true, EventType.BASE_INSERTION, new int[]{0,1,0,0,0}});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "CalculateIsIndelData")
    public void testCalculateIsIndel(final String cigar, final boolean negStrand, final EventType mode, final int[] expected) {
        final GATKRead read = ArtificialReadUtils.createArtificialRead(TextCigarCodec.decode(cigar));
        read.setIsReverseStrand(negStrand);
        // Fake reference data, since the indel calculation does not use the reference at all.
        final ReferenceDataSource refSource = new ReferenceMemorySource(new ReferenceBases(Utils.repeatBytes((byte)'A', read.getEnd() - read.getStart() + 1), new SimpleInterval(read)), ArtificialReadUtils.createArtificialSamHeader().getSequenceDictionary());

        int[] isSNP = new int[read.getLength()];
        int[] isInsertion = new int[isSNP.length];
        int[] isDeletion = new int[isSNP.length];
        BaseRecalibrationEngine.calculateIsSNPOrIndel(read, refSource, isSNP, isInsertion, isDeletion);
        final int[] actual = (mode == EventType.BASE_INSERTION ? isInsertion : isDeletion);
        Assert.assertEquals(actual, expected, "calculateIsSNPOrIndel() failed with " + mode + " and cigar " + cigar + " Expected " + Arrays.toString(expected) + " but got " + Arrays.toString(actual));
    }
}