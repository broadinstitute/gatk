package org.broadinstitute.hellbender.utils.read;

import com.google.api.services.genomics.model.CigarUnit;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class CigarConversionUtilsUnitTest extends BaseTest {

    private CigarUnit createCigarUnit( final Long operationLength, final String operation ) {
        final CigarUnit unit = new CigarUnit();
        unit.setOperationLength(operationLength);
        unit.setOperation(operation);
        return unit;
    }

    @DataProvider(name = "MatchingCigarUnitsAndElementsData")
    public Object[][] getMatchingCigarUnitsAndElementsData() {
        return new Object[][] {
                { createCigarUnit(1l, "ALIGNMENT_MATCH"), new CigarElement(1, CigarOperator.M) },
                { createCigarUnit(1l, "INSERT"), new CigarElement(1, CigarOperator.I) },
                { createCigarUnit(1l, "DELETE"), new CigarElement(1, CigarOperator.D) },
                { createCigarUnit(1l, "SKIP"), new CigarElement(1, CigarOperator.N) },
                { createCigarUnit(1l, "CLIP_SOFT"), new CigarElement(1, CigarOperator.S) },
                { createCigarUnit(1l, "CLIP_HARD"), new CigarElement(1, CigarOperator.H) },
                { createCigarUnit(1l, "PAD"), new CigarElement(1, CigarOperator.P) },
                { createCigarUnit(1l, "SEQUENCE_MATCH"), new CigarElement(1, CigarOperator.EQ) },
                { createCigarUnit(1l, "SEQUENCE_MISMATCH"), new CigarElement(1, CigarOperator.X) }
        };
    }

    @Test(dataProvider = "MatchingCigarUnitsAndElementsData")
    public void testConvertCigarUnitToSAMCigarElement( final CigarUnit cigarUnit, final CigarElement cigarElement ) {
        final CigarElement convertedElement = CigarConversionUtils.convertCigarUnitToSAMCigarElement(cigarUnit);
        Assert.assertEquals(convertedElement, cigarElement, "CigarUnit -> CigarElement conversion failed");
    }

    @Test(dataProvider = "MatchingCigarUnitsAndElementsData")
    public void testConvertSAMCigarElementToCigarUnit( final CigarUnit cigarUnit, final CigarElement cigarElement ) {
        final CigarUnit convertedUnit = CigarConversionUtils.convertSAMCigarElementToCigarUnit(cigarElement);
        Assert.assertEquals(convertedUnit.getOperationLength(), cigarUnit.getOperationLength(), "CigarElement -> CigarUnit conversion failed: wrong operation length");
        Assert.assertEquals(convertedUnit.getOperation(), cigarUnit.getOperation(), "CigarElement -> CigarUnit conversion failed: wrong operation");
    }

    @DataProvider(name = "InvalidCigarUnitData")
    public Object[][] getInvalidCigarUnitData() {
        return new Object[][] {
                { null, IllegalArgumentException.class },
                { createCigarUnit(null, "ALIGNMENT_MATCH"), IllegalArgumentException.class },
                { createCigarUnit(1l, null), IllegalArgumentException.class },
                { createCigarUnit(1l, "NON_EXISTENT_OPERATION"), GATKException.class },
                { createCigarUnit(-1l, "ALIGNMENT_MATCH"), GATKException.class },
                { createCigarUnit(1l + Integer.MAX_VALUE, "ALIGNMENT_MATCH"), GATKException.class }
        };
    }

    @Test(dataProvider = "InvalidCigarUnitData")
    public void testHandleInvalidCigarUnit( final CigarUnit invalidCigarUnit, final Class<? extends Exception> exceptionClass ) {
        boolean sawException = false;

        try {
            final CigarElement convertedElement = CigarConversionUtils.convertCigarUnitToSAMCigarElement(invalidCigarUnit);
        }
        catch ( Exception e ) {
            sawException = true;
            Assert.assertEquals(e.getClass(), exceptionClass, "Wrong class of exception thrown for invalid CigarUnit");
        }

        Assert.assertTrue(sawException, "Expected exception of class " + exceptionClass + ", but no exception was thrown");
    }

    @DataProvider(name = "InvalidCigarElementData")
    public Object[][] getInvalidCigarElementData() {
        return new Object[][] {
                { null, IllegalArgumentException.class },
                { new CigarElement(1, null), GATKException.class },
                { new CigarElement(-1, CigarOperator.M), GATKException.class }
        };
    }

    @Test(dataProvider = "InvalidCigarElementData")
    public void testHandleInvalidCigarElement( final CigarElement invalidCigarElement, final Class<? extends Exception> exceptionClass ) {
        boolean sawException = false;

        try {
            final CigarUnit convertedUnit = CigarConversionUtils.convertSAMCigarElementToCigarUnit(invalidCigarElement);
        }
        catch ( Exception e ) {
            sawException = true;
            Assert.assertEquals(e.getClass(), exceptionClass, "Wrong class of exception thrown for invalid CigarElement");
        }

        Assert.assertTrue(sawException, "Expected exception of class " + exceptionClass + ", but no exception was thrown");
    }

    @DataProvider(name = "MatchingSAMCigarAndCigarUnitListData")
    public Object[][] getMatchingSAMCigarAndCigarUnitListData() {
        return new Object[][] {
                { Collections.<CigarUnit>emptyList(), new Cigar() },
                { Arrays.asList(createCigarUnit(1l, "ALIGNMENT_MATCH")), TextCigarCodec.decode("1M") },
                { Arrays.asList(createCigarUnit(1l, "ALIGNMENT_MATCH"), createCigarUnit(2l, "INSERT")), TextCigarCodec.decode("1M2I") },
                { Arrays.asList(createCigarUnit(1l, "ALIGNMENT_MATCH"), createCigarUnit(2l, "INSERT"), createCigarUnit(3l, "CLIP_HARD")), TextCigarCodec.decode("1M2I3H") },
        };
    }

    @Test(dataProvider = "MatchingSAMCigarAndCigarUnitListData")
    public void testConvertCigarUnitListToSAMCigar( final List<CigarUnit> cigarUnitList, final Cigar samCigar ) {
        final Cigar convertedCigar = CigarConversionUtils.convertCigarUnitListToSAMCigar(cigarUnitList);
        Assert.assertEquals(convertedCigar, samCigar, "Wrong Cigar after CigarUnit list -> SAM Cigar conversion");
    }

    @Test(dataProvider = "MatchingSAMCigarAndCigarUnitListData")
    public void testConvertSAMCigarToCigarUnitList( final List<CigarUnit> cigarUnitList, final Cigar samCigar ) {
        final List<CigarUnit> convertedCigar = CigarConversionUtils.convertSAMCigarToCigarUnitList(samCigar);

        Assert.assertEquals(convertedCigar.size(), cigarUnitList.size(), "Wrong number of elements in CigarUnit list after conversion from SAM cigar");
        for ( int i = 0; i < convertedCigar.size(); ++i ) {
            final CigarUnit actualCigarUnit = convertedCigar.get(i);
            final CigarUnit expectedCigarUnit = cigarUnitList.get(i);

            Assert.assertEquals(actualCigarUnit.getOperationLength(), expectedCigarUnit.getOperationLength(), "Wrong operation length in CigarUnit after conversion from SAM Cigar");
            Assert.assertEquals(actualCigarUnit.getOperation(), expectedCigarUnit.getOperation(), "Wrong operation in CigarUnit after conversion from SAM Cigar");
        }
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testHandleNullCigarUnitList() {
        CigarConversionUtils.convertCigarUnitListToSAMCigar(null);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testHandleNullCigar() {
        CigarConversionUtils.convertSAMCigarToCigarUnitList(null);
    }
}
