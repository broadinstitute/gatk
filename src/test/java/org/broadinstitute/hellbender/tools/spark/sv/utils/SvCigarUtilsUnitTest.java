package org.broadinstitute.hellbender.tools.spark.sv.utils;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.stream.Stream;


public class SvCigarUtilsUnitTest extends GATKBaseTest {

    @Test(expectedExceptions=IllegalArgumentException.class, groups = "sv")
    public void testCigarChecker_emptyCigarElementList(){
        final List<CigarElement> emptyList = Collections.emptyList();
        SvCigarUtils.validateCigar(emptyList);
    }

    @Test(expectedExceptions=IllegalArgumentException.class, groups = "sv")
    public void testCigarChecker_deletionNeighboringClipping(){
        SvCigarUtils.validateCigar(TextCigarCodec.decode("10S10D10M").getCigarElements());
        SvCigarUtils.validateCigar(TextCigarCodec.decode("10H10S10D10M").getCigarElements());
        SvCigarUtils.validateCigar(TextCigarCodec.decode("10M10D10S").getCigarElements());
        SvCigarUtils.validateCigar(TextCigarCodec.decode("10M10D10S10H").getCigarElements());
    }

    @Test(expectedExceptions=IllegalArgumentException.class, groups = "sv")
    public void testCigarChecker_only1NonAlignment(){
        SvCigarUtils.validateCigar(TextCigarCodec.decode("10S").getCigarElements());
    }

    @Test(expectedExceptions=IllegalArgumentException.class, groups = "sv")
    public void testCigarChecker_noAlignment(){
        SvCigarUtils.validateCigar(TextCigarCodec.decode("10H10S10I10S10H").getCigarElements());
    }

    @Test(groups = "sv")
    public void testGetIndexOfFirstNonClippingBase(){
        Assert.assertEquals(SvCigarUtils.findIndexOfFirstNonClippingOperation(TextCigarCodec.decode("151M").getCigarElements(), true), 0);
        Assert.assertEquals(SvCigarUtils.findIndexOfFirstNonClippingOperation(TextCigarCodec.decode("151M").getCigarElements(), false), 0);
        Assert.assertEquals(SvCigarUtils.findIndexOfFirstNonClippingOperation(TextCigarCodec.decode("10S10D10M").getCigarElements(), true), 1);
        Assert.assertEquals(SvCigarUtils.findIndexOfFirstNonClippingOperation(TextCigarCodec.decode("10H10S10D10M").getCigarElements(), true), 2);
        Assert.assertEquals(SvCigarUtils.findIndexOfFirstNonClippingOperation(TextCigarCodec.decode("10M10D10S").getCigarElements(), false), 1);
        Assert.assertEquals(SvCigarUtils.findIndexOfFirstNonClippingOperation(TextCigarCodec.decode("10M10D10S10H").getCigarElements(), false), 1);
    }

}
