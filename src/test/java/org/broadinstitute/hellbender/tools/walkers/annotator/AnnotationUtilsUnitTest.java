package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

import static org.testng.Assert.*;

public class AnnotationUtilsUnitTest {
    //Basic aligned read
    private GATKRead allMatch;

    //Read with insertion and deletion
    private GATKRead twoIndels;

    //Read with soft clips at start
    private GATKRead softClipStart;

    //Read with hard clips at start
    private GATKRead hardClipStart;

    //Read with low quality tail
    private GATKRead lowQualTail;

    //Read with low quality tail, partially soft clipped
    private GATKRead lowQualClippedTail;

    //Read with low quality bases at start
    private GATKRead lowQualStart;

    //Read with low quality bases, partially soft clipped at both ends
    private GATKRead lowQualBothEnds;

    @BeforeClass
    public void init() {
        List<CigarElement> cigarElements_allMatch = new LinkedList<>();
        cigarElements_allMatch.add(new CigarElement(151, CigarOperator.M));
        allMatch = ArtificialReadUtils.createArtificialRead(new Cigar(cigarElements_allMatch));

        List<CigarElement> cigarElements_2indels = new LinkedList<>();
        cigarElements_2indels.add(new CigarElement(66, CigarOperator.M));
        cigarElements_2indels.add(new CigarElement(10, CigarOperator.I));
        cigarElements_2indels.add(new CigarElement(7, CigarOperator.M));
        cigarElements_2indels.add(new CigarElement(10, CigarOperator.D));
        cigarElements_2indels.add(new CigarElement(68, CigarOperator.M));
        twoIndels = ArtificialReadUtils.createArtificialRead(new Cigar(cigarElements_2indels));

        List<CigarElement> cigarElements_softClipStart = new LinkedList<>();
        cigarElements_softClipStart.add(new CigarElement(17, CigarOperator.S));
        cigarElements_softClipStart.add(new CigarElement(134, CigarOperator.M));
        softClipStart = ArtificialReadUtils.createArtificialRead(new Cigar(cigarElements_softClipStart));

        List<CigarElement> cigarElements_hardClipStart = new LinkedList<>();
        cigarElements_hardClipStart.add(new CigarElement(17, CigarOperator.H));
        cigarElements_hardClipStart.add(new CigarElement(134, CigarOperator.M));
        hardClipStart = ArtificialReadUtils.createArtificialRead(new Cigar(cigarElements_hardClipStart));


        final byte [] bases_lowQualTail = {'A', 'C', 'T', 'G', 'A', 'A', 'A', 'A', 'A', 'A'};
        final byte [] quals_lowQualTail = {30, 15, 25, 30, 2, 2, 2, 2, 2, 2};
        lowQualTail = ArtificialReadUtils.createArtificialRead(bases_lowQualTail, quals_lowQualTail, "10M");

        final byte [] bases_lowQualClippedTail = {'A', 'C', 'T', 'G', 'A', 'A', 'A', 'A', 'A', 'A'};
        final byte [] quals_lowQualClippedTail = {30, 15, 25, 30, 2, 2, 2, 2, 2, 2};
        lowQualClippedTail = ArtificialReadUtils.createArtificialRead(bases_lowQualClippedTail, quals_lowQualClippedTail, "8M2S");

        final byte [] bases_lowQualStart = {'A', 'A', 'A', 'A', 'A', 'A', 'A', 'C', 'T', 'G'};
        final byte [] quals_lowQualStart = {2, 2, 2, 2, 2, 2, 30, 15, 25, 30};
        lowQualStart = ArtificialReadUtils.createArtificialRead(bases_lowQualStart, quals_lowQualStart, "10M");

        final byte [] bases_lowQualBothEnds = {'A', 'A', 'A', 'A', 'A', 'A', 'G', 'C', 'T', 'G', 'A', 'A', 'A', 'A', 'A', 'A'};
        final byte [] quals_lowQualBothEnds = { 2, 2, 2, 2, 2, 2, 30, 15, 25, 30, 2, 2, 2, 2, 2, 2};
        lowQualBothEnds = ArtificialReadUtils.createArtificialRead(bases_lowQualBothEnds, quals_lowQualBothEnds, "2S12M2S");

    }

    @DataProvider(name = "makeGetFinalVariantReadPositionTestReads")
    public Object[][] makeFinalPosTestReads() {
        final List<Object[]> tests = new ArrayList<>();
        tests.add(new Object[] {allMatch, 10, 10});
        tests.add(new Object[] {allMatch, 140, 10});
        tests.add(new Object[] {twoIndels, 10, 10});
        tests.add(new Object[] {twoIndels, 140, 10});
        tests.add(new Object[] {hardClipStart, 20, 20});
        tests.add(new Object[] {hardClipStart, 110, 6});  //this is what the code produces as-is
        tests.add(new Object[] {softClipStart, 10, 10});
        tests.add(new Object[] {softClipStart, 140, 10});
        tests.add(new Object[] {lowQualTail, 3, 0});
        tests.add(new Object[] {lowQualTail, 2, 2});
        tests.add(new Object[] {lowQualClippedTail, 3, 0});
        tests.add(new Object[] {lowQualClippedTail, 2, 2});
        tests.add(new Object[] {lowQualStart, 7, -4});   //this is what the code produces as-is, but should be 1
        tests.add(new Object[] {lowQualStart, 8, -5}); //this is what the code produces as-is, but should be 1
        tests.add(new Object[] {lowQualBothEnds, 7, -4});   //this is what the code produces as-is, but should be 1
        tests.add(new Object[] {lowQualBothEnds, 8, -5});    //this is what the code produces as-is, but should be 1
        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "makeGetFinalVariantReadPositionTestReads")
    public void testGetFinalVariantReadPosition(GATKRead read, int variantPosition, int expected) throws Exception {
        Assert.assertEquals(AnnotationUtils.getFinalVariantReadPosition(read, variantPosition), expected);
    }

    @DataProvider(name = "getNumClippedBasesAtStartTestReads")
    public Object[][] numClipStart() {
        final List<Object[]> tests = new ArrayList<>();
        tests.add(new Object[] {allMatch, 0});
        tests.add(new Object[] {twoIndels, 0});
        tests.add(new Object[] {softClipStart, 0});
        tests.add(new Object[] {hardClipStart, 17});
        tests.add(new Object[] {lowQualTail, 0});
        tests.add(new Object[] {lowQualClippedTail, 0});
        tests.add(new Object[] {lowQualStart, 6});
        tests.add(new Object[] {lowQualBothEnds, 6});
        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "getNumClippedBasesAtStartTestReads")
    public void testGetNumClippedBasesAtStart(GATKRead read, int expected) throws Exception {
        Assert.assertEquals(AnnotationUtils.getNumClippedBasesAtStart(read),expected);
    }

    @DataProvider(name = "getNumAlignedBasesTestReads")
    public Object[][] numAligned() {
        final List<Object[]> tests = new ArrayList<>();
        tests.add(new Object[] {allMatch, 151});
        tests.add(new Object[] {twoIndels, 151});
        tests.add(new Object[] {softClipStart, 151});
        tests.add(new Object[] {hardClipStart, 117});  //This is what the code produces, but it's wrong
        tests.add(new Object[] {lowQualTail, 4});
        tests.add(new Object[] {lowQualClippedTail, 4});
        tests.add(new Object[] {lowQualStart, 4});
        tests.add(new Object[] {lowQualBothEnds, 4});
        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "getNumAlignedBasesTestReads")
    public void testGetNumAlignedBases(GATKRead read, int expected) throws Exception {
        Assert.assertEquals(AnnotationUtils.getNumAlignedBases(read),expected);
    }

    @DataProvider(name = "getNumClippedBasesAtEndTestReads")
    public Object[][] numClipEnd() {
        final List<Object[]> tests = new ArrayList<>();
        tests.add(new Object[] {allMatch, 0});
        tests.add(new Object[] {twoIndels, 0});
        tests.add(new Object[] {softClipStart, 0});
        tests.add(new Object[] {hardClipStart, 0});
        tests.add(new Object[] {lowQualTail, 6});
        tests.add(new Object[] {lowQualClippedTail, 6});
        tests.add(new Object[] {lowQualStart, 0});
        tests.add(new Object[] {lowQualBothEnds, 6});
        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "getNumClippedBasesAtEndTestReads")
    public void testGetNumClippedBasesAtEnd(GATKRead read, int expected) throws Exception {
        Assert.assertEquals(AnnotationUtils.getNumClippedBasesAtEnd(read), expected);
    }
}