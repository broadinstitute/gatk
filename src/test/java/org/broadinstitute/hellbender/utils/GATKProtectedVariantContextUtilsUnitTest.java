package org.broadinstitute.hellbender.utils;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * Created by David Benjamin on 2/15/17.
 */
public class GATKProtectedVariantContextUtilsUnitTest extends GATKBaseTest {

    public static final int LENGTH_OF_ARTIFICIAL_READ = 50;

    @Test(dataProvider = "variantTypes")
    public void testVariantTypesAndIsComplex(final String ref, final String alt, final VariantContext.Type gtType, boolean isComplexIndel) {
        Assert.assertEquals(GATKProtectedVariantContextUtils.typeOfVariant(Allele.create(ref), Allele.create(alt)), gtType);
        Assert.assertEquals(GATKProtectedVariantContextUtils.isComplexIndel(Allele.create(ref), Allele.create(alt)), isComplexIndel);
    }
    @Test(expectedExceptions = IllegalStateException.class)
    public void testSymbolicRef() {
        GATKProtectedVariantContextUtils.typeOfVariant(Allele.NON_REF_ALLELE, Allele.create("C"));
    }

    @DataProvider(name = "variantTypes")
    public Object[][] variantTypes() {
        return new Object[][]{
                // ref, alt, type, isComplex?
                {"CCTTGGCTTATTCCA", "C", VariantContext.Type.INDEL, false},
                {"C", "CCTTGGCTTATTCCA", VariantContext.Type.INDEL, false},
                {"ACTAG", "A", VariantContext.Type.INDEL, false},
                {"ATT", "AT", VariantContext.Type.INDEL, false},
                {"AT", "ATT", VariantContext.Type.INDEL, false},
                {"CT", "CAGG", VariantContext.Type.INDEL, true},
                {"CTTT", "CAGG", VariantContext.Type.MNP, false},
                {"CTTT", "CAGGG", VariantContext.Type.INDEL, true},
                {"T", "T", VariantContext.Type.NO_VARIATION, false},
                {"CTAG", "CTAG", VariantContext.Type.NO_VARIATION, false},
                {"A", "AAGAAGCATGC", VariantContext.Type.INDEL, false},
                {"A", "C", VariantContext.Type.SNP, false},
                {"AG", "CA", VariantContext.Type.MNP, false},
                {"AGAAGG", "CATTCC", VariantContext.Type.MNP, false},
                {"GC", "GA", VariantContext.Type.SNP, false},
                {"GA", "<NON_REF>", VariantContext.Type.SYMBOLIC, false},
                {"GA", "*", VariantContext.Type.NO_VARIATION, false},

                // There are two MNPs here
                {"AGAAGG", "CATACC", VariantContext.Type.MNP, false},

                // Note that this is technically a simple AT insertion, but the isComplex cannot handle this properly.
                {"CT", "CATT", VariantContext.Type.INDEL, true},
        };
    }

    @Test(dataProvider = "bestAlleleSNP")
    public void testChooseBestAlleleSNP(final Allele referenceAllele, final List<Allele> altAlleles, int offsetIntoRead, int minBaseQualityCutoff, int gtChosenInAlt, boolean isGtRef) {

        // No indels
        if (isGtRef) {
            final PileupElement pileupElement = ArtificialReadUtils.createNonIndelPileupElement(offsetIntoRead, referenceAllele, LENGTH_OF_ARTIFICIAL_READ);
            final Allele chosen = GATKProtectedVariantContextUtils.chooseAlleleForRead(pileupElement, referenceAllele, altAlleles, minBaseQualityCutoff);
            Assert.assertEquals(chosen, referenceAllele);
        } else {
            final PileupElement pileupElement = ArtificialReadUtils.createNonIndelPileupElement(offsetIntoRead, altAlleles.get(gtChosenInAlt), LENGTH_OF_ARTIFICIAL_READ);
            final Allele chosen = GATKProtectedVariantContextUtils.chooseAlleleForRead(pileupElement, referenceAllele, altAlleles, minBaseQualityCutoff);
            Assert.assertEquals(chosen, altAlleles.get(gtChosenInAlt));
        }
    }

    @DataProvider(name = "bestAlleleSNP")
    public Object[][] createSNPAlleles() {
        return new Object[][] {
                // final Allele referenceAllele, final List<Allele> altAlleles, int offsetIntoRead, int minBaseQualityCutoff, int gtChosenInAlt, boolean isGtRef
                {
                    Allele.create("A", true), Arrays.asList(Allele.create("T", false), Allele.create("CA", false)), 0, 0,
                        0, false
                }, {
                    Allele.create("A", true), Arrays.asList(Allele.create("T", false), Allele.create("CA", false)), 0,  0,
                        0, true
                }, {
                    Allele.create("A", true), Arrays.asList(Allele.create("C", false), Allele.create("T", false)), 1, 0,
                        0, false
                }, {
                    Allele.create("AA", true), Arrays.asList(Allele.create("CC", false), Allele.create("T", false)), 1, 0,
                        0, false
                }, {
                    Allele.create("AA", true), Arrays.asList(Allele.create("T", false), Allele.create("CC", false)), 1, 0,
                        1, false
                }, {
                    Allele.create("GA", true), Arrays.asList(Allele.create("T", false), Allele.create("CC", false)), 1, 0,
                    1, true
                }
        };
    }

    @Test(dataProvider = "testInsertions")
    public void testChooseBestAlleleInsertion(final int offsetIntoRead, final String refBases, final String altBases,
                                              final boolean isSpliceInAlt) {

        // We assume that we are on the base right before the insertion.  We are pretending each read is 50 bases long.
        final Allele referenceAllele = Allele.create(refBases, true);
        final Allele insertionAllele = Allele.create(altBases, false);

        if (isSpliceInAlt) {
            final PileupElement pileupElement = ArtificialReadUtils.createSplicedInsertionPileupElement(offsetIntoRead, insertionAllele, LENGTH_OF_ARTIFICIAL_READ);
            Assert.assertEquals(
                    GATKProtectedVariantContextUtils.chooseAlleleForRead(pileupElement, referenceAllele,
                            Collections.singletonList(insertionAllele), 0),
                    insertionAllele);
        } else {
            final PileupElement pileupElement = ArtificialReadUtils.createNonIndelPileupElement(offsetIntoRead, referenceAllele, LENGTH_OF_ARTIFICIAL_READ);
            Assert.assertEquals(
                    GATKProtectedVariantContextUtils.chooseAlleleForRead(pileupElement, referenceAllele,
                            Collections.singletonList(insertionAllele), 0),
                    referenceAllele);
        }
    }

    @Test(dataProvider = "testDeletions")
    public void testChooseBestAlleleDeletion(final int offsetIntoRead, final String refBases, final String altBases,
                                              final boolean isSpliceInAlt) {

        // We assume that we are on the base right before the insertion.  We are pretending each read is 50 bases long.
        final Allele referenceAllele = Allele.create(refBases, true);
        final Allele deletionAllele = Allele.create(altBases, false);

        if (isSpliceInAlt) {
            final PileupElement pileupElement = ArtificialReadUtils.createSplicedDeletionPileupElement(offsetIntoRead, referenceAllele, LENGTH_OF_ARTIFICIAL_READ);
            Assert.assertEquals(
                    GATKProtectedVariantContextUtils.chooseAlleleForRead(pileupElement, referenceAllele,
                            Collections.singletonList(deletionAllele), 0),
                    deletionAllele);
        } else {
            final PileupElement pileupElement = ArtificialReadUtils.createNonIndelPileupElement(offsetIntoRead, referenceAllele, LENGTH_OF_ARTIFICIAL_READ);
            Assert.assertEquals(
                    GATKProtectedVariantContextUtils.chooseAlleleForRead(pileupElement, referenceAllele,
                            Collections.singletonList(deletionAllele), 0),
                    referenceAllele);
        }
    }

    @Test
    public void testChooseBestAlleleNull() {

        final int offsetIntoRead = 10;

        // We assume that we are on the base right before the insertion.  We are pretending each read is 50 bases long.
        final Allele referenceAllele = Allele.create("T", true);
        final Allele deletionAllele = Allele.create("A", false);

        final PileupElement pileupElement = ArtificialReadUtils.createNonIndelPileupElement(offsetIntoRead, referenceAllele, LENGTH_OF_ARTIFICIAL_READ);
        final byte[] newBases = pileupElement.getRead().getBases();
        newBases[offsetIntoRead] = (byte) 'C';
        pileupElement.getRead().setBases(newBases);

        Assert.assertEquals(
                GATKProtectedVariantContextUtils.chooseAlleleForRead(pileupElement, referenceAllele,
                        Collections.singletonList(deletionAllele), 0),
                null);
    }

    @DataProvider(name="testInsertions")
    public Object[][] createTestInsertions() {
        return new Object[][] {
                {0, "A", "ATT", true},
                {0, "A", "ATT", false},
                {10, "A", "ATT", true},
                {10, "A", "ATT", false}
        };
    }
    @DataProvider(name="testDeletions")
    public Object[][] createTestDeletions() {
        return new Object[][] {
                {0, "ATT", "A", true},
                {0, "ATT", "A", false},
                {10, "ATT", "A", true},
                {10, "ATT", "A", false}
        };
    }

    @Test(dataProvider = "doesReadContainAllele")
    public void testDoesReadContainAllele(final int offsetIntoRead, final String actualBases, final String searchAlleleBases,
                                          final Trilean gt) {

        final PileupElement pileupElement = ArtificialReadUtils.createNonIndelPileupElement(offsetIntoRead, Allele.create(actualBases, true), LENGTH_OF_ARTIFICIAL_READ);
        Assert.assertEquals(GATKProtectedVariantContextUtils.doesReadContainAllele(pileupElement, Allele.create(searchAlleleBases, false)),
                gt);
        Assert.assertEquals(GATKProtectedVariantContextUtils.doesReadContainAllele(pileupElement, Allele.create(actualBases, false)),
                Trilean.TRUE);
    }

    @DataProvider(name="doesReadContainAllele")
    public Object[][] createDoesReadContainAlelle() {
        return new Object[][] {
                {10, "ATT", "C", Trilean.FALSE},
                {10, "AT", "AT", Trilean.TRUE},
                {49, "A", "ATT", Trilean.UNKNOWN},
                {48, "AT", "ATT", Trilean.UNKNOWN},
        };
    }
}