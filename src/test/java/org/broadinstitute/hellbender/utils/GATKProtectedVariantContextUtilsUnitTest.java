package org.broadinstitute.hellbender.utils;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;

/**
 * Created by David Benjamin on 2/15/17.
 */
public class GATKProtectedVariantContextUtilsUnitTest {
    @Test
    public void testGetPileup() {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader(1, 1, 1000);
        final Locatable loc = new SimpleInterval("chr1", 10, 10);
        final int readLength = 3;

        //this read doesn't overlap {@code loc}
        final GATKRead read1 = ArtificialReadUtils.createArtificialRead(header, "read1", 0, 1, readLength);
        read1.setBases(Utils.dupBytes((byte) 'A', readLength));
        read1.setBaseQualities(Utils.dupBytes((byte) 30, readLength));
        read1.setCigar("3M");

        //this read overlaps {@code loc} with a Q30 'A'
        final GATKRead read2 = ArtificialReadUtils.createArtificialRead(header, "read2", 0, 10, readLength);
        read2.setBases(Utils.dupBytes((byte) 'A', readLength));
        read2.setBaseQualities(Utils.dupBytes((byte) 30, readLength));
        read2.setCigar("3M");

        //this read overlaps {@code loc} with a Q20 'C' (due to the deletion)
        final GATKRead read3 = ArtificialReadUtils.createArtificialRead(header, "read3", 0, 7, readLength);
        read3.setBases(Utils.dupBytes((byte) 'C', readLength));
        read3.setBaseQualities(Utils.dupBytes((byte) 20, readLength));
        read3.setCigar("1M1D2M");

        //this read doesn't overlap {@code loc} due to the deletion
        final GATKRead read4 = ArtificialReadUtils.createArtificialRead(header, "read4", 0, 7, readLength);
        read4.setBases(Utils.dupBytes((byte) 'C', readLength));
        read4.setBaseQualities(Utils.dupBytes((byte) 20, readLength));
        read4.setCigar("1M5D2M");

        final ReadPileup pileup = GATKProtectedVariantContextUtils.getPileup(loc, Arrays.asList(read1, read2, read3, read4));

        // the pileup should contain a Q30 'A' and a Q20 'C'
        final int[] counts = pileup.getBaseCounts();
        Assert.assertEquals(counts, new int[]{1, 1, 0, 0});

    }

    @Test(dataProvider = "variantTypes")
    public void testVariantTypesAndIsComplex(final String ref, final String alt, final VariantContext.Type gtType, boolean isComplexIndel) {
        Assert.assertEquals(GATKProtectedVariantContextUtils.typeOfVariant(Allele.create(ref), Allele.create(alt)), gtType);
        Assert.assertEquals(GATKProtectedVariantContextUtils.isComplexIndel(Allele.create(ref), Allele.create(alt)), isComplexIndel);
    }
    @Test(expectedExceptions = IllegalStateException.class)
    public void testSymbolicRef() {
        GATKProtectedVariantContextUtils.typeOfVariant(GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE, Allele.create("C"));
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

}