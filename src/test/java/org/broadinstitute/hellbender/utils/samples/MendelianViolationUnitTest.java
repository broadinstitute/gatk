package org.broadinstitute.hellbender.utils.samples;

import htsjdk.variant.variantcontext.*;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.junit.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

public final class MendelianViolationUnitTest extends BaseTest {

    private static final int GQ30 = 30;

    @DataProvider(name = "isViolation")
    public Object[][] isViolation() {
        final List<Object[]> tests = new ArrayList<>();

        final Allele refAllele = Allele.create("A", true);
        final Allele altAllele = Allele.create("T");
        final Allele noCallAllele = Allele.NO_CALL;

        final List<Allele> homRef = Arrays.asList(refAllele, refAllele);
        final List<Allele> homVar = Arrays.asList(altAllele, altAllele);
        final List<Allele> het = Arrays.asList(refAllele, altAllele);
        final List<Allele> notCall = Arrays.asList(noCallAllele, noCallAllele);

        final Genotype g00 = new GenotypeBuilder("2", homRef).DP(10).AD(new int[]{10,0}).GQ(GQ30).make();
        final Genotype g01 = new GenotypeBuilder("1", het).DP(10).AD(new int[]{5, 5}).GQ(GQ30).make();
        final Genotype g11 = new GenotypeBuilder("3", homVar).DP(10).AD(new int[]{0, 10}).GQ(GQ30).make();
        final Genotype gNo = new GenotypeBuilder("4", notCall).make();

        tests.add(new Object[]{g00, g00, g00, false});
        tests.add(new Object[]{g00, g00, g01, true});
        tests.add(new Object[]{g00, g00, g11, true});
        tests.add(new Object[]{g00, g00, gNo, false});

        tests.add(new Object[]{g00, g01, g00, false});
        tests.add(new Object[]{g00, g01, g01, false});
        tests.add(new Object[]{g00, g01, g11, true});
        tests.add(new Object[]{g00, g01, gNo, false});

        tests.add(new Object[]{g00, g11, g00, true});
        tests.add(new Object[]{g00, g11, g01, false});
        tests.add(new Object[]{g00, g11, g11, true});
        tests.add(new Object[]{g00, g11, gNo, false});

        tests.add(new Object[]{g00, gNo, g00, false});
        tests.add(new Object[]{g00, gNo, g01, false});
        tests.add(new Object[]{g00, gNo, g11, true});
        tests.add(new Object[]{g00, gNo, gNo, false});


        tests.add(new Object[]{g01, g00, g00, false});
        tests.add(new Object[]{g01, g00, g01, false});
        tests.add(new Object[]{g01, g00, g11, true});
        tests.add(new Object[]{g01, g00, gNo, false});

        tests.add(new Object[]{g01, g01, g00, false});
        tests.add(new Object[]{g01, g01, g01, false});
        tests.add(new Object[]{g01, g01, g11, false});
        tests.add(new Object[]{g01, g01, gNo, false});

        tests.add(new Object[]{g01, g11, g00, true});
        tests.add(new Object[]{g01, g11, g01, false});
        tests.add(new Object[]{g01, g11, g11, false});
        tests.add(new Object[]{g01, g11, gNo, false});

        tests.add(new Object[]{g01, gNo, g00, false});
        tests.add(new Object[]{g01, gNo, g01, false});
        tests.add(new Object[]{g01, gNo, g11, false});
        tests.add(new Object[]{g01, gNo, gNo, false});


        tests.add(new Object[]{g11, g00, g00, true});
        tests.add(new Object[]{g11, g00, g01, false});
        tests.add(new Object[]{g11, g00, g11, true});
        tests.add(new Object[]{g11, g00, gNo, false});

        tests.add(new Object[]{g11, g01, g00, true});
        tests.add(new Object[]{g11, g01, g01, false});
        tests.add(new Object[]{g11, g01, g11, false});
        tests.add(new Object[]{g11, g01, gNo, false});

        tests.add(new Object[]{g11, g11, g00, true});
        tests.add(new Object[]{g11, g11, g01, true});
        tests.add(new Object[]{g11, g11, g11, false});
        tests.add(new Object[]{g11, g11, gNo, false});

        tests.add(new Object[]{g11, gNo, g00, true});
        tests.add(new Object[]{g11, gNo, g01, false});
        tests.add(new Object[]{g11, gNo, g11, false});
        tests.add(new Object[]{g11, gNo, gNo, false});


        tests.add(new Object[]{gNo, g00, g00, false});
        tests.add(new Object[]{gNo, g00, g01, false});
        tests.add(new Object[]{gNo, g00, g11, true});
        tests.add(new Object[]{gNo, g00, gNo, false});

        tests.add(new Object[]{gNo, g01, g00, false});
        tests.add(new Object[]{gNo, g01, g01, false});
        tests.add(new Object[]{gNo, g01, g11, false});
        tests.add(new Object[]{gNo, g01, gNo, false});

        tests.add(new Object[]{gNo, g11, g00, true});
        tests.add(new Object[]{gNo, g11, g01, false});
        tests.add(new Object[]{gNo, g11, g11, false});
        tests.add(new Object[]{gNo, g11, gNo, false});

        tests.add(new Object[]{gNo, gNo, g00, false});
        tests.add(new Object[]{gNo, gNo, g01, false});
        tests.add(new Object[]{gNo, gNo, g11, false});
        tests.add(new Object[]{gNo, gNo, gNo, false});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "isViolation")
    public void testIsViolation(final Genotype gMom, final Genotype gDad, final Genotype gChild, final boolean expected) throws Exception {
        final boolean actual = MendelianViolation.isViolation(gMom, gDad, gChild);
        Assert.assertEquals(expected, actual);
    }

    @Test(dataProvider = "isViolation")
    public void testIsViolationFromSamples(final Genotype gMom, final Genotype gDad, final Genotype gChild, final boolean expected) throws Exception {

        final Sample sMom = new Sample("mom", "fam", null, null, Sex.FEMALE);
        final Sample sDad = new Sample("dad",     sMom.getFamilyID(), null, null, Sex.MALE);
        final Sample sChild = new Sample("child", sMom.getFamilyID(), sDad.getID(), sMom.getID(), Sex.MALE);

        final Allele refAllele = Allele.create("A", true);
        final Allele noCallAllele = Allele.NO_CALL;

        final ArrayList<Genotype> genotypes = new ArrayList<>(Arrays.<Genotype>asList(gMom, gDad, gChild));
        final Map<String, Integer> sampleNameToOffset= new LinkedHashMap<>(3);
        sampleNameToOffset.put("mom", 0);
        sampleNameToOffset.put("dad", 1);
        sampleNameToOffset.put("child", 2);
        final List<String> sampleNamesInOrder=Arrays.asList("mom", "dad", "child");
        final GenotypesContext context = GenotypesContext.create(genotypes, sampleNameToOffset, sampleNamesInOrder);

        final Collection<Allele> alleles= new HashSet<>(2);
        alleles.addAll(gMom.getAlleles());
        alleles.addAll(gDad.getAlleles());
        alleles.addAll(gChild.getAlleles());
        alleles.add(refAllele);
        alleles.remove(noCallAllele);
        final VariantContext vc = new VariantContextBuilder("test", "20", 10, 10, alleles).genotypes(context).make();

        final MendelianViolation mv = new MendelianViolation(0, false, false);
        final boolean actual = mv.isViolation(sMom, sDad, sChild, vc);
        Assert.assertEquals(expected, actual);
        Assert.assertEquals(gMom.isHomRef() && gDad.isHomRef() && gChild.isHet(), mv.getParentsRefRefChildHet() > 0);

        final boolean actualFiltered = new MendelianViolation(GQ30+1, false, false).isViolation(sMom, sDad, sChild, vc);
        Assert.assertEquals(false, actualFiltered);

    }
}
