package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import com.google.common.collect.Maps;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.haplotype.EventMap;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.test.VariantContextTestUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

public class AssemblyBasedCallerGenotypingEngineUnitTest extends GATKBaseTest {

    @DataProvider(name = "getVcsAtThisLocation")
    public Object[][] getVcsAtThisLocationData() {
        final List<Object[]> tests = new ArrayList<>();

        tests.add(new Object[]{new ArrayList<>(), 1000, new ArrayList<>(), new ArrayList<>()});

        final Haplotype snpHaplotype = new Haplotype("ACTGGTCAACTGGTCAACTGGTCAACTGGTCA".getBytes());
        final List<Allele> snpAlleles = Arrays.asList(Allele.create("A", true), Allele.create("G"));
        final VariantContextBuilder snpVCBuilder = new VariantContextBuilder("a", "20", 1000, 1000, snpAlleles);
        final VariantContext snpVc = snpVCBuilder.make();
        snpHaplotype.setEventMap(new EventMap(Arrays.asList(snpVc)));

        // this one matches the snp haplotype above (to test duplicate removal)
        final Haplotype snpHaplotypeDuplicate = new Haplotype("ACTGGTCAACTGGTCAACTGGTCAACTGGACA".getBytes());
        final List<Allele> snpAlleles2 = Arrays.asList(Allele.create("A", true), Allele.create("G"));
        final VariantContextBuilder svpVC2Builder = new VariantContextBuilder("a", "20", 1000, 1000, snpAlleles2);
        final VariantContext snpVc2 = svpVC2Builder.make();
        final List<Allele> snpAlleles3 = Arrays.asList(Allele.create("T", true), Allele.create("A"));
        final VariantContextBuilder snpVC3Builder = new VariantContextBuilder("a", "20", 1020, 1020, snpAlleles3);
        final VariantContext snpVc3 = snpVC3Builder.make();
        snpHaplotypeDuplicate.setEventMap(new EventMap(Arrays.asList(snpVc2, snpVc3)));


        final Haplotype deletionHaplotype = new Haplotype("ACTGGTCAGGTCAACTGGTCA".getBytes());
        final List<Allele> deletionAlleles = Arrays.asList(Allele.create("ACTGGTCAACT", true), Allele.create("A"));
        final VariantContextBuilder deletionVCBuilder = new VariantContextBuilder("a", "20", 995, 1005, deletionAlleles);
        final VariantContext deletionVc = deletionVCBuilder.make();
        deletionHaplotype.setEventMap(new EventMap(Arrays.asList(deletionVc)));

        // matches the deletion alleles above but at a different position (to catch an edge case in duplicate removal)
        final Haplotype deletionHaplotypeFalseDuplicate = new Haplotype("ACTGGTCAGGTCAACTGGTCA".getBytes());
        final List<Allele> deletionAllelesFalseDuplicate = Arrays.asList(Allele.create("ACTGGTCAACT", true), Allele.create("A"));
        final VariantContextBuilder deletionFalseDuplicateBuilder = new VariantContextBuilder("a", "20", 998, 1008, deletionAllelesFalseDuplicate);
        final VariantContext deletionVcFalseDuplicate = deletionFalseDuplicateBuilder.make();
        deletionHaplotypeFalseDuplicate.setEventMap(new EventMap(Arrays.asList(deletionVcFalseDuplicate)));

        // doesn't overlap 1000
        final Haplotype deletionHaplotypeNoSpan = new Haplotype("CAACTGGTCAACTGGTCAACTGGTCAACTGGTCAACTGGTCA".getBytes());
        final List<Allele> deletionAllelesNoSpan = Arrays.asList(Allele.create("GTCAA", true), Allele.create("G"));
        final VariantContextBuilder deletionVcNoSpanBuilder = new VariantContextBuilder("a", "20", 990, 994, deletionAllelesNoSpan);
        final VariantContext deletionVcNoSpan = deletionVcNoSpanBuilder.make();
        deletionHaplotypeNoSpan.setEventMap(new EventMap(Arrays.asList(deletionVcNoSpan)));

        tests.add(new Object[]{Arrays.asList(snpHaplotype), 1000, new ArrayList<>(), Arrays.asList(snpVc)});
        tests.add(new Object[]{Arrays.asList(snpHaplotype, snpHaplotypeDuplicate), 1000, new ArrayList<>(), Arrays.asList(snpVc)});
        tests.add(new Object[]{Arrays.asList(deletionHaplotype), 995, new ArrayList<>(), Arrays.asList(deletionVc)});
        tests.add(new Object[]{Arrays.asList(deletionHaplotype), 1000, new ArrayList<>(), Arrays.asList(deletionVc)});
        tests.add(new Object[]{Arrays.asList(deletionHaplotype, deletionHaplotypeNoSpan), 1000, new ArrayList<>(), Arrays.asList(deletionVc)});
        tests.add(new Object[]{Arrays.asList(deletionHaplotype, deletionHaplotypeFalseDuplicate, deletionHaplotypeNoSpan), 1000, new ArrayList<>(), Arrays.asList(deletionVc, deletionVcFalseDuplicate)});

        tests.add(new Object[]{Arrays.asList(deletionHaplotype, snpHaplotype), 1000, new ArrayList<>(), Arrays.asList(deletionVc, snpVc)});

        final Haplotype sameLocDelHap1 = new Haplotype("AAAAAAAGAAA".getBytes());
        final List<Allele> sameLocDelAlleles1 = Arrays.asList(Allele.create("GTT", true), Allele.create("G"));
        final VariantContext sameLocDelVc1 = new VariantContextBuilder("a", "20", 10093568, 10093570, sameLocDelAlleles1).make();
        sameLocDelHap1.setEventMap(new EventMap(Arrays.asList(sameLocDelVc1)));

        final Haplotype sameLocDelHap2 = new Haplotype("AAAAAAAGTAAA".getBytes());
        final List<Allele> sameLocDelAlleles2 = Arrays.asList(Allele.create("GT", true), Allele.create("G"));
        final VariantContext sameLocDelVc2 = new VariantContextBuilder("a", "20", 10093568, 10093569, sameLocDelAlleles2).make();
        sameLocDelHap2.setEventMap(new EventMap(Arrays.asList(sameLocDelVc2)));

        final Haplotype sameLocDelHap3 = new Haplotype("AAAAAAAGTTTAAA".getBytes());
        final List<Allele> sameLocDelAlleles3 = Arrays.asList(Allele.create("G", true), Allele.create("GT"));
        final VariantContext sameLocDelVc3 = new VariantContextBuilder("a", "20", 10093568, 10093568, sameLocDelAlleles3).make();
        sameLocDelHap3.setEventMap(new EventMap(Arrays.asList(sameLocDelVc3)));

        tests.add(new Object[]{Arrays.asList(sameLocDelHap1, sameLocDelHap2, sameLocDelHap3), 10093568, new ArrayList<>(), Arrays.asList(sameLocDelVc1, sameLocDelVc2, sameLocDelVc3)});

        tests.add(new Object[]{new ArrayList<>(), 1000, Arrays.asList(snpVc), Arrays.asList(snpVCBuilder.source("Comp0Allele0").make())});
        tests.add(new Object[]{new ArrayList<>(), 995, Arrays.asList(deletionVc), Arrays.asList(deletionVCBuilder.source("Comp0Allele0").make())});
        tests.add(new Object[]{new ArrayList<>(), 1000, Arrays.asList(deletionVc), Arrays.asList(deletionVCBuilder.source("Comp0Allele0").make())});
        tests.add(new Object[]{new ArrayList<>(), 1000, Arrays.asList(deletionVc, snpVc),
                Arrays.asList(deletionVCBuilder.source("Comp0Allele0").make(), snpVCBuilder.source("Comp1Allele0").make())});
        tests.add(new Object[]{new ArrayList<>(), 1000, Arrays.asList(deletionVc, deletionVcNoSpan), Arrays.asList(deletionVCBuilder.source("Comp0Allele0").make())});
        tests.add(new Object[]{new ArrayList<>(), 1000, Arrays.asList(deletionVc, deletionVcFalseDuplicate, deletionVcNoSpan),
                Arrays.asList(deletionVCBuilder.source("Comp0Allele0").make(), deletionFalseDuplicateBuilder.source("Comp1Allele0").make())});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "getVcsAtThisLocation")
    public void testGetVCsAtThisLocation(final List<Haplotype> haplotypes,
                                         final int loc,
                                         final List<VariantContext> activeAllelesToGenotype,
                                         final List<VariantContext> expectedVcsAtThisLocation) {

        final List<VariantContext> vcsAtThisPosition = AssemblyBasedCallerGenotypingEngine.getVCsAtThisLocation(haplotypes, loc, activeAllelesToGenotype);
        Assert.assertEquals(vcsAtThisPosition.size(), expectedVcsAtThisLocation.size());
        for (int i = 0; i < expectedVcsAtThisLocation.size(); i++) {
            VariantContextTestUtils.assertVariantContextsAreEqual(vcsAtThisPosition.get(i), expectedVcsAtThisLocation.get(i), new ArrayList<>());
            Assert.assertEquals(vcsAtThisPosition.get(i).getSource(), expectedVcsAtThisLocation.get(i).getSource());
        }
    }


    @DataProvider(name = "getEventMapper")
    public Object[][] getEventMapperData() {

        final Haplotype refHaplotype = new Haplotype("ACTGGTCAACTAGTCAACTGGTCAACTGGTCA".getBytes());
        refHaplotype.setEventMap(new EventMap(new HashSet<>()));

        final Haplotype snpHaplotype = new Haplotype("ACTGGTCAACTGGTCAACTGGTCAACTGGTCA".getBytes());
        final Allele refAllele = Allele.create("A", true);
        final List<Allele> snpAlleles = Arrays.asList(refAllele, Allele.create("G"));
        final VariantContextBuilder snpVCBuilder = new VariantContextBuilder("a", "20", 1000, 1000, snpAlleles);
        final VariantContext snpVc = snpVCBuilder.make();
        snpHaplotype.setEventMap(new EventMap(Arrays.asList(snpVc)));

        final Haplotype snpHaplotypeNotPresentInEventsAtThisLoc = new Haplotype("ACTGGTCAACTTGTCAACTGGTCAACTGGTCA".getBytes());
        final List<Allele> snpAllelesNotPresentInEventsAtThisLoc = Arrays.asList(refAllele, Allele.create("T"));
        final VariantContextBuilder snpNotPresentInEventsAtThisLocVCBuilder = new VariantContextBuilder("a", "20", 1000, 1000, snpAllelesNotPresentInEventsAtThisLoc);
        final VariantContext snpVcNotPresentInEventsAtThisLoc = snpNotPresentInEventsAtThisLocVCBuilder.make();
        snpHaplotypeNotPresentInEventsAtThisLoc.setEventMap(new EventMap(Arrays.asList(snpVcNotPresentInEventsAtThisLoc)));

        final Haplotype deletionHaplotype = new Haplotype("ACTGGTCAGGTCAACTGGTCA".getBytes());
        final List<Allele> deletionAlleles = Arrays.asList(Allele.create("ACTGGTCAACT", true), Allele.create("A"));
        final VariantContextBuilder deletionVCBuilder = new VariantContextBuilder("a", "20", 995, 1005, deletionAlleles);
        final VariantContext deletionVc = deletionVCBuilder.make();
        deletionHaplotype.setEventMap(new EventMap(Arrays.asList(deletionVc)));

        final VariantContext spandDelVc = new VariantContextBuilder("a", "20", 1000, 1000, Arrays.asList(refAllele, Allele.SPAN_DEL)).make();

        final VariantContext mergedSnpAndDelVC = new VariantContextBuilder("a", "20", 1000, 1000,
                Arrays.asList(refAllele,
                        Allele.SPAN_DEL,
                        Allele.create("G"))).make();



        final List<Object[]> tests = new ArrayList<>();
        tests.add(new Object[]{
                Arrays.asList(snpVc),
                snpVc,
                snpVc.getStart(),
                Arrays.asList(snpHaplotype, refHaplotype),
                Maps.asMap(new HashSet<>(snpAlleles),
                (key) -> {
                    if (snpAlleles.get(1).equals(key)) return Arrays.asList(snpHaplotype);
                    return Arrays.asList(refHaplotype);
                })
        });
        tests.add(new Object[]{
                Arrays.asList(spandDelVc, snpVc),
                mergedSnpAndDelVC,
                mergedSnpAndDelVC.getStart(),
                Arrays.asList(snpHaplotype, refHaplotype, deletionHaplotype),
                Maps.asMap(new HashSet<>(mergedSnpAndDelVC.getAlleles()),
                        (key) -> {
                            if (snpAlleles.get(1).equals(key)) return Arrays.asList(snpHaplotype);
                            if (Allele.SPAN_DEL.equals(key)) return Arrays.asList(deletionHaplotype);
                            return Arrays.asList(refHaplotype);
                        })
        });
        // includes a SNP haplotype not present in events at this loc (which might happen in GGA mode)
        tests.add(new Object[]{
                Arrays.asList(snpVc),
                snpVc,
                snpVc.getStart(),
                Arrays.asList(snpHaplotype, refHaplotype, snpHaplotypeNotPresentInEventsAtThisLoc),
                Maps.asMap(new HashSet<>(snpVc.getAlleles()),
                        (key) -> {
                            if (snpAlleles.get(1).equals(key)) return Arrays.asList(snpHaplotype);
                            return Arrays.asList(refHaplotype);
                        })
        });


        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "getEventMapper")
    public void testGetEventMapper(final List<VariantContext> eventsAtThisLoc,
                                   final VariantContext mergedVc,
                                   final int loc,
                                   final List<Haplotype> haplotypes,
                                   final Map<Allele, List<Haplotype>> expectedEventMap) {
        final Map<Allele, List<Haplotype>> actualEventMap = AssemblyBasedCallerGenotypingEngine.createAlleleMapper(eventsAtThisLoc, mergedVc, loc, haplotypes, new ArrayList<>());
        Assert.assertEquals(actualEventMap.size(), expectedEventMap.size());
        for (final Allele key : actualEventMap.keySet()) {
            Assert.assertTrue(expectedEventMap.containsKey(key), "Got unexpected allele " + key + " with values " + actualEventMap.get(key));
            Assert.assertEquals(actualEventMap.get(key), expectedEventMap.get(key), "Lists don't match for key " + key);
        }

        for (final Allele key : expectedEventMap.keySet()) {
            Assert.assertTrue(actualEventMap.containsKey(key), "Didn't get back allele " + key);
        }

    }
}