package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AlleleFiltering;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AlleleFilteringHC;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerGenotypingEngine;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.genotyper.*;
import org.broadinstitute.hellbender.utils.haplotype.EventMap;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.*;

public class AlleleFilteringUnitTest {

    @Test
    public void testNoNeedToFilter(){
        List<Haplotype> haplotypeList = new ArrayList<>();
        final byte[] fullReferenceWithPadding = "CATGCATG".getBytes();
        Haplotype haplotype = new Haplotype(fullReferenceWithPadding, true, 0, TextCigarCodec.decode("8M"));
        haplotype.setGenomeLocation(new SimpleInterval("chr", 100, 108));
        haplotype.setEventMap(new EventMap(haplotype, fullReferenceWithPadding,
                new SimpleInterval("chr", (int)haplotype.getStartPosition(),(int)haplotype.getStopPosition()), "test", 0));


        haplotypeList.add(haplotype);
        haplotype = new Haplotype("CAGGCATG".getBytes(), false, 0, TextCigarCodec.decode("8M"));
        haplotype.setGenomeLocation(new SimpleInterval("chr", 100, 108));

        haplotypeList.add(haplotype);
        haplotype.setEventMap(new EventMap(haplotype, fullReferenceWithPadding,
                new SimpleInterval("chr", (int)haplotype.getStartPosition(), (int)haplotype.getStopPosition()), "test", 0));

        AlleleList<Haplotype> haplotypes = new IndexedAlleleList<>(haplotypeList);
        SampleList samples = new IndexedSampleList(Arrays.asList("sm1"));

        List<GATKRead> readList = new ArrayList<>(30);
        Map<String, List<GATKRead>>ebs = new HashMap<>();
        ebs.put("sm1", readList);

        for (int i = 0 ; i < 30; i++) {
            readList.add(ArtificialReadUtils.createArtificialRead("20M"));
        }

        double[][] values = {{0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3},
                                {3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0}};

        AlleleLikelihoods<GATKRead, Haplotype> lks = new AlleleLikelihoods<>(samples, haplotypes, ebs);
        LikelihoodMatrix<GATKRead, Haplotype> lkm = lks.sampleMatrix(0);
        for (int i = 0; i < lks.numberOfAlleles(); i++){
            for (int j = 0 ; j < lkm.evidenceCount(); j++) {
                lkm.set(i,j,values[i][j]);
            }
        }

        HaplotypeCallerArgumentCollection hcArgs = new HaplotypeCallerArgumentCollection();
        HaplotypeCallerGenotypingEngine genotypingEngine = new HaplotypeCallerGenotypingEngine(hcArgs, samples, ! hcArgs.doNotRunPhysicalPhasing, false);


        AlleleFiltering alleleFiltering = new AlleleFilteringHC(hcArgs, null, genotypingEngine);
        AlleleLikelihoods<GATKRead, Haplotype> filtered_lks = alleleFiltering.filterAlleles(lks, 0, new HashSet<>());
        Assert.assertEquals(filtered_lks.alleles(), lks.alleles());
    }

    @Test
    public void testFilterCloseMis(){

        // create haplotypes
        List<Haplotype> haplotypeList = new ArrayList<>();
        final byte[] fullReferenceWithPadding = "CATGCATG".getBytes();
        Haplotype haplotype = new Haplotype(fullReferenceWithPadding, true, 0, TextCigarCodec.decode("8M"));
        haplotype.setGenomeLocation(new SimpleInterval("chr", 100, 108));
        haplotype.setEventMap(new EventMap(haplotype, fullReferenceWithPadding,
                new SimpleInterval("chr", (int)haplotype.getStartPosition(),(int)haplotype.getStopPosition()), "test", 0));
        haplotypeList.add(haplotype);

        haplotype = new Haplotype("CAGGCATG".getBytes(), false, 0, TextCigarCodec.decode("8M"));
        haplotype.setGenomeLocation(new SimpleInterval("chr", 100, 108));

        haplotype.setEventMap(new EventMap(haplotype, fullReferenceWithPadding,
                new SimpleInterval("chr", (int)haplotype.getStartPosition(), (int)haplotype.getStopPosition()), "test", 0));
        haplotypeList.add(haplotype);

        haplotype = new Haplotype("CAGTCATG".getBytes(), false, 0, TextCigarCodec.decode("8M"));
        haplotype.setGenomeLocation(new SimpleInterval("chr", 100, 108));

        haplotype.setEventMap(new EventMap(haplotype, fullReferenceWithPadding,
                new SimpleInterval("chr", (int)haplotype.getStartPosition(), (int)haplotype.getStopPosition()), "test", 0));
        haplotypeList.add(haplotype);

        AlleleList<Haplotype> haplotypes = new IndexedAlleleList<>(haplotypeList);
        SampleList samples = new IndexedSampleList(Arrays.asList("sm1"));



        List<GATKRead> readList = new ArrayList<>(30);
        Map<String, List<GATKRead>>ebs = new HashMap<>();
        ebs.put("sm1", readList);

        for (int i = 0 ; i < 30; i++) {
            readList.add(ArtificialReadUtils.createArtificialRead("20M"));
        }

        double[][] values = {{0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3},
                {3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0},
                {2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
        };

        AlleleLikelihoods<GATKRead, Haplotype> lks = new AlleleLikelihoods<>(samples, haplotypes, ebs);
        LikelihoodMatrix<GATKRead, Haplotype> lkm = lks.sampleMatrix(0);
        for (int i = 0; i < lks.numberOfAlleles(); i++){
            for (int j = 0 ; j < lkm.evidenceCount(); j++) {
                lkm.set(i,j,values[i][j]);
            }
        }

        HaplotypeCallerArgumentCollection hcArgs = new HaplotypeCallerArgumentCollection();
        HaplotypeCallerGenotypingEngine genotypingEngine = new HaplotypeCallerGenotypingEngine(hcArgs, samples, ! hcArgs.doNotRunPhysicalPhasing, false);


        AlleleFiltering alleleFiltering = new AlleleFilteringHC(hcArgs, null, genotypingEngine);
        AlleleLikelihoods<GATKRead, Haplotype> filtered_lks = alleleFiltering.filterAlleles(lks, 0, new HashSet<>());
        Assert.assertEquals(filtered_lks.alleles(), haplotypeList.subList(0,2));
    }

    @Test
    public void testFilterDistantHindel(){

        // create haplotypes
        List<Haplotype> haplotypeList = new ArrayList<>();
        final byte[] fullReferenceWithPadding = "CATGCATG".getBytes();
        Haplotype haplotype = new Haplotype(fullReferenceWithPadding, true, 0, TextCigarCodec.decode("8M"));
        haplotype.setGenomeLocation(new SimpleInterval("chr", 100, 108));
        haplotype.setEventMap(new EventMap(haplotype, fullReferenceWithPadding,
                new SimpleInterval("chr", (int)haplotype.getStartPosition(),(int)haplotype.getStopPosition()), "test", 0));
        haplotypeList.add(haplotype);

        haplotype = new Haplotype("CAGGCATG".getBytes(), false, 0, TextCigarCodec.decode("8M"));
        haplotype.setGenomeLocation(new SimpleInterval("chr", 100, 108));

        haplotype.setEventMap(new EventMap(haplotype, fullReferenceWithPadding,
                new SimpleInterval("chr", (int)haplotype.getStartPosition(), (int)haplotype.getStopPosition()), "test", 0));
        haplotypeList.add(haplotype);

        haplotype = new Haplotype("CAGGCATTG".getBytes(), false, 0, TextCigarCodec.decode("7M1I1M"));
        haplotype.setGenomeLocation(new SimpleInterval("chr", 100, 109));

        haplotype.setEventMap(new EventMap(haplotype, fullReferenceWithPadding,
                new SimpleInterval("chr", (int)haplotype.getStartPosition(), (int)haplotype.getStopPosition()), "test", 0));
        haplotypeList.add(haplotype);

        AlleleList<Haplotype> haplotypes = new IndexedAlleleList<>(haplotypeList);
        SampleList samples = new IndexedSampleList(Arrays.asList("sm1"));



        List<GATKRead> readList = new ArrayList<>(30);
        Map<String, List<GATKRead>>ebs = new HashMap<>();
        ebs.put("sm1", readList);

        for (int i = 0 ; i < 30; i++) {
            readList.add(ArtificialReadUtils.createArtificialRead("20M"));
        }

        double[][] values = {{0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3},
                {3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0},
                {2, 0, 1, 0, 2, 0, 2, 0, 2, 0, 1, 0, 4, 0, 2, 0, 2, 0, 2, 0, 2, 0, 4, 0, 2, 0, 2, 0, 2, 0}
        };

        AlleleLikelihoods<GATKRead, Haplotype> lks = new AlleleLikelihoods<>(samples, haplotypes, ebs);
        LikelihoodMatrix<GATKRead, Haplotype> lkm = lks.sampleMatrix(0);
        for (int i = 0; i < lks.numberOfAlleles(); i++){
            for (int j = 0 ; j < lkm.evidenceCount(); j++) {
                lkm.set(i,j,values[i][j]);
            }
        }

        HaplotypeCallerArgumentCollection hcArgs = new HaplotypeCallerArgumentCollection();
        HaplotypeCallerGenotypingEngine genotypingEngine = new HaplotypeCallerGenotypingEngine(hcArgs, samples, ! hcArgs.doNotRunPhysicalPhasing, false);


        AlleleFiltering alleleFiltering = new AlleleFilteringHC(hcArgs, null, genotypingEngine);
        AlleleLikelihoods<GATKRead, Haplotype> filtered_lks = alleleFiltering.filterAlleles(lks, 0, new HashSet<>());
        Assert.assertEquals(filtered_lks.alleles(), haplotypeList.subList(0,2));
    }

    @Test
    public void testNotFilterDistantM(){

        // create haplotypes
        List<Haplotype> haplotypeList = new ArrayList<>();
        final byte[] fullReferenceWithPadding = "CATGCATG".getBytes();
        Haplotype haplotype = new Haplotype(fullReferenceWithPadding, true, 0, TextCigarCodec.decode("8M"));
        haplotype.setGenomeLocation(new SimpleInterval("chr", 100, 108));
        haplotype.setEventMap(new EventMap(haplotype, fullReferenceWithPadding,
                new SimpleInterval("chr", (int)haplotype.getStartPosition(),(int)haplotype.getStopPosition()), "test", 0));
        haplotypeList.add(haplotype);

        haplotype = new Haplotype("CAGGCATG".getBytes(), false, 0, TextCigarCodec.decode("8M"));
        haplotype.setGenomeLocation(new SimpleInterval("chr", 100, 108));

        haplotype.setEventMap(new EventMap(haplotype, fullReferenceWithPadding,
                new SimpleInterval("chr", (int)haplotype.getStartPosition(), (int)haplotype.getStopPosition()), "test", 0));
        haplotypeList.add(haplotype);

        haplotype = new Haplotype("CATGCATC".getBytes(), false, 0, TextCigarCodec.decode("8M"));
        haplotype.setGenomeLocation(new SimpleInterval("chr", 100, 108));

        haplotype.setEventMap(new EventMap(haplotype, fullReferenceWithPadding,
                new SimpleInterval("chr", (int)haplotype.getStartPosition(), (int)haplotype.getStopPosition()), "test", 0));
        haplotypeList.add(haplotype);

        AlleleList<Haplotype> haplotypes = new IndexedAlleleList<>(haplotypeList);
        SampleList samples = new IndexedSampleList(Arrays.asList("sm1"));



        List<GATKRead> readList = new ArrayList<>(30);
        Map<String, List<GATKRead>>ebs = new HashMap<>();
        ebs.put("sm1", readList);

        for (int i = 0 ; i < 30; i++) {
            readList.add(ArtificialReadUtils.createArtificialRead("20M"));
        }

        double[][] values = {{0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3},
                {3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0},
                {3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0}
        };

        AlleleLikelihoods<GATKRead, Haplotype> lks = new AlleleLikelihoods<>(samples, haplotypes, ebs);
        LikelihoodMatrix<GATKRead, Haplotype> lkm = lks.sampleMatrix(0);
        for (int i = 0; i < lks.numberOfAlleles(); i++){
            for (int j = 0 ; j < lkm.evidenceCount(); j++) {
                lkm.set(i,j,values[i][j]);
            }
        }

        HaplotypeCallerArgumentCollection hcArgs = new HaplotypeCallerArgumentCollection();
        HaplotypeCallerGenotypingEngine genotypingEngine = new HaplotypeCallerGenotypingEngine(hcArgs, samples, ! hcArgs.doNotRunPhysicalPhasing, false);


        AlleleFiltering alleleFiltering = new AlleleFilteringHC(hcArgs, null, genotypingEngine);
        AlleleLikelihoods<GATKRead, Haplotype> filtered_lks = alleleFiltering.filterAlleles(lks, 100, new HashSet<>());
        Assert.assertEquals(filtered_lks.alleles(), haplotypeList.subList(0,1));
    }

    @Test
    public void testNotFilterLoneWeakAllele(){
        //if hcArgs.filterLoneAlleles is false - weak lone alleles should be kept

        // create haplotypes
        List<Haplotype> haplotypeList = new ArrayList<>();
        final byte[] fullReferenceWithPadding = "CATGCATG".getBytes();
        Haplotype haplotype = new Haplotype(fullReferenceWithPadding, true, 0, TextCigarCodec.decode("8M"));
        haplotype.setGenomeLocation(new SimpleInterval("chr", 100, 108));
        haplotype.setEventMap(new EventMap(haplotype, fullReferenceWithPadding,
                new SimpleInterval("chr", (int)haplotype.getStartPosition(),(int)haplotype.getStopPosition()), "test", 0));
        haplotypeList.add(haplotype);

        haplotype = new Haplotype("CAGGCATG".getBytes(), false, 0, TextCigarCodec.decode("8M"));
        haplotype.setGenomeLocation(new SimpleInterval("chr", 100, 108));

        haplotype.setEventMap(new EventMap(haplotype, fullReferenceWithPadding,
                new SimpleInterval("chr", (int)haplotype.getStartPosition(), (int)haplotype.getStopPosition()), "test", 0));
        haplotypeList.add(haplotype);

        AlleleList<Haplotype> haplotypes = new IndexedAlleleList<>(haplotypeList);
        SampleList samples = new IndexedSampleList(Arrays.asList("sm1"));

        List<GATKRead> readList = new ArrayList<>(30);
        Map<String, List<GATKRead>>ebs = new HashMap<>();
        ebs.put("sm1", readList);

        for (int i = 0 ; i < 30; i++) {
            readList.add(ArtificialReadUtils.createArtificialRead("20M"));
        }

        double[][] values = {{3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0},
                {3, 0, 4, 0, 3, 0, 4, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0}
        };

        AlleleLikelihoods<GATKRead, Haplotype> lks = new AlleleLikelihoods<>(samples, haplotypes, ebs);
        LikelihoodMatrix<GATKRead, Haplotype> lkm = lks.sampleMatrix(0);
        for (int i = 0; i < lks.numberOfAlleles(); i++){
            for (int j = 0 ; j < lkm.evidenceCount(); j++) {
                lkm.set(i,j,values[i][j]);
            }
        }

        HaplotypeCallerArgumentCollection hcArgs = new HaplotypeCallerArgumentCollection();
        HaplotypeCallerGenotypingEngine genotypingEngine = new HaplotypeCallerGenotypingEngine(hcArgs, samples, ! hcArgs.doNotRunPhysicalPhasing, false);


        AlleleFiltering alleleFiltering = new AlleleFilteringHC(hcArgs, null, genotypingEngine);
        AlleleLikelihoods<GATKRead, Haplotype> filtered_lks = alleleFiltering.filterAlleles(lks, 100, new HashSet<>());
        //hcArgs.filterLoneAlleles is false, so we keep the weak lone allele
        Assert.assertEquals(filtered_lks.alleles(), haplotypeList);

        hcArgs.filterLoneAlleles = true;
        alleleFiltering = new AlleleFilteringHC(hcArgs, null, genotypingEngine);
        filtered_lks = alleleFiltering.filterAlleles(lks, 100, new HashSet<>());

        //hcArgs.filterLoneAlleles is true, so we keep the remove the weak lone allele
        Assert.assertEquals(filtered_lks.alleles(), haplotypeList.subList(0,1));

    }


}
