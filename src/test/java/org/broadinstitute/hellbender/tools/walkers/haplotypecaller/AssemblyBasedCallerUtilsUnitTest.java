package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import com.google.common.collect.Maps;
import htsjdk.samtools.*;
import htsjdk.variant.variantcontext.*;

import htsjdk.samtools.util.Locatable;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.engine.AssemblyRegion;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.genotyper.*;
import org.broadinstitute.hellbender.utils.haplotype.EventMap;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.apache.commons.lang3.StringUtils;

import java.util.*;
import java.util.stream.Collectors;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAligner;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

public class AssemblyBasedCallerUtilsUnitTest extends GATKBaseTest {
    final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader(1, 1, 100000000);
    final SAMLineParser parser = new SAMLineParser(header);

    // In finalizeRegion(), the base qualities of overlapped read clips pairs are adjusted.
    // Most of read clips are clipped/copy from original reads, and the base qualities of original reads are not affected.
    // However, GATK 4.0.0.x has a bug : in some cases, the clipping procedure directly returns the original reads.
    // So the base qualities of original reads are changed.
    // This test is added to make sure the unexpected behavior is fixed.
    @Test
    public void testfinalizeRegion() {
        SAMFileHeader header;
        header = ArtificialReadUtils.createArtificialSamHeader(1, 1, 100000000);
        SAMReadGroupRecord rg = new SAMReadGroupRecord("tumor");
        rg.setSample("tumor");
        header.addReadGroup(rg);

        Assert.assertEquals(header.getSequenceIndex("1"), 0);
        final AssemblyRegion activeRegion = new AssemblyRegion(new SimpleInterval("1",42596728,42598843), 100, header);
        Assert.assertTrue(activeRegion.isActive());
        Assert.assertFalse(activeRegion.isFinalized());

        SAMLineParser parser = new SAMLineParser(header);
        List<GATKRead> reads = new LinkedList<GATKRead>();
        SAMRecord orgRead0 = parser.parseLine("HWI-ST807:461:C2P0JACXX:4:2204:18080:5857\t83\t1\t42596803\t39\t1S95M5S\t=\t42596891\t-7\tGAATCATCATCAAATGGAATCTAATGGAATCATTGAACAGAATTGAATGGAATCGTCATCGAATGAATTGAATGCAATCATCGAATGGTCTCGAATAGAAT\tDAAAEDCFCCGEEDDBEDDDGCCDEDECDDFDCEECCFEECDCEDBCDBDBCC>DCECC>DBCDDBCBDDBCDDEBCCECC>DBCDBDBGC?FCCBDB>>?\tRG:Z:tumor");
        SAMRecord orgRead1 = parser.parseLine("HWI-ST807:461:C2P0JACXX:4:2204:18080:5857\t163\t1\t42596891\t39\t101M\t=\t42596803\t7\tCTCGAATGGAATCATTTTCTACTGGAAAGGAATGGAATCATCGCATAGAATCGAATGGAATTAACATGGAATGGAATCGAATGTAATCATCATCAAATGGA\t>@>:ABCDECCCEDCBBBDDBDDEBCCBEBBCBEBCBCDDCD>DECBGCDCF>CCCFCDDCBABDEDFCDCDFFDDDG?DDEGDDFDHFEGDDGECB@BAA\tRG:Z:tumor");

        // use a copy of the original reads so the original reads are not touched
        reads.add(new SAMRecordToGATKReadAdapter(orgRead0.deepCopy()));
        reads.add(new SAMRecordToGATKReadAdapter(orgRead1.deepCopy()));

        // add reads into active region and call finalizeRegion
        activeRegion.addAll(reads);
        SampleList sampleList = SampleList.singletonSampleList("tumor");
        Byte minbq = 9;
        AssemblyBasedCallerUtils.finalizeRegion(activeRegion, false, false, minbq, header, sampleList, false);

        // make sure reads are not changed due to finalizeRegion()
        Assert.assertTrue(reads.get(0).convertToSAMRecord(header).equals(orgRead0));
        Assert.assertTrue(reads.get(1).convertToSAMRecord(header).equals(orgRead1));
    }

    // ------------------------------------------------------------------------
    //
    //  Test annotation of reads for bamout
    //
    // ------------------------------------------------------------------------
    @DataProvider(name = "testAnnotateReadLikelihoodsWithRegionsDataProvider")
    public Object[][] testAnnotateReadLikelihoodsWithRegionsDataProvider() {

        final String hap1String = "ACGTGGCGTTGCACTTCAGATCGATCGGATCGATCGGCTAGTCGTCGCACTTCGCTAGGCTAG";
        final String contig = "1";
        final int start = 13763;
        final SimpleInterval loc = new SimpleInterval(contig, start, start + hap1String.length());
        final SimpleInterval loc2 = new SimpleInterval(contig, loc.getStart() + 4, loc.getEnd() - 4);
        List<Integer> snpIndices = Arrays.asList(null, 12, 23);
        List<Character> snpReplacements = Arrays.asList(null, 'G', 'C');
        List<Integer> delIndices = Arrays.asList(null, 27);
        List<Integer> delLenghts = Arrays.asList(null, 4);
        List<Integer> insIndices = Arrays.asList(null, 34);
        List<String> insStrings = Arrays.asList(null, "GGCTGGATCGAG");

        List<String> haplotypeStrings = new ArrayList<String>();

        for (int iSnp = 0; iSnp < snpIndices.size(); iSnp++) {
            //create set of haplotype sequences with different combinations of snps, deletions, indels
            Integer snpIndex = snpIndices.get(iSnp);
            Character snpReplacement = snpReplacements.get(iSnp);
            for (int iDel = 0; iDel < delIndices.size(); iDel++) {
                Integer delIndex = delIndices.get(iDel);
                Integer delLength = delLenghts.get(iDel);
                for (int iIns = 0; iIns < insIndices.size(); iIns++) {
                    Integer insIndex = insIndices.get(iIns);
                    String insString = insStrings.get(iIns);
                    String haplotypeString = hap1String;
                    if (snpIndex != null) {
                        haplotypeString = applySNP(haplotypeString, snpIndex, snpReplacement);
                    }
                    if (insIndex != null) {
                        haplotypeString = applyInsertion(haplotypeString, insIndex, insString);
                    }
                    if (delIndex != null) {
                        haplotypeString = applyDeletion(haplotypeString, delIndex, delLength);
                    }
                    haplotypeStrings.add(haplotypeString);
                }
            }
        }
        final List<Haplotype> haplotypesList = haplotypeStrings.stream().map(h -> new Haplotype(h.getBytes(), loc)).collect(Collectors.toList());

        int qNameIndex = 0;
        List<GATKRead> reads = new ArrayList<>();
        Map<GATKRead, Haplotype> readHaplotypeMap = new HashMap<>();
        for (Haplotype haplotype : haplotypesList) {
            //create a bunch of reads for different haplotypes
            final String hapString = haplotype.getBaseString();
            for (int readStart = 0; readStart < hapString.length() - 6; readStart += 3) {
                for (int readEnd = Math.min(hapString.length(), readStart + 20); readEnd < Math.min(hapString.length(), readStart + 40); readEnd += 4) {
                    final String readString = hapString.substring(readStart, readEnd);
                    final String cigar = readString.length() + "M"; //actual cigar is unimportant for this test
                    final String qname = "r" + qNameIndex;
                    qNameIndex++;
                    final SAMRecordToGATKReadAdapter read = buildRead(readString,cigar,start+readStart,qname,contig);
                    reads.add(read);
                    readHaplotypeMap.put(read, haplotype);
                }
            }
        }
        Map<String, List<GATKRead>> sampleReadMap = new HashMap<>();
        sampleReadMap.put("sample1", reads);
        final AlleleList<Haplotype> haplotypes = new IndexedAlleleList<>(haplotypesList);
        final SampleList samples = new IndexedSampleList("sample1");

        final ReadLikelihoods<Haplotype> readLikelihoods = new ReadLikelihoods<>(samples, haplotypes, sampleReadMap);
        LikelihoodMatrix<Haplotype> sampleMatrix = readLikelihoods.sampleMatrix(0);
        for (GATKRead read : reads) {
            //set likelihoods, -1.0 for haplotype read assigned to, -8.0 for all other haplotypes
            final int readIndex = sampleMatrix.indexOfRead(read);
            for (Haplotype haplotype : haplotypesList) {
                final int haplotypeIndex = sampleMatrix.indexOfAllele(haplotype);
                if (readHaplotypeMap.get(read) == haplotype) {
                    sampleMatrix.set(haplotypeIndex, readIndex, -1.0);
                } else {
                    sampleMatrix.set(haplotypeIndex, readIndex, -8.0);
                }
            }
        }

        return new Object[][]{{readLikelihoods, loc, loc},
                {readLikelihoods, loc, loc2}
        };

    }


    private String applySNP(final String initialSeq, final int iReplace, final char replaceWith) {
        if (initialSeq.length() <= iReplace) {
            return initialSeq;
        }
        String retSeq = initialSeq.substring(0, iReplace) + replaceWith;
        if (initialSeq.length() > iReplace + 1) {
            retSeq += initialSeq.substring(iReplace + 1);
        }
        return retSeq;
    }

    private String applyDeletion(final String initialSeq, final int iDelete, final int lengthDelete) {
        if (initialSeq.length() <= iDelete) {
            return initialSeq;
        }
        String retSeq = initialSeq.substring(0, iDelete);
        if (initialSeq.length() > iDelete + lengthDelete) {
            retSeq += initialSeq.substring(iDelete + lengthDelete);
        }
        return retSeq;
    }

    private String applyInsertion(final String initialSeq, final int iInsert, final String insertionSeq) {
        if (initialSeq.length() <= iInsert) {
            return initialSeq;
        }
        String retSeq = initialSeq.substring(0, iInsert) + insertionSeq;
        if (initialSeq.length() > iInsert + 1) {
            retSeq += initialSeq.substring(iInsert + 1);
        }
        return retSeq;
    }

    private SAMRecordToGATKReadAdapter buildRead(final String seq, final String cigar, final int pos, final String qName,final String rName) {
        final String baseQuals = StringUtils.repeat("<", seq.length());
        final int mapq = 39;
        final String rnext = "=";
        final int pnext = pos + 100;
        final int tlen = 200;
        final int flag = 83;
        final String samLine = qName + "\t" + flag + "\t" + rName + "\t" + pos + "\t" + mapq + "\t" + cigar + "\t" + rnext + "\t" +
                pnext + "\t" + tlen + "\t" + seq + "\t" + baseQuals;
        final SAMRecord samRead = parser.parseLine(samLine);
        final SAMRecordToGATKReadAdapter read = new SAMRecordToGATKReadAdapter(samRead);

        return read;
    }

    @Test(dataProvider = "testAnnotateReadLikelihoodsWithRegionsDataProvider")
    public void testAnnotateReadLikelihoodsWithRegions(ReadLikelihoods<Haplotype> readLikelihoods, final Locatable loc, final Locatable callableLoc) {
        AssemblyBasedCallerUtils.annotateReadLikelihoodsWithRegions(readLikelihoods, callableLoc);
        for (GATKRead read : readLikelihoods.sampleReads(0)) {
            Assert.assertEquals(read.getAttributeAsString(AssemblyBasedCallerUtils.ALIGNMENT_REGION_TAG), loc.toString());
            Assert.assertEquals(read.getAttributeAsString(AssemblyBasedCallerUtils.CALLABLE_REGION_TAG), callableLoc.toString());
        }
    }

    @DataProvider(name = "testAnnotateReadLikelihoodsWithSupportedAllelesDataProvider")
    public Object[][] testAnnotateReadLikelihoodsWithSupportedAllelesDataProvider() {
        final String refString = "ACGTGGCGTTGCACTTCAGATCGATCGGATCGATCGGCTAGTCGTCGCACTTCGCTAGGCTAG";
        final String contig = "1";
        final int start = 13763;
        final SimpleInterval loc = new SimpleInterval(contig, start, start + refString.length());
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader(1, 1, 100000000);
        SAMLineParser parser = new SAMLineParser(header);

        final Haplotype refHaplotype = new Haplotype(refString.getBytes(), loc);
        final Cigar refCigar = new Cigar(Arrays.asList(new CigarElement(refString.length(), CigarOperator.M)));
        refHaplotype.setCigar(refCigar);

        //add three variants (2 snps and an indel)
        List<VariantContext> vcs = new ArrayList<>();
        vcs.add(new VariantContextBuilder("source", contig, start + 3, start + 3, Arrays.asList(
                Allele.create("T".getBytes(), true), Allele.create("G".getBytes(), false))).make());
        vcs.add(new VariantContextBuilder("source", contig, start + 38, start + 38, Arrays.asList(
                Allele.create("T".getBytes(), true), Allele.create("G".getBytes(), false), Allele.create("C".getBytes(), false))).make());
        vcs.add(new VariantContextBuilder("source", contig, start + 20, start + 23, Arrays.asList(
                Allele.create("TCGA", true), Allele.create("T", false))).make());

        Map<VariantContext, List<GATKRead>> vcReadsMap = new HashMap<>();
        Map<GATKRead, Integer> supportedAlleleMap = new HashMap<>();
        for (VariantContext vc : vcs) {
            vcReadsMap.put(vc, new ArrayList<>());
        }

        int qNameIndex = 0;
        //three true haplotypes.  All alt for first snp and indel, cycle through 3 alleles for second snp
        for (int readStart = loc.getStart(); readStart < loc.getEnd() - 6; readStart += 4) {
            for (int readEnd = Math.min(readStart + 20, loc.getEnd()); readEnd < Math.min(readStart + 50, loc.getEnd()); readEnd += 3) {
                final SimpleInterval readLoc = new SimpleInterval(contig, readStart, readEnd);
                final Haplotype refRead = refHaplotype.trim(readLoc);
                for (int phase = 0; phase < 3; phase++) {
                    Haplotype readHaplotype = applyVariant(vcs.get(0), refRead, 0);
                    if (phase > 0) {
                        readHaplotype = applyVariant(vcs.get(1), readHaplotype, phase - 1);
                    }
                    readHaplotype = applyVariant(vcs.get(2), readHaplotype, 0);
                    final String readString = readHaplotype.getBaseString();
                    final String cigar = readString.length() + "M"; //actual cigar is unimportant for this test
                    final String qname = "r" + qNameIndex;
                    qNameIndex++;
                    final SAMRecordToGATKReadAdapter read = buildRead(readString,cigar,readStart,qname,contig);
                    supportedAlleleMap.put(read, phase);
                    for (final VariantContext vc : vcs) {
                        if (read.getStart() <= vc.getStart() && read.getEnd() >= vc.getEnd()) {
                            vcReadsMap.get(vc).add(read);
                        }
                    }

                }
            }
        }
        List<ReadLikelihoods<Allele>> readLikelihoodsList = new ArrayList<>();
        List<List<String>> readAttributeListList = new ArrayList<>();
        for (VariantContext vc : vcs) {
            Map<String, List<GATKRead>> sampleReadMap = new HashMap<>();
            sampleReadMap.put("sample1", vcReadsMap.get(vc));
            final SampleList samples = new IndexedSampleList("sample1");
            final AlleleList<Allele> alleles = new IndexedAlleleList<>(vc.getAlleles());
            final ReadLikelihoods<Allele> readLikelihoods = new ReadLikelihoods<>(samples, alleles, sampleReadMap);
            LikelihoodMatrix<Allele> sampleMatrix = readLikelihoods.sampleMatrix(0);
            List<String> readAttributeList = new ArrayList<>();
            for (GATKRead read : vcReadsMap.get(vc)) {
                final int readIndex = sampleMatrix.indexOfRead(read);
                String attribute = contig + ":" + vc.getStart() + "=";
                if (vc == vcs.get(1)) {
                    attribute += supportedAlleleMap.get(read);
                } else {
                    attribute += "1";
                }
                readAttributeList.add(attribute);
                for (Allele allele : vc.getAlleles()) {
                    final int alleleIndex = sampleMatrix.indexOfAllele(allele);
                    if (vc == vcs.get(1)) {
                        if (vc.getAlleleIndex(allele) == supportedAlleleMap.get(read)) {
                            sampleMatrix.set(alleleIndex, readIndex, -1.0);
                        } else {
                            sampleMatrix.set(alleleIndex, readIndex, -8.0);
                        }
                    } else {
                        if (vc.getAlleleIndex(allele) == 1) {
                            sampleMatrix.set(alleleIndex, readIndex, -1.0);
                        } else {
                            sampleMatrix.set(alleleIndex, readIndex, -8.0);
                        }
                    }
                }
            }
            readLikelihoodsList.add(readLikelihoods);
            readAttributeListList.add(readAttributeList);
        }
        return new Object[][]{{readLikelihoodsList, vcs, readAttributeListList}};
    }

    private Haplotype applyVariant(final VariantContext vc, final Haplotype refHaplotype, final int altIndex) {
        Haplotype retHaplotype = refHaplotype.insertAllele(vc.getReference(), vc.getAlternateAllele(altIndex), vc.getStart() - (int) refHaplotype.getStartPosition(), vc.getStart());
        if (retHaplotype == null) {
            return refHaplotype;
        }
        retHaplotype.setGenomeLocation(refHaplotype.getGenomeLocation());
        final Cigar cigar = new Cigar(Arrays.asList(new CigarElement(retHaplotype.length(), CigarOperator.M)));
        retHaplotype.setCigar(cigar);
        return retHaplotype;
    }

    @Test(dataProvider = "testAnnotateReadLikelihoodsWithSupportedAllelesDataProvider")
    public void testAnnotateReadLikelihoodsWithSupportedAlleles(List<ReadLikelihoods<Allele>> readLikelihoodsList, final List<VariantContext> vcs, final List<List<String>> readAttributeListList) {
        for (int i = 0; i < readLikelihoodsList.size(); i++) {
            ReadLikelihoods<Allele> readLikelihoods = readLikelihoodsList.get(i);
            VariantContext vc = vcs.get(i);
            List<String> readAttributeList = readAttributeListList.get(i);


            List<String> initReadAttributes = new ArrayList<>();
            for (GATKRead read : readLikelihoods.sampleReads(0)) {
                initReadAttributes.add(read.getAttributeAsString(AssemblyBasedCallerUtils.SUPPORTED_ALLELES_TAG));
            }
            AssemblyBasedCallerUtils.annotateReadLikelihoodsWithSupportedAlleles(vc, readLikelihoods);
            for (int j = 0; j < readLikelihoods.sampleReadCount(0); j++) {
                GATKRead read = readLikelihoods.sampleReads(0).get(j);

                String expectedAttribute = (initReadAttributes.get(j) != null ? initReadAttributes.get(j) + ", " : "") + readAttributeList.get(j);
                Assert.assertEquals(read.getAttributeAsString(AssemblyBasedCallerUtils.SUPPORTED_ALLELES_TAG), expectedAttribute);
            }
        }
    }

    @Test(dataProvider = "getVariantContextsFromGivenAlleles")
    public void testGetVariantContextsFromGivenAlleles(final int loc,
                                                       final List<VariantContext> activeAllelesToGenotype,
                                                       final List<VariantContext> expectedVcsAtThisLocation) {

        final List<VariantContext> vcsAtThisPosition = AssemblyBasedCallerUtils.getVariantContextsFromGivenAlleles(loc, activeAllelesToGenotype, true);
        Assert.assertEquals(vcsAtThisPosition.size(), expectedVcsAtThisLocation.size());
        for (int i = 0; i < expectedVcsAtThisLocation.size(); i++) {
            VariantContextTestUtils.assertVariantContextsAreEqual(vcsAtThisPosition.get(i), expectedVcsAtThisLocation.get(i), new ArrayList<>());
            Assert.assertEquals(vcsAtThisPosition.get(i).getSource(), expectedVcsAtThisLocation.get(i).getSource());
        }
    }

    @DataProvider(name = "getVariantContextsFromGivenAlleles")
    public Object[][] getVcsAtThisLocationFromGivenAllelesData() {
        final List<Object[]> tests = new ArrayList<>();

        tests.add(new Object[]{1000, new ArrayList<>(), new ArrayList<>()});

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

        final Haplotype sameLocDelHap1 = new Haplotype("AAAAAAAGAAA".getBytes());
        final List<Allele> sameLocDelAlleles1 = Arrays.asList(Allele.create("GTT", true), Allele.create("G"));
        final VariantContext sameLocDelVc1 = new VariantContextBuilder("a", "20", 10093568, 10093570, sameLocDelAlleles1).make();
        sameLocDelHap1.setEventMap(new EventMap(Arrays.asList(sameLocDelVc1)));

        final Haplotype sameLocDelHap2 = new Haplotype("AAAAAAAGTAAA".getBytes());
        final List<Allele> sameLocDelAlleles2 = Arrays.asList(Allele.create("GT", true), Allele.create("G"));
        final VariantContext sameLocDelVc2 = new VariantContextBuilder("a", "20", 10093568, 10093569, sameLocDelAlleles2).make();
        sameLocDelHap2.setEventMap(new EventMap(Arrays.asList(sameLocDelVc2)));

        final Haplotype sameLocInsHap1 = new Haplotype("AAAAAAAGTTTAAA".getBytes());
        final List<Allele> sameLocInsAlleles1 = Arrays.asList(Allele.create("G", true), Allele.create("GT"));
        final VariantContext sameLocInsVc1 = new VariantContextBuilder("a", "20", 10093568, 10093568, sameLocInsAlleles1).make();
        sameLocInsHap1.setEventMap(new EventMap(Arrays.asList(sameLocInsVc1)));

        final VariantContextBuilder deletionVCBuilderWithGts = new VariantContextBuilder("a", "20", 995, 1005, deletionAlleles)
                .genotypes(new GenotypeBuilder("TEST", Arrays.asList(deletionAlleles.get(0), deletionAlleles.get(1))).make());
        final VariantContext deletionVcWithGts = deletionVCBuilderWithGts.make();

        tests.add(new Object[]{1000, Arrays.asList(snpVc), Arrays.asList(snpVCBuilder.source("Comp0Allele0").make())});
        tests.add(new Object[]{995, Arrays.asList(deletionVc), Arrays.asList(deletionVCBuilder.source("Comp0Allele0").make())});
        tests.add(new Object[]{1000, Arrays.asList(deletionVc), Arrays.asList(deletionVCBuilder.source("Comp0Allele0").make())});
        tests.add(new Object[]{1000, Arrays.asList(deletionVc, snpVc),
                Arrays.asList(deletionVCBuilder.source("Comp0Allele0").make(), snpVCBuilder.source("Comp1Allele0").make())});
        tests.add(new Object[]{1000, Arrays.asList(deletionVc, deletionVcNoSpan), Arrays.asList(deletionVCBuilder.source("Comp0Allele0").make())});
        tests.add(new Object[]{1000, Arrays.asList(deletionVc, deletionVcFalseDuplicate, deletionVcNoSpan),
                Arrays.asList(deletionVCBuilder.source("Comp0Allele0").make(), deletionFalseDuplicateBuilder.source("Comp1Allele0").make())});

        tests.add(new Object[]{1000, Arrays.asList(deletionVcWithGts, snpVc),
                Arrays.asList(deletionVCBuilder.source("Comp0Allele0").make(), snpVCBuilder.source("Comp1Allele0").make())});

        return tests.toArray(new Object[][]{});
    }

    @DataProvider(name = "getVariantContextsFromActiveHaplotypes")
    public Object[][] getVariantContextsFromActiveHaplotypesData() {
        final List<Object[]> tests = new ArrayList<>();

        tests.add(new Object[]{new ArrayList<>(), 1000, new ArrayList<>()});

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

        tests.add(new Object[]{Arrays.asList(snpHaplotype), 1000, Arrays.asList(snpVc)});
        tests.add(new Object[]{Arrays.asList(snpHaplotype, snpHaplotypeDuplicate), 1000, Arrays.asList(snpVc)});
        tests.add(new Object[]{Arrays.asList(deletionHaplotype), 995, Arrays.asList(deletionVc)});
        tests.add(new Object[]{Arrays.asList(deletionHaplotype), 1000, Arrays.asList(deletionVc)});
        tests.add(new Object[]{Arrays.asList(deletionHaplotype, deletionHaplotypeNoSpan), 1000, Arrays.asList(deletionVc)});
        tests.add(new Object[]{Arrays.asList(deletionHaplotype, deletionHaplotypeFalseDuplicate, deletionHaplotypeNoSpan), 1000, Arrays.asList(deletionVc, deletionVcFalseDuplicate)});

        tests.add(new Object[]{Arrays.asList(deletionHaplotype, snpHaplotype), 1000, Arrays.asList(deletionVc, snpVc)});

        final Haplotype sameLocDelHap1 = new Haplotype("AAAAAAAGAAA".getBytes());
        final List<Allele> sameLocDelAlleles1 = Arrays.asList(Allele.create("GTT", true), Allele.create("G"));
        final VariantContext sameLocDelVc1 = new VariantContextBuilder("a", "20", 10093568, 10093570, sameLocDelAlleles1).make();
        sameLocDelHap1.setEventMap(new EventMap(Arrays.asList(sameLocDelVc1)));

        final Haplotype sameLocDelHap2 = new Haplotype("AAAAAAAGTAAA".getBytes());
        final List<Allele> sameLocDelAlleles2 = Arrays.asList(Allele.create("GT", true), Allele.create("G"));
        final VariantContext sameLocDelVc2 = new VariantContextBuilder("a", "20", 10093568, 10093569, sameLocDelAlleles2).make();
        sameLocDelHap2.setEventMap(new EventMap(Arrays.asList(sameLocDelVc2)));

        final Haplotype sameLocInsHap1 = new Haplotype("AAAAAAAGTTTAAA".getBytes());
        final List<Allele> sameLocInsAlleles1 = Arrays.asList(Allele.create("G", true), Allele.create("GT"));
        final VariantContext sameLocInsVc1 = new VariantContextBuilder("a", "20", 10093568, 10093568, sameLocInsAlleles1).make();
        sameLocInsHap1.setEventMap(new EventMap(Arrays.asList(sameLocInsVc1)));

        tests.add(new Object[]{Arrays.asList(sameLocDelHap1, sameLocDelHap2, sameLocInsHap1), 10093568, Arrays.asList(sameLocDelVc1, sameLocDelVc2, sameLocInsVc1)});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "getVariantContextsFromActiveHaplotypes")
    public void testGetVariantContextsFromActiveHaplotypes(final List<Haplotype> haplotypes,
                                                           final int loc,
                                                           final List<VariantContext> expectedVcsAtThisLocation) {

        final List<VariantContext> vcsAtThisPosition = AssemblyBasedCallerUtils.getVariantContextsFromActiveHaplotypes(loc, haplotypes, true);
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

        final Haplotype deletionHaplotype2 = new Haplotype("ACTGGTCAGGTCAAGGTCA".getBytes());
        final List<Allele> deletionAlleles2 = Arrays.asList(Allele.create("ACTGGTCAACTCT", true), Allele.create("A"));
        final VariantContextBuilder deletionVCBuilder2 = new VariantContextBuilder("b", "20", 995, 1007, deletionAlleles2);
        final VariantContext deletionVc2 = deletionVCBuilder2.make();
        deletionHaplotype2.setEventMap(new EventMap(Arrays.asList(deletionVc2)));

        final VariantContext spandDelVc2 = new VariantContextBuilder("b", "20", 1000, 1000, Arrays.asList(refAllele, Allele.SPAN_DEL)).make();

        final Haplotype deletionStartingAtLocHaplotype = new Haplotype("ACTGGTCAGGTCAAGGTCA".getBytes());
        final Allele deletionStartingAtLocRefAllele = Allele.create("ACTGGTCAACTCT", true);
        final List<Allele> deletionStartingAtLocAlleles = Arrays.asList(deletionStartingAtLocRefAllele, Allele.create("A"));
        final VariantContextBuilder deletionStartingAtLocVCBuilder = new VariantContextBuilder("b", "20", 1000, 1012, deletionStartingAtLocAlleles);
        final VariantContext deletionStartingAtLocVc = deletionStartingAtLocVCBuilder.make();
        deletionStartingAtLocHaplotype.setEventMap(new EventMap(Arrays.asList(deletionStartingAtLocVc)));

        final Allele remappedSNPAllele = Allele.create("GCTGGTCAACTCT");
        final VariantContext mergedSnpAndDelStartingAtLocVC = new VariantContextBuilder("a", "20", 1000, 1012,
                Arrays.asList(deletionStartingAtLocRefAllele,
                        Allele.create("A"), // for the deletion,
                        remappedSNPAllele // for the SNP
                )).make();


        final List<VariantContext> emptyGivenAllelesList = new ArrayList<>();

        final VariantContext mergedSnpAndDelVC = new VariantContextBuilder("a", "20", 1000, 1000,
                Arrays.asList(refAllele,
                        Allele.SPAN_DEL,
                        Allele.create("G"))).make();



        final List<Object[]> tests = new ArrayList<>();
        tests.add(new Object[]{
                snpVc,
                snpVc.getStart(),
                Arrays.asList(snpHaplotype, refHaplotype),
                emptyGivenAllelesList,
                Maps.asMap(new HashSet<>(snpAlleles),
                        (key) -> {
                            if (snpAlleles.get(1).equals(key)) return Arrays.asList(snpHaplotype);
                            return Arrays.asList(refHaplotype);
                        })
        });
        tests.add(new Object[]{
                mergedSnpAndDelVC,
                mergedSnpAndDelVC.getStart(),
                Arrays.asList(snpHaplotype, refHaplotype, deletionHaplotype),
                emptyGivenAllelesList,
                Maps.asMap(new HashSet<>(mergedSnpAndDelVC.getAlleles()),
                        (key) -> {
                            if (snpAlleles.get(1).equals(key)) return Arrays.asList(snpHaplotype);
                            if (Allele.SPAN_DEL.equals(key)) return Arrays.asList(deletionHaplotype);
                            return Arrays.asList(refHaplotype);
                        })
        });
        // includes a SNP haplotype not present in events at this loc (which might happen in GGA mode)
        tests.add(new Object[]{
                snpVc,
                snpVc.getStart(),
                Arrays.asList(snpHaplotype, refHaplotype, snpHaplotypeNotPresentInEventsAtThisLoc),
                Arrays.asList(snpVc),
                Maps.asMap(new HashSet<>(snpVc.getAlleles()),
                        (key) -> {
                            if (snpAlleles.get(1).equals(key)) return Arrays.asList(snpHaplotype);
                            return Arrays.asList(refHaplotype);
                        })
        });

        // two spanning deletions, no given alleles -> both dels should be in event map for span del
        tests.add(new Object[]{
                mergedSnpAndDelVC,
                snpVc.getStart(),
                Arrays.asList(snpHaplotype, refHaplotype, deletionHaplotype, deletionHaplotype2),
                emptyGivenAllelesList,
                Maps.asMap(new HashSet<>(mergedSnpAndDelVC.getAlleles()),
                        (key) -> {
                            if (snpAlleles.get(1).equals(key)) return Arrays.asList(snpHaplotype);
                            if (Allele.SPAN_DEL.equals(key)) return Arrays.asList(deletionHaplotype, deletionHaplotype2);
                            return Arrays.asList(refHaplotype);
                        })
        });

        // two spanning deletions, one in given alleles -> only the matching deletion should be in the event map for the span del
        tests.add(new Object[]{
                mergedSnpAndDelVC,
                snpVc.getStart(),
                Arrays.asList(snpHaplotype, refHaplotype, deletionHaplotype, deletionHaplotype2),
                Arrays.asList(snpVc, deletionVc2),
                Maps.asMap(new HashSet<>(mergedSnpAndDelVC.getAlleles()),
                        (key) -> {
                            if (snpAlleles.get(1).equals(key)) return Arrays.asList(snpHaplotype);
                            if (Allele.SPAN_DEL.equals(key)) return Arrays.asList(deletionHaplotype2);
                            return Arrays.asList(refHaplotype);
                        })
        });

        // A deletion starting at the loc in the given alleles, the snp not in the given alleles
        tests.add(new Object[]{
                deletionStartingAtLocVc,
                deletionStartingAtLocVc.getStart(),
                Arrays.asList(snpHaplotype, refHaplotype, deletionStartingAtLocHaplotype),
                Arrays.asList(deletionStartingAtLocVc),
                Maps.asMap(new HashSet<>(deletionStartingAtLocVc.getAlleles()),
                        (key) -> {
                            if (deletionStartingAtLocAlleles.get(1).equals(key)) return Arrays.asList(deletionStartingAtLocHaplotype);
                            return Arrays.asList(refHaplotype);
                        })
        });

        // A deletion starting at the loc not in the given alleles, the snp in the given alleles
        tests.add(new Object[]{
                snpVc,
                snpVc.getStart(),
                Arrays.asList(snpHaplotype, refHaplotype, deletionStartingAtLocHaplotype),
                Arrays.asList(snpVc),
                Maps.asMap(new HashSet<>(snpVc.getAlleles()),
                        (key) -> {
                            if (snpAlleles.get(1).equals(key)) return Arrays.asList(snpHaplotype);
                            return Arrays.asList(refHaplotype);
                        })
        });

        // A deletion starting at the loc and the SNP in the given alleles
        tests.add(new Object[]{
                mergedSnpAndDelStartingAtLocVC,
                snpVc.getStart(),
                Arrays.asList(snpHaplotype, refHaplotype, deletionStartingAtLocHaplotype),
                Arrays.asList(deletionStartingAtLocVc, snpVc),
                Maps.asMap(new HashSet<>(mergedSnpAndDelStartingAtLocVC.getAlleles()),
                        (key) -> {
                            if (deletionStartingAtLocAlleles.get(1).equals(key)) return Arrays.asList(deletionStartingAtLocHaplotype);
                            if (remappedSNPAllele.equals(key)) return Arrays.asList(snpHaplotype);
                            return Arrays.asList(refHaplotype);
                        })
        });


        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "getEventMapper")
    public void testGetEventMapper(final VariantContext mergedVc,
                                   final int loc,
                                   final List<Haplotype> haplotypes,
                                   final List<VariantContext> activeAllelesToGenotype,
                                   final Map<Allele, List<Haplotype>> expectedEventMap) {
        final Map<Allele, List<Haplotype>> actualEventMap = AssemblyBasedCallerUtils.createAlleleMapper(mergedVc, loc, haplotypes, activeAllelesToGenotype);
        Assert.assertEquals(actualEventMap.size(), expectedEventMap.size());
        for (final Allele key : actualEventMap.keySet()) {
            Assert.assertTrue(expectedEventMap.containsKey(key), "Got unexpected allele " + key + " with values " + actualEventMap.get(key));
            Assert.assertEquals(actualEventMap.get(key), expectedEventMap.get(key), "Lists don't match for key " + key);
        }

        for (final Allele key : expectedEventMap.keySet()) {
            Assert.assertTrue(actualEventMap.containsKey(key), "Didn't get back allele " + key);
        }
    }

    @DataProvider(name = "ConstructPhaseGroupsProvider")
    public Object[][] makeConstructPhaseGroupsData() {
        List<Object[]> tests = new ArrayList<>();

        final Allele ref = Allele.create("A", true);
        final Allele altC = Allele.create("C", false);

        final Genotype g1 = new GenotypeBuilder().alleles(Arrays.asList(ref, altC)).make();
        final VariantContext vc1 = new VariantContextBuilder().chr("20").start(1).stop(1).alleles(Arrays.asList(ref, altC)).genotypes(g1).make();
        final Genotype g2 = new GenotypeBuilder().alleles(Arrays.asList(ref, altC)).make();
        final VariantContext vc2 = new VariantContextBuilder().chr("20").start(2).stop(2).alleles(Arrays.asList(ref, altC)).genotypes(g2).make();
        final Genotype g3 = new GenotypeBuilder().alleles(Arrays.asList(ref, altC)).make();
        final VariantContext vc3 = new VariantContextBuilder().chr("20").start(3).stop(3).alleles(Arrays.asList(ref, altC)).genotypes(g3).make();
        final List<VariantContext> calls = Arrays.asList(vc1, vc2, vc3);

        // test no phased variants, empty map
        final Map<VariantContext, Pair<Integer, String>> nonePhased1 = new HashMap<>();
        tests.add(new Object[]{calls, nonePhased1, 0, 0, 0});

        // test no phased variants, full map, exception expected
        final Map<VariantContext, Pair<Integer, String>> nonePhased2 = new HashMap<>();
        nonePhased2.put(vc1, Pair.of(0, "0/1"));
        nonePhased2.put(vc2, Pair.of(1, "0/1"));
        nonePhased2.put(vc3, Pair.of(2, "0/1"));
        tests.add(new Object[]{calls, nonePhased2, 3, -1, -1});

        // test 2 phased variants
        final Map<VariantContext, Pair<Integer, String>> twoPhased = new HashMap<>();
        twoPhased.put(vc1, Pair.of(0, "0/1"));
        twoPhased.put(vc2, Pair.of(0, "0/1"));
        tests.add(new Object[]{calls, twoPhased, 1, 1, 2});

        // test all phased variants
        final Map<VariantContext, Pair<Integer, String>> allPhased = new HashMap<>();
        allPhased.put(vc1, Pair.of(0, "0/1"));
        allPhased.put(vc2, Pair.of(0, "0/1"));
        allPhased.put(vc3, Pair.of(0, "0/1"));
        tests.add(new Object[]{calls, allPhased, 1, 1, 3});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider="ConstructPhaseGroupsProvider")
    public void testConstructPhaseGroups(final List<VariantContext> calls,
                                         final Map<VariantContext, Pair<Integer, String>> phaseMap,
                                         final int endIndex,
                                         final int expectedNumGroups,
                                         final int expectedGroupSize) {
        final List<VariantContext> actualPhasedCalls;
        try {
            actualPhasedCalls = AssemblyBasedCallerUtils.constructPhaseGroups(calls, phaseMap, endIndex);
        } catch (IllegalStateException e) {
            Assert.assertEquals(-1, expectedNumGroups);
            return;
        }

        final Set<String> uniqueGroups = new HashSet<>();
        int counter = 0;
        for ( final VariantContext call : actualPhasedCalls ) {
            for ( final Genotype g : call.getGenotypes() ) {
                if ( g.hasExtendedAttribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_ID_KEY) ) {
                    uniqueGroups.add(g.getExtendedAttribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_ID_KEY).toString());
                    counter++;
                }
            }
        }

        Assert.assertEquals(uniqueGroups.size(), expectedNumGroups);
        Assert.assertEquals(counter, expectedGroupSize);
    }

    @DataProvider(name = "ConstructPhaseSetMappingProvider")
    public Object[][] makeConstructPhaseSetMappingData() {
        List<Object[]> tests = new ArrayList<Object[]>();

        final Allele ref = Allele.create("A", true);
        final Allele altC = Allele.create("C", false);
        final Allele altT = Allele.create("T", false);

        final VariantContext vc1 = new VariantContextBuilder().chr("20").start(1).stop(1).alleles(Arrays.asList(ref, altC)).make();
        final VariantContext vc2 = new VariantContextBuilder().chr("20").start(2).stop(2).alleles(Arrays.asList(ref, altC)).make();
        final VariantContext vc3 = new VariantContextBuilder().chr("20").start(3).stop(3).alleles(Arrays.asList(ref, altT)).make();
        final VariantContext vc4 = new VariantContextBuilder().chr("20").start(4).stop(4).alleles(Arrays.asList(ref, altC)).make();
        final List<VariantContext> calls = Arrays.asList(vc2, vc3, vc4);

        final Haplotype pos1 = new Haplotype("CAAAA".getBytes());
        pos1.setEventMap(new EventMap(Arrays.asList(vc1)));
        pos1.getEventMap().put(1, vc1);
        final Haplotype pos2 = new Haplotype("ACAAA".getBytes());
        pos2.setEventMap(new EventMap(Arrays.asList(vc2)));
        pos2.getEventMap().put(2, vc2);
        final Haplotype pos3 = new Haplotype("AACAA".getBytes());
        pos3.setEventMap(new EventMap(Arrays.asList(vc3)));
        pos3.getEventMap().put(3, vc3);
        final Haplotype pos4 = new Haplotype("AAACA".getBytes());
        pos4.setEventMap(new EventMap(Arrays.asList(vc4)));
        pos4.getEventMap().put(4, vc4);
        final Haplotype pos24 = new Haplotype("ACACA".getBytes());
        pos24.setEventMap(new EventMap(Arrays.asList(vc2, vc4)));
        pos24.getEventMap().put(2, vc2);
        pos24.getEventMap().put(4, vc4);
        final Haplotype pos34 = new Haplotype("AACCA".getBytes());
        pos34.setEventMap(new EventMap(Arrays.asList(vc3, vc4)));
        pos34.getEventMap().put(3, vc3);
        pos34.getEventMap().put(4, vc4);
        final Haplotype pos234 = new Haplotype("ACCCA".getBytes());
        pos234.setEventMap(new EventMap(Arrays.asList(vc2, vc3, vc4)));
        pos234.getEventMap().put(2, vc2);
        pos234.getEventMap().put(3, vc3);
        pos234.getEventMap().put(4, vc4);

        final Map<VariantContext, Set<Haplotype>> haplotypeMap = new HashMap<>();

        // test no phased variants #1
        final Set<Haplotype> haplotypes2 = new HashSet<>();
        haplotypes2.add(pos2);
        haplotypeMap.put(vc2, haplotypes2);
        tests.add(new Object[]{Arrays.asList(vc2), new HashMap<>(haplotypeMap), 2, 0, 0, 0, 0});

        // test no phased variants #2
        final Set<Haplotype> haplotypes3 = new HashSet<>();
        haplotypes3.add(pos3);
        haplotypeMap.put(vc3, haplotypes3);
        tests.add(new Object[]{Arrays.asList(vc2, vc3), new HashMap<>(haplotypeMap), 3, 0, 0, 0, 0});

        // test opposite phase
        tests.add(new Object[]{Arrays.asList(vc2, vc3), new HashMap<>(haplotypeMap), 2, 2, 1, 1, 1});

        // test no phased variants #3
        final Set<Haplotype> haplotypes4 = new HashSet<>();
        haplotypes4.add(pos4);
        haplotypeMap.put(vc4, haplotypes4);
        tests.add(new Object[]{calls, new HashMap<>(haplotypeMap), 3, 0, 0, 0, 0});

        // test mixture
        final Set<Haplotype> haplotypes24 = new HashSet<>();
        haplotypes24.add(pos24);
        haplotypeMap.put(vc2, haplotypes24);
        haplotypeMap.put(vc4, haplotypes24);
        tests.add(new Object[]{calls, new HashMap<>(haplotypeMap), 2, 3, 1, 2, 1});

        // test 2 hets
        haplotypeMap.remove(vc3);
        tests.add(new Object[]{Arrays.asList(vc2, vc4), new HashMap<>(haplotypeMap), 1, 2, 1, 2, 0});

        // test 2 with opposite phase
        final Set<Haplotype> haplotypes1 = new HashSet<>();
        haplotypes1.add(pos1);
        haplotypeMap.put(vc1, haplotypes1);
        tests.add(new Object[]{Arrays.asList(vc1, vc2, vc4), new HashMap<>(haplotypeMap), 2, 3, 1, 1, 2});

        // test homs around a het
        final Set<Haplotype> haplotypes2hom = new HashSet<>();
        haplotypes2hom.add(pos24);
        haplotypes2hom.add(pos234);
        final Set<Haplotype> haplotypes4hom = new HashSet<>();
        haplotypes4hom.add(pos24);
        haplotypes4hom.add(pos234);
        final Set<Haplotype> haplotypes3het = new HashSet<>();
        haplotypes3het.add(pos234);
        haplotypeMap.put(vc2, haplotypes2hom);
        haplotypeMap.put(vc3, haplotypes3het);
        haplotypeMap.put(vc4, haplotypes4hom);
        tests.add(new Object[]{calls, new HashMap<>(haplotypeMap), 2, 3, 1, 3, 0});

        // test hets around a hom
        final Set<Haplotype> haplotypes2het = new HashSet<>();
        haplotypes2het.add(pos234);
        final Set<Haplotype> haplotypes4het = new HashSet<>();
        haplotypes4het.add(pos234);
        final Set<Haplotype> haplotypes3hom = new HashSet<>();
        haplotypes3hom.add(pos3);
        haplotypes3hom.add(pos234);
        haplotypeMap.put(vc2, haplotypes2het);
        haplotypeMap.put(vc3, haplotypes3hom);
        haplotypeMap.put(vc4, haplotypes4het);
        tests.add(new Object[]{calls, new HashMap<>(haplotypeMap), 2, 3, 1, 3, 0});

        // test no phased variants around a hom
        final Set<Haplotype> haplotypes2incomplete = new HashSet<>();
        haplotypes2incomplete.add(pos24);
        final Set<Haplotype> haplotypes3incomplete = new HashSet<>();
        haplotypes3incomplete.add(pos34);
        final Set<Haplotype> haplotypes4complete = new HashSet<>();
        haplotypes4complete.add(pos24);
        haplotypes4complete.add(pos34);
        haplotypes4complete.add(pos234);
        haplotypeMap.put(vc2, haplotypes2incomplete);
        haplotypeMap.put(vc3, haplotypes3incomplete);
        haplotypeMap.put(vc4, haplotypes4complete);
        tests.add(new Object[]{calls, new HashMap<>(haplotypeMap), 0, 0, 0, 0, 0});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider="ConstructPhaseSetMappingProvider")
    public void testConstructPhaseSetMapping(final List<VariantContext> calls,
                                             final Map<VariantContext, Set<Haplotype>> haplotypeMap,
                                             final int totalHaplotypes,
                                             final int expectedMapSize,
                                             final int expectedNumGroups,
                                             final int expectedNum01,
                                             final int expectedNum10) {
        final Map<VariantContext, Pair<Integer, String>> actualPhaseSetMapping = new HashMap<>();
        final int actualNumGroups = AssemblyBasedCallerUtils.constructPhaseSetMapping(calls, haplotypeMap, totalHaplotypes, actualPhaseSetMapping);
        Assert.assertEquals(actualNumGroups, expectedNumGroups);
        Assert.assertEquals(actualPhaseSetMapping.size(), expectedMapSize);

        int num01 = 0, num10 = 0;
        for ( final Pair<Integer, String> phase : actualPhaseSetMapping.values() ) {
            if ( phase.getRight().equals("0|1") )
                num01++;
            else if ( phase.getRight().equals("1|0") )
                num10++;
        }
        Assert.assertEquals(num01, expectedNum01);
        Assert.assertEquals(num10, expectedNum10);
    }

    @DataProvider(name = "CreateHaplotypeMappingProvider")
    public Object[][] makeCreateHaplotypeMappingData() {
        List<Object[]> tests = new ArrayList<Object[]>();

        final Set<Haplotype> haplotypes = new HashSet<>();
        final Allele ref = Allele.create("A", true);
        final Allele altC = Allele.create("C", false);
        final Allele altT = Allele.create("T", false);

        final Haplotype AtoC1 = new Haplotype("AACAA".getBytes());
        final VariantContext vc1 = new VariantContextBuilder().chr("20").start(3).stop(3).alleles(Arrays.asList(ref, altC)).make();
        AtoC1.setEventMap(new EventMap(Arrays.asList(vc1)));
        AtoC1.getEventMap().put(3, vc1);
        haplotypes.add(AtoC1);

        final Haplotype AtoC2 = new Haplotype("AAACA".getBytes());
        final VariantContext vc2 = new VariantContextBuilder().chr("20").start(4).stop(4).alleles(Arrays.asList(ref, altT)).make();
        AtoC2.setEventMap(new EventMap(Arrays.asList(vc2)));
        AtoC2.getEventMap().put(4, vc2);
        haplotypes.add(AtoC2);

        tests.add(new Object[]{vc1, haplotypes, AtoC1});
        tests.add(new Object[]{vc2, haplotypes, AtoC2});
        tests.add(new Object[]{new VariantContextBuilder().chr("20").start(1).stop(1).alleles(Arrays.asList(ref, altT)).make(), haplotypes, null});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider="CreateHaplotypeMappingProvider")
    public void testCreateHaplotypeMapping(final VariantContext vc, final Set<Haplotype> haplotypes, final Haplotype expected) {
        final Map<VariantContext, Set<Haplotype>> mapping = AssemblyBasedCallerUtils.constructHaplotypeMapping(Arrays.asList(vc), haplotypes);
        final Set<Haplotype> actual = mapping.get(vc);
        if ( expected == null )
            Assert.assertTrue(actual.isEmpty(), actual.toString());
        else {
            Assert.assertEquals(actual.size(), 1);
            Assert.assertEquals(actual.iterator().next(), expected);
        }
    }

    @Test
    public void testAddGivenHaplotypes() {
        final int assemblyRegionStart = 1;
        final int maxMnpDistance = 0;
        final SmithWatermanAligner aligner = SmithWatermanAligner.getAligner(SmithWatermanAligner.Implementation.FASTEST_AVAILABLE);
        final AssemblyResultSet assemblyResultSet = new AssemblyResultSet();

        final Haplotype refHaplotype = new Haplotype("AAAACCCCGGGGTTTT".getBytes(), true);
        final byte[] fullReferenceWithPadding = ("A" + refHaplotype.getBaseString()).getBytes();
        refHaplotype.setAlignmentStartHapwrtRef(assemblyRegionStart);
        refHaplotype.setCigar(new Cigar(Collections.singletonList(new CigarElement(refHaplotype.length(), CigarOperator.M))));
        refHaplotype.setGenomeLocation(new SimpleInterval("chr", assemblyRegionStart, assemblyRegionStart + refHaplotype.length()));
        assemblyResultSet.setPaddedReferenceLoc(new SimpleInterval("chr", 1, assemblyRegionStart + refHaplotype.length()));
        assemblyResultSet.add(refHaplotype);
        assemblyResultSet.setFullReferenceWithPadding(fullReferenceWithPadding);

        // add a SNP
        final VariantContext givenVC = new VariantContextBuilder("test", "chr", 2, 2,
                Arrays.asList(Allele.create((byte) 'A', true), Allele.create((byte) 'C', false))).make();

        AssemblyBasedCallerUtils.addGivenHaplotypes(assemblyRegionStart, Collections.singletonList(givenVC), maxMnpDistance,
                aligner, refHaplotype, assemblyResultSet);
        Assert.assertEquals(assemblyResultSet.getHaplotypeCount(), 2);
        Assert.assertEquals(assemblyResultSet.getHaplotypeList().get(1).getBaseString(), "ACAACCCCGGGGTTTT");


        // adding the same VC should have no effect
        AssemblyBasedCallerUtils.addGivenHaplotypes(assemblyRegionStart, Collections.singletonList(givenVC), maxMnpDistance,
                aligner, refHaplotype, assemblyResultSet);
        Assert.assertEquals(assemblyResultSet.getHaplotypeCount(), 2);

        // add another SNP
        final VariantContext givenVC2 = new VariantContextBuilder("test", "chr", 5, 5,
                Arrays.asList(Allele.create((byte) 'C', true), Allele.create((byte) 'G', false))).make();
        AssemblyBasedCallerUtils.addGivenHaplotypes(assemblyRegionStart, Collections.singletonList(givenVC2), maxMnpDistance,
                aligner, refHaplotype, assemblyResultSet);
        // SNP is not found in existing variation, so it's added to the ref and the first SNP
        Assert.assertEquals(assemblyResultSet.getHaplotypeCount(), 4);
        Assert.assertEquals(assemblyResultSet.getHaplotypeList().get(2).getBaseString(), "AAAAGCCCGGGGTTTT");
        Assert.assertEquals(assemblyResultSet.getHaplotypeList().get(3).getBaseString(), "ACAAGCCCGGGGTTTT");

        // add a deletion that overlaps the second SNP.  This variant gets added to the ref and first SNP haplotypes but not either
        // haplotype that contains the overlapping 2nd SNP
        final VariantContext givenVC3 = new VariantContextBuilder("test", "chr", 5, 7,
                Arrays.asList(Allele.create("CCC".getBytes(), true), Allele.create((byte) 'C', false))).make();
        AssemblyBasedCallerUtils.addGivenHaplotypes(assemblyRegionStart, Collections.singletonList(givenVC3), maxMnpDistance,
                aligner, refHaplotype, assemblyResultSet);
        Assert.assertEquals(assemblyResultSet.getHaplotypeCount(), 6);
        Assert.assertEquals(assemblyResultSet.getHaplotypeList().get(4).getBaseString(), "AAAACCGGGGTTTT");
        Assert.assertEquals(assemblyResultSet.getHaplotypeList().get(5).getBaseString(), "ACAACCGGGGTTTT");

        // adding an equivalent deletion should do nothing
        final VariantContext givenVC4 = new VariantContextBuilder("test", "chr", 5, 8,
                Arrays.asList(Allele.create("CCCC".getBytes(), true), Allele.create("CC".getBytes(), false))).make();
        AssemblyBasedCallerUtils.addGivenHaplotypes(assemblyRegionStart, Collections.singletonList(givenVC4), maxMnpDistance,
                aligner, refHaplotype, assemblyResultSet);
        Assert.assertEquals(assemblyResultSet.getHaplotypeCount(), 6);

        // finally, add a haplotype with two new phased SNPs, after which adding an allele with one of these SNPs does nothing
        final Haplotype phasedHaplotype = new Haplotype("AAAACCTCGAGGTTTT".getBytes(), false);
        phasedHaplotype.setAlignmentStartHapwrtRef(assemblyRegionStart);
        phasedHaplotype.setCigar(new Cigar(Collections.singletonList(new CigarElement(refHaplotype.length(), CigarOperator.M))));
        phasedHaplotype.setGenomeLocation(new SimpleInterval("chr", assemblyRegionStart, assemblyRegionStart + refHaplotype.length()));
        assemblyResultSet.add(phasedHaplotype);
        assemblyResultSet.regenerateVariationEvents(maxMnpDistance);

        Assert.assertEquals(assemblyResultSet.getHaplotypeCount(), 7);


        final VariantContext givenVC5 = new VariantContextBuilder("test", "chr", 8, 8,
                Arrays.asList(Allele.create((byte) 'C', true), Allele.create((byte) 'T', false))).make();
        AssemblyBasedCallerUtils.addGivenHaplotypes(assemblyRegionStart, Collections.singletonList(givenVC5), maxMnpDistance,
                aligner, refHaplotype, assemblyResultSet);
        Assert.assertEquals(assemblyResultSet.getHaplotypeCount(), 7);

    }
}
