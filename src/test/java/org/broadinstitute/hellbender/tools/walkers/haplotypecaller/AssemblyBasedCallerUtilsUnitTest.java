package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.samtools.*;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContextBuilder;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.engine.AssemblyRegion;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.genotyper.*;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.apache.commons.lang3.StringUtils;

import java.util.*;
import java.util.stream.Collectors;

import org.broadinstitute.hellbender.GATKBaseTest;
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
}
