package org.broadinstitute.hellbender.tools.walkers.variantutils;

import static org.testng.Assert.*;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.ReferenceMemorySource;
import org.broadinstitute.hellbender.engine.VariantWalker;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.stream.*;

public class LeftAlignAndTrimVariantsUnitTest extends GATKBaseTest{
    final String refBases1 = "ACAGAGCTGACCCTCCCTCCCCTCTCCCAGTGCAACAGCACGGGCGGCGACTGCTTTTACCGAGGCTACACGTCAGGCGTGGCGGCTGTCCAGGACTGGTACCACTTCCACTATGTGGATCTCTGCTGAGGACCAGGAAAGCCAGCACCCGCAGAGACTCTTCCCCAGTGCTCCATACGATCACCATTCTCTGCAGAAGG";
    final List<String> longPieces = Arrays.asList("AAAAAAAAAAAAAAAAAAAAAAAAAAAA","TCTCTCTCTCTCTC"); // where we'll perform tests
    final List<String> strs=Arrays.asList("A","TC"); //
    final List<String> refBasesStrings = longPieces.stream().map(l->refBases1 + l + refBases1).collect(Collectors.toList());

    final List<Integer> contigStops = refBasesStrings.stream().map(String::length).collect(Collectors.toList());
    final List<SAMFileHeader> headers=contigStops.stream().map(c->
            ArtificialReadUtils.createArtificialSamHeader(1, 1, c )).collect(Collectors.toList());
    final String artificialContig = "1";
    final int repeatStart=refBases1.length();
    final List<SimpleInterval> refIntervals=contigStops.stream().map(c->
            new SimpleInterval(artificialContig,1,c)).collect(Collectors.toList());
    final List<ReferenceBases> refBases=IntStream.range(0,longPieces.size()).mapToObj(i->new ReferenceBases(refBasesStrings.get(i).getBytes(),refIntervals.get(i)))
            .collect(Collectors.toList());
    //final ReferenceMemorySource refSource=new ReferenceMemorySource(refBases,header.getSequenceDictionary());
    final List<ReferenceMemorySource> refSources=IntStream.range(0,longPieces.size()).mapToObj(i->new ReferenceMemorySource(refBases.get(i),headers.get(i).getSequenceDictionary()))
            .collect(Collectors.toList());

    @DataProvider(name = "LeftAlignDataProvider")
    public Object[][] LeftAlignTestData() {
        List<Object[]> tests = new ArrayList<Object[]>();
        for (int iLongPiece=0; iLongPiece < longPieces.size();iLongPiece++) {
            final String longPiece=longPieces.get(iLongPiece);
            final ReferenceMemorySource refSource=refSources.get(iLongPiece);
            final String str=strs.get(iLongPiece);
            final int nStrs=longPiece.length()/str.length();
            for (int offset = 1; offset < nStrs; offset++) {
                for (int indelSize = -nStrs + offset; indelSize < nStrs - offset; indelSize++) {
                    if (indelSize == 0) {
                        continue;
                    }
                    final List<Allele> alleles = new ArrayList<Allele>();
                    if (indelSize < 0) { // deletion
                        alleles.add(Allele.create(Utils.dupString(str, Math.abs(indelSize) + 1), true));
                        alleles.add(Allele.create(str, false));
                    } else {
                        alleles.add(Allele.create(str, true));
                        alleles.add(Allele.create(Utils.dupString(str, Math.abs(indelSize) + 1), false));
                    }
                    final SimpleInterval interval = new SimpleInterval(artificialContig, repeatStart + offset*str.length(), repeatStart + offset*str.length());
                    ReferenceContext ref = new ReferenceContext(refSource, interval);

                    final VariantContext vc = new VariantContextBuilder("test", artificialContig, repeatStart + offset*str.length(), repeatStart + offset*str.length() + alleles.get(0).length() - 1, alleles).make();
                    tests.add(new Object[]{vc, ref, offset != 0, repeatStart});
                }
            }
        }
        return tests.toArray(new Object[][]{});

    }

    @Test(dataProvider = "LeftAlignDataProvider")
    public void testLeftAlign(final VariantContext vc, final ReferenceContext ref,final boolean expectRealigned, final int expectedStart) {
        final LeftAlignAndTrimVariants leftAligner=new LeftAlignAndTrimVariants();
        VariantContext vcTrim=GATKVariantContextUtils.trimAlleles(vc, true, true);
        final VariantContext realignedV=leftAligner.leftAlign(vcTrim,ref);
        Assert.assertEquals(realignedV!=vcTrim,expectRealigned);
        Assert.assertEquals(realignedV.getStart(),expectedStart);
    }
}