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

import java.io.UnsupportedEncodingException;
import java.util.*;
import java.util.stream.*;

public class LeftAlignAndTrimVariantsUnitTest extends GATKBaseTest{
    final String refBases1 = "ACAGAGCTGACCCTCCCTCCCCTCTCCCAGTGCAACAGCACGGGCGGCGACTGCTTTTACCGAGGCTACACGTCAGGCGTGGCGGCTGTCCAGGACTGGTACCACTTCCACTATGTGGATCTCTGCTGAGGACCAGGAAAGCCAGCACCCGCAGAGACTCTTCCCCAGTGCTCCATACGATCACCATTCTCTGCAGAAGG";
    final List<String> longPieces = Arrays.asList("AAAAAAAAAAAAAAAAAAAAAAAAAAAA","TCTCTCTCTCTCTC"); // where we'll perform tests
    final List<Integer> strLengths=Arrays.asList(1,2); //
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
    public Object[][] LeftAlignTestData() throws UnsupportedEncodingException {
        List<Object[]> tests = new ArrayList<Object[]>();
        for (int iLongPiece=0; iLongPiece < longPieces.size();iLongPiece++) {
            final String longPiece=longPieces.get(iLongPiece);
            final ReferenceMemorySource refSource=refSources.get(iLongPiece);
            final int strLength=strLengths.get(iLongPiece);
            //final int nStrs=longPiece.length()/strLength; //number of strs in the repeat
            final ReferenceBases theseRefBases=refBases.get(iLongPiece);
            for (int indelIndex = repeatStart-1; indelIndex < repeatStart+longPiece.length(); indelIndex++) {
                for (int indelRepeats = -(longPiece.length()-(indelIndex-repeatStart))/strLength; indelRepeats < 10; indelRepeats++) {
                    if (indelRepeats == 0) {
                        continue;
                    }
                    final List<Allele> alleles = new ArrayList<Allele>();
                    if (indelRepeats < 0) { // deletion
                        byte [] basesRef= new byte[Math.abs(indelRepeats)*strLength+1];
                        byte [] basesAlt = new byte[1];
                        System.arraycopy(theseRefBases,repeatStart,basesRef,0,Math.abs(indelRepeats)*strLength+1);
                        System.arraycopy(theseRefBases,repeatStart,basesAlt,0,1);
                        alleles.add(Allele.create(new String(basesRef, "UTF-8"), true));
                        alleles.add(Allele.create(new String(basesAlt, "UTF-8"), false));
                    } else {
                        alleles.add(Allele.create(str, true));
                        alleles.add(Allele.create(Utils.dupString(str, Math.abs(indelSize) + 1), false));
                    }
                    final SimpleInterval interval = new SimpleInterval(artificialContig, repeatStart + offset*str.length()+1, repeatStart + offset*str.length()+1);
                    ReferenceContext ref = new ReferenceContext(refSource, interval);

                    final VariantContext vc = new VariantContextBuilder("test", artificialContig, repeatStart + offset*str.length()+1, repeatStart + offset*str.length() + alleles.get(0).length(), alleles).make();
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