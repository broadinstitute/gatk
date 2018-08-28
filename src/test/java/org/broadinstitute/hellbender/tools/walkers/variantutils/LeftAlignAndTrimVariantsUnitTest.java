package org.broadinstitute.hellbender.tools.walkers.variantutils;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.ReferenceMemorySource;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
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

public class LeftAlignAndTrimVariantsUnitTest extends GATKBaseTest {
    final String refBases1 = "ACAGAGCTGACCCTCCCTCCCCTCTCCCAGTGCAACAGCACGGGCGGCGACTGCTTTTACCGAGGCTACACGTCAGGCGTGGCGGCTGTCCAGGACTGGTACCACTTCCACTATGTGGATCTCTGCTGAGGACCAGGAAAGCCAGCACCCGCAGAGACTCTTCCCCAGTGCTCCATACGATCACCATTCTCTGCAGAAGG";
    final List<String> longPieces = Arrays.asList("AAAAAAAAAAAAAAAAAAAAAAAAAAAA", "TCTCTCTCTCTCTC", "CAGTCCAGTCCAGTCCAGTCCAGTCCAGTCCAGTCCAGTCCAGTCCAGTCCAGTC", Utils.dupString("A", 120)); // where we'll perform tests
    final List<Integer> strLengths = Arrays.asList(1, 2, 5, 1); //
    final List<String> refBasesStrings = longPieces.stream().map(l -> refBases1 + l + refBases1).collect(Collectors.toList());

    final List<Integer> contigStops = refBasesStrings.stream().map(String::length).collect(Collectors.toList());
    final List<SAMFileHeader> headers = contigStops.stream().map(c ->
            ArtificialReadUtils.createArtificialSamHeader(1, 1, c)).collect(Collectors.toList());
    final String artificialContig = "1";
    final int repeatStart = refBases1.length();
    final List<SimpleInterval> refIntervals = contigStops.stream().map(c ->
            new SimpleInterval(artificialContig, 1, c)).collect(Collectors.toList());
    final List<ReferenceBases> refBases = IntStream.range(0, longPieces.size()).mapToObj(i -> new ReferenceBases(refBasesStrings.get(i).getBytes(), refIntervals.get(i)))
            .collect(Collectors.toList());
    //final ReferenceMemorySource refSource=new ReferenceMemorySource(refBases,header.getSequenceDictionary());
    final List<ReferenceMemorySource> refSources = IntStream.range(0, longPieces.size()).mapToObj(i -> new ReferenceMemorySource(refBases.get(i), headers.get(i).getSequenceDictionary()))
            .collect(Collectors.toList());

    @DataProvider(name = "LeftAlignDataProvider")
    public Object[][] LeftAlignTestData() throws UnsupportedEncodingException {
        List<Object[]> tests = new ArrayList<Object[]>();
        for (int iLongPiece = 0; iLongPiece < longPieces.size(); iLongPiece++) {
            final String longPiece = longPieces.get(iLongPiece);
            final ReferenceMemorySource refSource = refSources.get(iLongPiece);
            final int strLength = strLengths.get(iLongPiece);
            //final int nStrs=longPiece.length()/strLength; //number of strs in the repeat
            final byte[] theseRefBases = refBases.get(iLongPiece).getBases();
            Integer indelIndexStart, indelIndexStop;
            if (longPiece.length() < 90) {
                indelIndexStart = repeatStart - 1;
                indelIndexStop = repeatStart + Math.min(longPiece.length(), 20);
            } else {
                indelIndexStart = repeatStart + 60;
                indelIndexStop = indelIndexStart + 1;
            }

            for (int indelIndex = indelIndexStart; indelIndex < indelIndexStop; indelIndex++) {
                for (int indelRepeats = Math.max(-(longPiece.length() - (indelIndex - repeatStart + 1)), -10) / strLength; indelRepeats < 10; indelRepeats++) {
                    if (indelRepeats == 0) {
                        continue;
                    }
                    final List<Allele> alleles = new ArrayList<Allele>();
                    if (indelRepeats < 0) { // deletion
                        byte[] basesRef = new byte[Math.abs(indelRepeats) * strLength + 1];
                        byte[] basesAlt = new byte[1];
                        System.arraycopy(theseRefBases, indelIndex, basesRef, 0, Math.abs(indelRepeats) * strLength + 1);
                        System.arraycopy(theseRefBases, indelIndex, basesAlt, 0, 1);
                        alleles.add(Allele.create(new String(basesRef, "UTF-8"), true));
                        alleles.add(Allele.create(new String(basesAlt, "UTF-8"), false));
                    } else {
                        byte[] basesRef = new byte[1];
                        byte[] basesAlt = new byte[Math.abs(indelRepeats) * strLength + 1];
                        byte[] basesRepeat = new byte[strLength];
                        System.arraycopy(theseRefBases, indelIndex, basesRef, 0, 1);
                        System.arraycopy(theseRefBases, indelIndex, basesAlt, 0, 1);
                        if (indelIndex < repeatStart + longPiece.length() - strLength) {
                            //look at next strLength bases to find str
                            System.arraycopy(theseRefBases, indelIndex + 1, basesRepeat, 0, strLength);
                        } else {
                            //look at previous strLength bases (including this one) to find str
                            System.arraycopy(theseRefBases, indelIndex - strLength + 1, basesRepeat, 0, strLength);
                        }
                        System.arraycopy(Utils.dupString(new String(basesRepeat, "UTF-8"), indelRepeats).getBytes(), 0, basesAlt, 1, Math.abs(indelRepeats) * strLength);
                        alleles.add(Allele.create(new String(basesRef, "UTF-8"), true));
                        alleles.add(Allele.create(new String(basesAlt, "UTF-8"), false));
                    }
                    final SimpleInterval interval = new SimpleInterval(artificialContig, indelIndex + 1, indelIndex + alleles.get(0).length());
                    ReferenceContext ref = new ReferenceContext(refSource, interval);

                    final VariantContext vc = new VariantContextBuilder("test", artificialContig, indelIndex + 1, indelIndex + alleles.get(0).length(), alleles).make();
                    tests.add(new Object[]{vc, ref, indelIndex != repeatStart - 1, repeatStart});
                }
            }
        }
        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "LeftAlignDataProvider")
    public void testLeftAlign(final VariantContext vc, final ReferenceContext ref, final boolean expectRealigned, final int expectedStart) {
        final LeftAlignAndTrimVariants leftAligner = new LeftAlignAndTrimVariants();
        final VariantContext realignedV = leftAligner.leftAlign(vc, ref);
        Assert.assertEquals(realignedV != vc, expectRealigned);
        Assert.assertEquals(realignedV.getStart(), expectedStart);
    }
}