package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.variant.variantcontext.*;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.genotyper.*;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.testng.Assert;
import org.testng.annotations.Test;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;

import java.util.*;

public class StrandBiasBySampleUnitTest extends GATKBaseTest {

    @Test
    public void testGetContingencyArray(){
        final int[][] t = new int[2][2];
        t[0][0] = 1; t[0][1] = 2; t[1][0] = 3; t[1][1] = 4;
        final List<Integer> tList = StrandBiasBySample.getContingencyArray(t);
        final List<Integer> truthList = new ArrayList<>();
        truthList.add(1); truthList.add(2); truthList.add(3); truthList.add(4);
        Assert.assertEquals(tList, truthList);
    }

    @Test
    public void testOverlappingReads(){
        final SAMFileHeader SAM_HEADER = ArtificialReadUtils.createArtificialSamHeader(10, 0, 1000);

        final List<Allele> alleles = Arrays.asList(Allele.create((byte) 'A', true), Allele.create((byte) 'C', false));
        final AlleleList<Allele> alleleList = new IndexedAlleleList<>(alleles);

        // variant is a SNP at position 20
        final int chromosomeIndex = 5;
        final int variantSite = 20;
        final VariantContext vc = new VariantContextBuilder("source", Integer.toString(chromosomeIndex), variantSite, variantSite, alleles).make();
        final SampleList sampleList = new IndexedSampleList("SAMPLE");

        final Map<String, List<GATKRead>> readMap = new LinkedHashMap<>();
        final List<GATKRead> reads = new ArrayList<>();
        final byte[] quals = new byte[]{30, 30, 30, 30, 25, 25, 25};
        for (int r = 0; r < quals.length; r++) {
            final GATKRead read = ArtificialReadUtils.createArtificialRead(SAM_HEADER,
                    "Read" + r, chromosomeIndex, 19, "ACG".getBytes(), new byte[]{4, quals[r], 55}, "3M");
            final GATKRead mate = ArtificialReadUtils.createArtificialRead(SAM_HEADER,
                    "Read" + r, chromosomeIndex, 19, "ACG".getBytes(), new byte[]{4, quals[r], 55}, "3M");
            read.setMappingQuality(60);
            read.setIsReverseStrand(false);
            mate.setIsReverseStrand(true);
            read.setIsFirstOfPair();
            mate.setIsSecondOfPair();
            reads.add(read);
            reads.add(mate);
        }

        readMap.put("SAMPLE", reads);

        final ReadLikelihoods<Allele> likelihoods = new ReadLikelihoods<>(sampleList, alleleList, readMap);

        //we will make the first four reads ref (median position = 2) and the last three alt (median position 10, hence
        // median distance from end = 1)
        final LikelihoodMatrix<Allele> matrix = likelihoods.sampleMatrix(0);

        // log likelihoods are initialized to 0, so we can "turn on" a read for a particular allele by setting the
        // (allele, read) entry to 10
        for (int r = 0; r < quals.length; r++) {
            if (r % 2 == 0){
                matrix.set(0, 2*r, 10);
                matrix.set(0, 2*r + 1, 10);
            } else {
                matrix.set(1, 2*r, 10);
                matrix.set(1, 2*r + 1, 10);
            }
        }

        final StrandBiasBySample strandBiasBySample = new StrandBiasBySample();
        final GenotypeBuilder gb = new GenotypeBuilder("SAMPLE").alleles(alleles);

        strandBiasBySample.annotate(null, vc, gb.make(), gb, likelihoods);
        final Genotype g = gb.make();

        final List<Integer> contingencyTable = (List<Integer>)g.getExtendedAttribute(GATKVCFConstants.STRAND_BIAS_BY_SAMPLE_KEY);

        Assert.assertEquals(contingencyTable.get(0).intValue(), 4);
        Assert.assertEquals(contingencyTable.get(1).intValue(), 4);
        Assert.assertEquals(contingencyTable.get(2).intValue(), 3);
        Assert.assertEquals(contingencyTable.get(3).intValue(), 3);

    }
}
