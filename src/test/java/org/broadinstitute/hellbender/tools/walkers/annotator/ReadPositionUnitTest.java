package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.variant.variantcontext.*;
import org.broadinstitute.hellbender.utils.genotyper.*;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.*;

/**
 * Created by davidben on 3/23/17.
 */
public class ReadPositionUnitTest extends GATKBaseTest {
    final static private String SAMPLE = "SAMPLE";
    final static private Genotype DUMMY_GENOTYPE = new GenotypeBuilder(SAMPLE).make();

    @Test
    public void test() {
        final SAMFileHeader SAM_HEADER = ArtificialReadUtils.createArtificialSamHeader(10, 0, 1000);

        final List<Allele> alleles = Arrays.asList(Allele.create((byte) 'A', true), Allele.create((byte) 'C', false));
        final AlleleList<Allele> alleleList = new IndexedAlleleList<>(alleles);

        final int chromosomeIndex = 5;

        // variant is a SNP at position 20
        final int variantSite = 20;
        final VariantContext vc = new VariantContextBuilder("source", Integer.toString(chromosomeIndex), variantSite, variantSite, alleles).make();
        final SampleList sampleList = new IndexedSampleList("SAMPLE");

        //7 length-12 reads
        final Map<String,List<GATKRead>> readMap = new LinkedHashMap<>();
        final List<GATKRead> reads = new ArrayList<>();
        final int[] positionsOfSiteWithinReads = new int[] {1, 2, 2, 3, 10, 10, 11};
        final int[] alignmentStarts = Arrays.stream(positionsOfSiteWithinReads).map(n -> variantSite - n).toArray();
        for (int r = 0; r < positionsOfSiteWithinReads.length; r++) {
            final GATKRead read = ArtificialReadUtils.createArtificialRead(SAM_HEADER,
                    "RRR00" + r, chromosomeIndex, alignmentStarts[r], "ACGTACGTACGT".getBytes(), new byte[]{30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30}, "12M");
            read.setMappingQuality(60);
            reads.add(read);
        }
        readMap.put("SAMPLE",reads);

        final AlleleLikelihoods<GATKRead, Allele> likelihoods = new AlleleLikelihoods<>(sampleList, alleleList, readMap);

        //we will make the first four reads ref (median position = 2) and the last three alt (median position 10, hence
        // median distance from end = 1)
        final LikelihoodMatrix<GATKRead, Allele> matrix = likelihoods.sampleMatrix(0);

        // log likelihoods are initialized to 0, so we can "turn on" a read for a particular allele by setting the
        // (allele, read) entry to 10
        matrix.set(0,0, 10);
        matrix.set(0,1, 10);
        matrix.set(0,2, 10);
        matrix.set(0,3, 10);
        matrix.set(1,4, 10);
        matrix.set(1,5, 10);
        matrix.set(1,6, 10);

        final ReadPosition rp = new ReadPosition();
        final GenotypeBuilder gb = new GenotypeBuilder(DUMMY_GENOTYPE);

        final Map<String, Object> annotation = rp.annotate(null, vc, likelihoods);

        final int[] medianAltPositions = (int[]) annotation.get(GATKVCFConstants.MEDIAN_READ_POSITON_KEY);

        Assert.assertEquals(medianAltPositions[0], 1);
    }

}