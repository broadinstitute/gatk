package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.genotyper.*;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.*;

import static org.testng.Assert.*;

/**
 * Created by davidben on 3/23/17.
 */
public class ReadPositionUnitTest extends BaseTest {
    @Test
    public void test() {
        final SAMFileHeader SAM_HEADER = ArtificialReadUtils.createArtificialSamHeader(10, 0, 1000);

        final List<Allele> alleles = Arrays.asList(Allele.create((byte) 'A', true), Allele.create((byte) 'C', false));
        final AlleleList<Allele> alleleList = new IndexedAlleleList<>(alleles);

        // variant is a SNP at position 10
        final int variantSite = 20;
        final VariantContext vc = new VariantContextBuilder("source", "5", variantSite, variantSite, alleles).make();
        final org.broadinstitute.hellbender.utils.genotyper.SampleList sampleList = new IndexedSampleList("SAMPLE");

        //7 length-12 reads
        final Map<String,List<GATKRead>> readMap = new LinkedHashMap<>();
        final List<GATKRead> reads = new ArrayList<>();
        final int[] positionsOfSiteWithinReads = new int[] {1, 2, 2, 3, 10, 10, 11};
        final int[] alignmentStarts = Arrays.stream(positionsOfSiteWithinReads).map(n -> variantSite - n).toArray();
        for (int r = 0; r < 7; r++) {
            final GATKRead read = ArtificialReadUtils.createArtificialRead(SAM_HEADER,
                    "RRR00" + r, 0, alignmentStarts[r], "ACGTACGTACGT".getBytes(), new byte[]{30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30}, "12M");
            read.setMappingQuality(60);
            reads.add(read);
        }
        readMap.put("SAMPLE",reads);

        final ReadLikelihoods<Allele> likelihoods = new ReadLikelihoods<>(sampleList, alleleList, readMap);

        //we will make the first four reads ref (median position = 2) and the last three alt (median position 10, hence
        // median distance from end = 1)
        final LikelihoodMatrix<Allele> matrix = likelihoods.sampleMatrix(0);

        // log likelihoods are initialized to 0, so we can "turn on" a read for a particular allele by setting the
        // (allele, read) entry to 10
        matrix.set(0,0, 10);
        matrix.set(0,1, 10);
        matrix.set(0,2, 10);
        matrix.set(0,3, 10);
        matrix.set(1,4, 10);
        matrix.set(1,5, 10);
        matrix.set(1,6, 10);

        final Map<String, Object> annotations = new ReadPosition().annotate(null, vc, likelihoods);

        Assert.assertEquals((int) annotations.get(ReadPosition.REFERENCE_MEDIAN_POSITION_KEY), 2);
        Assert.assertEquals((int) annotations.get(ReadPosition.ALT_MEDIAN_POSITION_KEY), 1);
    }

}