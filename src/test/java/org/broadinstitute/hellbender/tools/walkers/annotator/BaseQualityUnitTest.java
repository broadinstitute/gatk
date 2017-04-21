package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.variant.variantcontext.*;
import org.broadinstitute.hellbender.utils.genotyper.*;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.*;

import static org.testng.Assert.*;

/**
 * Created by davidben on 3/24/17.
 */
public class BaseQualityUnitTest {
    final static private String SAMPLE = "SAMPLE";
    final static private Genotype DUMMY_GENOTYPE = new GenotypeBuilder(SAMPLE).make();

    @Test
    public void test() {
        final SAMFileHeader SAM_HEADER = ArtificialReadUtils.createArtificialSamHeader(10, 0, 1000);

        final List<Allele> alleles = Arrays.asList(Allele.create((byte) 'A', true), Allele.create((byte) 'C', false));
        final AlleleList<Allele> alleleList = new IndexedAlleleList<>(alleles);

        // variant is a SNP at position 20
        final int chromosomeIndex = 5;
        final int variantSite = 20;
        final VariantContext vc = new VariantContextBuilder("source", Integer.toString(chromosomeIndex), variantSite, variantSite, alleles).make();
        final org.broadinstitute.hellbender.utils.genotyper.SampleList sampleList = new IndexedSampleList("SAMPLE");

        //7 reads of length = 3 -- the middle base is at the variant position
        final Map<String, List<GATKRead>> readMap = new LinkedHashMap<>();
        final List<GATKRead> reads = new ArrayList<>();
        final byte[] quals = new byte[]{30, 30, 30, 30, 25, 25, 25};
        for (int r = 0; r < quals.length; r++) {
            final GATKRead read = ArtificialReadUtils.createArtificialRead(SAM_HEADER,
                    "RRR00" + r, chromosomeIndex, 19, "ACG".getBytes(), new byte[]{4, quals[r], 55}, "3M");
            read.setMappingQuality(60);
            reads.add(read);
        }
        readMap.put("SAMPLE", reads);

        final ReadLikelihoods<Allele> likelihoods = new ReadLikelihoods<>(sampleList, alleleList, readMap);

        //we will make the first four reads ref (median position = 2) and the last three alt (median position 10, hence
        // median distance from end = 1)
        final LikelihoodMatrix<Allele> matrix = likelihoods.sampleMatrix(0);

        // log likelihoods are initialized to 0, so we can "turn on" a read for a particular allele by setting the
        // (allele, read) entry to 10
        matrix.set(0, 0, 10);
        matrix.set(0, 1, 10);
        matrix.set(0, 2, 10);
        matrix.set(0, 3, 10);
        matrix.set(1, 4, 10);
        matrix.set(1, 5, 10);
        matrix.set(1, 6, 10);

        final BaseQuality bq = new BaseQuality();
        final GenotypeBuilder gb = new GenotypeBuilder(DUMMY_GENOTYPE);

        bq.annotate(null, vc, DUMMY_GENOTYPE, gb, likelihoods);
        final Genotype g = gb.make();

        final int[] medianRefAndAltQuals = (int[]) g.getExtendedAttribute(BaseQuality.KEY);

        Assert.assertEquals( medianRefAndAltQuals[0], 30);
        Assert.assertEquals(medianRefAndAltQuals[1], 25);
    }
}