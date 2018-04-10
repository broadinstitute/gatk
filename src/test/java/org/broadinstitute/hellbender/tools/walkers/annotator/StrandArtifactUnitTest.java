package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.variant.variantcontext.*;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.utils.GATKProtectedVariantContextUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.*;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

import static org.broadinstitute.hellbender.tools.walkers.annotator.StrandArtifact.ReadStrand.FWD;
import static org.broadinstitute.hellbender.tools.walkers.annotator.StrandArtifact.ReadStrand.REV;

/**
 * Created by tsato on 4/6/18.
 */
public class StrandArtifactUnitTest {
    final int numChromosomes = 2; // create chromosome 0 and chromosome 1
    final int startingChromosome = 1;
    final int chromosomeSize = 1000;
    final SAMFileHeader SAM_HEADER = ArtificialReadUtils.createArtificialSamHeader(numChromosomes, startingChromosome, chromosomeSize);

    final int chromosomeIndex = 1;
    final long variantSite = 105;

    final int readLength = 9;

    private static final String NEITHER = "neither";

    final List<Allele> alleles = Arrays.asList(Allele.create((byte) 'C', true), Allele.create((byte) 'A', false));
    final AlleleList<Allele> alleleList = new IndexedAlleleList<>(alleles);
    final VariantContext vc = new VariantContextBuilder("source", Integer.toString(chromosomeIndex), variantSite, variantSite, alleles)
            .attribute(GATKVCFConstants.TUMOR_LOD_KEY, new double[]{ 6.4 }).make();

    final String sampleName = "SamThreetree";
    final org.broadinstitute.hellbender.utils.genotyper.SampleList sampleList = new IndexedSampleList(sampleName);


    @DataProvider
    public Object[][] failingCases(){
        return new Object[][]{
                {19, 61, 10, 2, NEITHER},
        };
    }

    @DataProvider
    public Object[][] snpTestData(){
        return new Object[][]{
                // REF+, REF-, ALT+, ALT-
                {25, 25, 10, 0, "FWD"}, // Classic case of FWD artifact
                {100, 100, 10, 0, "FWD"}, // Low allele fraction
                {0, 50, 0, 5, NEITHER}, // If the ref is biased then alt can be biased too
                {50, 0, 5, 0, NEITHER}, // Test the other way
                {25, 25, 0, 5, "REV"}, // This is border line. Currently we get LOD of 1.2 or so, which is reasonable
                {25, 25, 0, 2, NEITHER}, // Here, the likelihood for artifact should win, but not enough to be filtered
                {40, 40, 0, 7, "REV"}, // The old filter gave this case the posterior artifact probability of 0.3
                {0, 0, 0, 7, "REV"}, // HOM ALT site, unbalanced
                {0, 0, 4, 6, NEITHER}, // HOM ALT, balanced
                {1, 0, 4, 6, NEITHER}, // HOM ALT with a single REF read (possibly a contaminant)
                {15, 0, 5, 5, NEITHER}, // REF is biased, ALT is balanced
                {15, 0, 0, 5, "REV"}, // Opposite bias between REF and ALT
                {0, 15, 5, 0, "FWD"}, // Other way
                {30, 15, 15, 30, NEITHER}, // REF and ALT have different distributions
                {50, 20, 5, 4, NEITHER},
                {25, 25, 1, 0, NEITHER}, // Low alt count
                {10, 40, 1, 0, NEITHER}, // Low alt count with opposing imbalance in ref
                {0, 10, 1, 0, "FWD"}, // Ref and Alt biased in opposite directions - should it be filtered?
                {53, 28, 14, 2, NEITHER}, // Gray area here...is this strand artifact?
                {19, 61, 10, 2, NEITHER} // Another gray area
        };
    }

    // all of the evidence for alt base comes from one strand
    @Test(dataProvider = "failingCases")
    public void testSNP(final int numRefFwd, final int numRefRev,
                        final int numAltFwd, final int numAltRev,
                        final String expectedArtifact){

        final Map<String, List<GATKRead>> readMap = new LinkedHashMap<>();
        final List<GATKRead> reads = new ArrayList<>();

        final String refReadBases = "CCCCCCCCC";
        final String altReadBases = "CCCCACCCC";

        addReads(numRefFwd, reads, refReadBases, FWD);
        addReads(numRefRev, reads, refReadBases, REV);
        addReads(numAltFwd, reads, altReadBases, FWD);
        addReads(numAltRev, reads, altReadBases, REV);

        readMap.put(sampleName, reads);

        final ReadLikelihoods<Allele> likelihoods = new ReadLikelihoods<>(sampleList, alleleList, readMap);

        final int sampleIndex = 0;
        final LikelihoodMatrix<Allele> matrix = likelihoods.sampleMatrix(sampleIndex);

        final int refIndex = 0;
        final int altIndex = 1;

        // The order here must reflect the order in which we added corresponding reads
        List<Pair<Integer, Integer>> counts = Arrays.asList(
                new ImmutablePair<>(numRefFwd, refIndex),
                new ImmutablePair<>(numRefRev, refIndex),
                new ImmutablePair<>(numAltFwd, altIndex),
                new ImmutablePair<>(numAltRev, altIndex));
        setLikelihoods(matrix, counts);

        final StrandArtifact annotation = new StrandArtifact();
        final GenotypeBuilder genotypeBuilder = new GenotypeBuilder(sampleName);
        annotation.annotate(null, vc, genotypeBuilder.make(), genotypeBuilder, likelihoods);

        final Genotype genotype = genotypeBuilder.make();
        final String artifact = GATKProtectedVariantContextUtils.getAttributeAsString(genotype,
                StrandArtifact.STRAND_ARTIFACT_DIRECTION_KEY, NEITHER);
        // This is more of a sanity check - make sure the artifact probability exists and is greater than 0
        final double artifactProbability = GATKProtectedVariantContextUtils.getAttributeAsDouble(genotype,
                StrandArtifact.STRAND_ARTIFACT_PROBABILITY_KEY, -1.0);
        if (! expectedArtifact.equals(NEITHER)){
            Assert.assertTrue(artifactProbability > 0);
        }

        Assert.assertEquals(artifact, expectedArtifact);
    }

    private void setLikelihoods(final LikelihoodMatrix<Allele> matrix, final List<Pair<Integer, Integer>> alleleCounts) {
        final double logLikelihoodOfBestAllele = 10.0;
        int i = 0;
        for (Pair<Integer, Integer> counts : alleleCounts){
            final int numReads = counts.getLeft();
            final int alleleIndex = counts.getRight();
            for (int j = 0; j < numReads; j++){
                matrix.set(alleleIndex, i, logLikelihoodOfBestAllele);
                i++;
            }
        }

    }

    private void addReads(final int numReads, final List<GATKRead> reads, final String readBases, final StrandArtifact.ReadStrand strand){
        Utils.nonNull(readBases);
        Utils.nonNull(reads);

        final int alignmentStart = 100;
        final byte[] quals = new byte[readLength];
        final byte baseq = 30;
        final int mapq = 60;
        Arrays.fill(quals, baseq);

        // Assume no indels
        final String cigar = readLength + "M";

        // ALT reads are all reverse strand
        for (int i = 0; i < numReads; i++) {
            // Duplicated read names across multiple calls of this method should not be an issue within this test
            final GATKRead read = ArtificialReadUtils.createArtificialRead(SAM_HEADER, "Read" + i, chromosomeIndex,
                    alignmentStart, readBases.getBytes(), quals, cigar);
            read.setMappingQuality(mapq);
            read.setIsReverseStrand(strand == REV);
            reads.add(read);
        }
    }

}