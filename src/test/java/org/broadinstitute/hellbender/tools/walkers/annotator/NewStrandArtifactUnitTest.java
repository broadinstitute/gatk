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
import org.broadinstitute.hellbender.utils.test.ArtificialAnnotationUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.*;

import static org.broadinstitute.hellbender.tools.walkers.annotator.NewStrandArtifact.ReadStrand.FWD;
import static org.broadinstitute.hellbender.tools.walkers.annotator.NewStrandArtifact.ReadStrand.REV;

/**
 * Created by tsato on 4/6/18.
 */
public class NewStrandArtifactUnitTest {
    final int numChromosomes = 2; // create chromosome 0 and chromosome 1
    final int startingChromosome = 1;
    final int chromosomeSize = 1000;
    final SAMFileHeader SAM_HEADER = ArtificialReadUtils.createArtificialSamHeader(numChromosomes, startingChromosome, chromosomeSize);

    final int chromosomeIndex = 1;
    final long variantSite = 105;

    final int readLength = 9;

    final List<Allele> alleles = Arrays.asList(Allele.create((byte) 'C', true), Allele.create((byte) 'A', false));
    final AlleleList<Allele> alleleList = new IndexedAlleleList<>(alleles);
    final VariantContext vc = new VariantContextBuilder("source", Integer.toString(chromosomeIndex), variantSite, variantSite, alleles)
            .attribute(GATKVCFConstants.TUMOR_LOD_KEY, new double[]{ 6.4 }).make();

    final String sampleName = "SamThreetree";
    final org.broadinstitute.hellbender.utils.genotyper.SampleList sampleList = new IndexedSampleList(sampleName);


    @DataProvider
    public Object[][] snpTestData(){
        return new Object[][]{
                // REF+, REF-, ALT+, ALT-
                {25, 25, 10, 0, "FWD"},
                {100, 100, 10, 0, "FWD"}, // Low allele fraction
                {0, 50, 0, 5, "neither"}, // If the ref is biased then alt can be biased too
                {50, 0, 5, 0, "neither"}, // Test the other way
                {25, 25, 0, 5, "neither"}, // This is border line. Currently we get LOD of 1.2 or so, which is reasonable
                {25, 25, 0, 2, "neither"}, // Here, the likelihood for artifact should win, but not enough to be filtered
                {40, 40, 0, 7, "REV"}, // The old filter gave this case the posterior artifact probability of 0.3
        };
    }

    // all of the evidence for alt base comes from one strand
    @Test(dataProvider = "snpTestData")
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

        final NewStrandArtifact annotation = new NewStrandArtifact();
        final GenotypeBuilder genotypeBuilder = new GenotypeBuilder(sampleName);
        annotation.annotate(null, vc, genotypeBuilder.make(), genotypeBuilder, likelihoods);

        final Genotype genotype = genotypeBuilder.make();
        final String artifact = GATKProtectedVariantContextUtils.getAttributeAsString(genotype,
                NewStrandArtifact.STRAND_ARTIFACT_DIRECTION_KEY, "neither");
        // This is more of a sanity check - make sure the lod exists and exceeds threshold
        final double lod = GATKProtectedVariantContextUtils.getAttributeAsDouble(genotype,
                NewStrandArtifact.STRAND_ARTIFACT_LOG10_ODDS_KEY, -1.0);

        Assert.assertEquals(artifact, expectedArtifact);
        final double lodThresholdForBorderLineCases = 1.0;
        if (! expectedArtifact.equals("neither")){
            Assert.assertTrue(lod > lodThresholdForBorderLineCases);
        }
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

    private void addReads(final int numReads, final List<GATKRead> reads, final String readBases, final NewStrandArtifact.ReadStrand strand){
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