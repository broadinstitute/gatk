package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.variant.variantcontext.*;
import org.broadinstitute.hellbender.tools.walkers.annotator.StrandArtifact;
import org.broadinstitute.hellbender.tools.walkers.annotator.StrandArtifact.StrandArtifactZ;
import org.broadinstitute.hellbender.utils.GATKProtectedVariantContextUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.genotyper.*;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.testng.Assert;
import org.testng.annotations.Test;


import java.io.IOException;
import java.util.*;

/**
 * Created by tsato on 4/19/17.
 */
public class StrandArtifactUnitTest {
    final int numChromosomes = 2; // create chromosome 0 and chromosome 1
    final int startingChromosome = 1;
    final int chromosomeSize = 1000;
    final SAMFileHeader SAM_HEADER = ArtificialReadUtils.createArtificialSamHeader(numChromosomes, startingChromosome, chromosomeSize);

    final int chromosomeIndex = 1;
    final long variantSite = 105;
    final int alignmentStart = 100;
    final int readLength = 9;

    final List<Allele> alleles = Arrays.asList(Allele.create((byte) 'C', true), Allele.create((byte) 'A', false));
    final AlleleList<Allele> alleleList = new IndexedAlleleList<>(alleles);
    final VariantContext vc = new VariantContextBuilder("source", Integer.toString(chromosomeIndex), variantSite, variantSite, alleles)
            .attribute(GATKVCFConstants.TUMOR_LOD_KEY, new double[]{ 6.4 }).make();

    final String sampleName = "SamThreetree";
    final SampleList sampleList = new IndexedSampleList(sampleName);

    final byte baseq = 30;
    final int mapq = 60;

    // all of the evidence for alt base comes from one strand
    @Test
    public void testSNP() throws IOException {
        // 25% ALT allele fraction at 100x coverage
        final int numAltReads = 25;
        final int numRefReads = 75;
        final int numReads = numAltReads + numRefReads;

        final Map<String, List<GATKRead>> readMap = new LinkedHashMap<>();
        final List<GATKRead> reads = new ArrayList<>(numReads);
        final byte[] quals = new byte[readLength];
        Arrays.fill(quals, baseq);

        // ALT reads are all reverse strand
        for (int i = 0; i < numAltReads; i++) {
            final GATKRead read = ArtificialReadUtils.createArtificialRead(SAM_HEADER, "Read" + i, chromosomeIndex, alignmentStart,
                    "CCCCACCCC".getBytes(), quals, "9M");
            read.setMappingQuality(mapq);
            read.setIsReverseStrand(true);
            reads.add(read);
        }

        // REF reads have good strand balance
        for (int i = numAltReads; i < numReads; i++) {
            final GATKRead read = ArtificialReadUtils.createArtificialRead(SAM_HEADER, "Read" + i, chromosomeIndex, alignmentStart,
                    "CCCCCCCCC".getBytes(), quals, "9M");
            read.setMappingQuality(mapq);
            read.setIsReverseStrand(i % 2 == 0);
            reads.add(read);
        }
        readMap.put(sampleName, reads);

        final ReadLikelihoods<Allele> likelihoods = new ReadLikelihoods<>(sampleList, alleleList, readMap);

        final int sampleIndex = 0;
        final LikelihoodMatrix<Allele> matrix = likelihoods.sampleMatrix(sampleIndex);

        final double logLikelihoodOfBestAllele = 10.0;
        final int refAlleleIndex = 0;
        final int altAlleleIndex = 1;

        // Log likelihoods are initialized to 0 for all alleles. For each read_i, set its likelihood for ALT allele to a positive value.
        // This makes read_i an ALT read
        for (int i = 0; i < numAltReads; i++) {
            matrix.set(altAlleleIndex, i, logLikelihoodOfBestAllele);
        }

        // Analogously, make read_j a REF read
        for (int j = numAltReads; j < numReads; j++) {
            matrix.set(refAlleleIndex, j, logLikelihoodOfBestAllele);
        }

        final StrandArtifact strandArtifactAnnotation = new StrandArtifact();
        final GenotypeBuilder genotypeBuilder = new GenotypeBuilder(sampleName);
        strandArtifactAnnotation.annotate(null, vc, genotypeBuilder.make(), genotypeBuilder, likelihoods);

        final Genotype genotype = genotypeBuilder.make();
        final double[] posteriorProbabilities =
                GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(genotype, StrandArtifact.POSTERIOR_PROBABILITIES_KEY, () -> null, -1);
        final double[] mapAlleleFractionEstimates =
                GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(genotype, StrandArtifact.MAP_ALLELE_FRACTIONS_KEY, () -> null, -1);

        final double epsilon = 1e-3;
        // Check that we correctly detect the artifact in reverse strands
        Assert.assertEquals(MathUtils.maxElementIndex(posteriorProbabilities), StrandArtifactZ.ART_REV.ordinal());
        Assert.assertEquals(posteriorProbabilities[StrandArtifactZ.ART_REV.ordinal()], 1.0, epsilon);

        // Check that, having taken strand artifact into consideration, we estimate that the true allele fraction is 0
        Assert.assertEquals(mapAlleleFractionEstimates[StrandArtifactZ.ART_REV.ordinal()], 0.0, epsilon);
    }

    // the underlying true allele fraction is non-zero
    @Test
    public void testHighAlleleFraction() throws IOException {
        // base allele fraction 0.2, but mix in 0.1 to one strand
        // test both high depth and low depth

        // 40% ALT allele fraction at 100x coverage; 20% is real, 20% is strand artifact
        // Since the underlying true allele fraction is high enough, the variant is real
        final int numAltReads = 40;
        final int numRefReads = 60;
        final int numReads = numAltReads + numRefReads;

        final Map<String, List<GATKRead>> readMap = new LinkedHashMap<>();
        final List<GATKRead> reads = new ArrayList<>(numReads);
        final byte[] quals = new byte[readLength];
        Arrays.fill(quals, baseq);

        // Half of ALT reads are all forward strand; these are the strand artifacts
        for (int i = 0; i < numAltReads/2; i++) {
            final GATKRead read = ArtificialReadUtils.createArtificialRead(SAM_HEADER, "Read" + i, chromosomeIndex, alignmentStart,
                    "CCCCACCCC".getBytes(), quals, "9M");
            read.setMappingQuality(mapq);
            read.setIsReverseStrand(false);
            reads.add(read);
        }

        // Half of ALT reads have good strand balance
        for (int i = numAltReads/2; i < numAltReads; i++) {
            final GATKRead read = ArtificialReadUtils.createArtificialRead(SAM_HEADER, "Read" + i, chromosomeIndex, alignmentStart,
                    "CCCCACCCC".getBytes(), quals, "9M");
            read.setMappingQuality(mapq);
            read.setIsReverseStrand(i % 2 == 0);
            reads.add(read);
        }

        // All of REF reads have good strand balance
        for (int i = numAltReads; i < numReads; i++) {
            final GATKRead read = ArtificialReadUtils.createArtificialRead(SAM_HEADER, "Read" + i, chromosomeIndex, alignmentStart,
                    "CCCCCCCCC".getBytes(), quals, "9M");
            read.setMappingQuality(mapq);
            read.setIsReverseStrand(i % 2 == 0);
            reads.add(read);
        }

        readMap.put(sampleName, reads);

        final ReadLikelihoods<Allele> likelihoods = new ReadLikelihoods<>(sampleList, alleleList, readMap);

        final int sampleIndex = 0;
        final LikelihoodMatrix<Allele> matrix = likelihoods.sampleMatrix(sampleIndex);

        final double logLikelihoodOfBestAllele = 10.0;
        final int refAlleleIndex = 0;
        final int altAlleleIndex = 1;

        // Log likelihoods are initialized to 0 for all alleles. For each read_i, set its likelihood for ALT allele to a positive value.
        // This makes read_i an ALT read
        for (int i = 0; i < numAltReads; i++) {
            matrix.set(altAlleleIndex, i, logLikelihoodOfBestAllele);
        }

        // Analogously, make read_j a REF read
        for (int j = numAltReads; j < numReads; j++) {
            matrix.set(refAlleleIndex, j, logLikelihoodOfBestAllele);
        }

        final StrandArtifact strandArtifactAnnotation = new StrandArtifact();
        final GenotypeBuilder genotypeBuilder = new GenotypeBuilder(sampleName);
        strandArtifactAnnotation.annotate(null, vc, genotypeBuilder.make(), genotypeBuilder, likelihoods);

        final Genotype genotype = genotypeBuilder.make();
        final double[] posteriorProbabilities =
                GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(genotype, StrandArtifact.POSTERIOR_PROBABILITIES_KEY, () -> null, -1);
        final double[] mapAlleleFractionEstimates =
                GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(genotype, StrandArtifact.MAP_ALLELE_FRACTIONS_KEY, () -> null, -1);

        final double epsilon = 0.1;
        // Check that, having taken strand artifact into consideration, we estimate that the true allele fraction is 0.2
        Assert.assertEquals(mapAlleleFractionEstimates[StrandArtifactZ.ART_FWD.ordinal()], 0.2, epsilon);

    }

}