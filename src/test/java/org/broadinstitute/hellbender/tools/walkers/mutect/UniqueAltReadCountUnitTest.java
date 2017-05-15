package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.variant.variantcontext.*;
import org.broadinstitute.hellbender.tools.walkers.annotator.UniqueAltReadCount;
import org.broadinstitute.hellbender.utils.GATKProtectedVariantContextUtils;
import org.broadinstitute.hellbender.utils.genotyper.*;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.*;

public class UniqueAltReadCountUnitTest {
    final String sampleName = "Mark";
    final int chromosomeIndex = 1;
    final long variantSite = 105;
    final int numAltReads = 12;

    final SampleList sampleList = new IndexedSampleList(sampleName);

    final List<Allele> alleles = Arrays.asList(Allele.create((byte) 'C', true), Allele.create((byte) 'A', false));
    final AlleleList<Allele> alleleList = new IndexedAlleleList<>(alleles);
    final VariantContext vc = new VariantContextBuilder("source", Integer.toString(chromosomeIndex), variantSite, variantSite, alleles).make();

    @Test
    public void testSingleDuplicate() throws IOException {
        final ReadLikelihoods<Allele> likelihoods = createTestLikelihoods(Optional.empty());
        final UniqueAltReadCount uniqueAltReadCountAnnotation = new UniqueAltReadCount();
        final GenotypeBuilder genotypeBuilder = new GenotypeBuilder(sampleName);
        uniqueAltReadCountAnnotation.annotate(null, vc, genotypeBuilder.make(), genotypeBuilder, likelihoods);

        final Genotype genotype = genotypeBuilder.make();

        final int uniqueReadSetCount = GATKProtectedVariantContextUtils.getAttributeAsInt(genotype, UniqueAltReadCount.UNIQUE_ALT_READ_SET_COUNT_KEY, -1);
        Assert.assertEquals(uniqueReadSetCount, 1);
    }

    @Test
    public void testMultipleDuplicateSets() throws IOException {
        final UniqueAltReadCount duplicateReadCountsAnnotation = new UniqueAltReadCount();

        // should get three unique sets of ALT reads
        final int numUniqueStarts1 = 3;
        final ReadLikelihoods<Allele> likelihoods1 = createTestLikelihoods(Optional.of(numUniqueStarts1));
        final GenotypeBuilder genotypeBuilder1 = new GenotypeBuilder(sampleName);
        duplicateReadCountsAnnotation.annotate(null, vc, genotypeBuilder1.make(), genotypeBuilder1, likelihoods1);

        final Genotype genotype1 = genotypeBuilder1.make();

        final int uniqueReadSetCount1 = GATKProtectedVariantContextUtils.getAttributeAsInt(genotype1, UniqueAltReadCount.UNIQUE_ALT_READ_SET_COUNT_KEY, -1);

        Assert.assertEquals(uniqueReadSetCount1, numUniqueStarts1);

        // here ALT reads are all distinct
        final int numUniqueStarts2 = numAltReads;
        final ReadLikelihoods<Allele> likelihoods2 = createTestLikelihoods(Optional.of(numUniqueStarts2));
        final GenotypeBuilder genotypeBuilder2 = new GenotypeBuilder(sampleName);
        duplicateReadCountsAnnotation.annotate(null, vc, genotypeBuilder2.make(), genotypeBuilder2, likelihoods2);

        final Genotype genotype2 = genotypeBuilder2.make();

        final int uniqueReadSetCount2 = GATKProtectedVariantContextUtils.getAttributeAsInt(genotype2, UniqueAltReadCount.UNIQUE_ALT_READ_SET_COUNT_KEY, -1);

        Assert.assertEquals(uniqueReadSetCount2, numUniqueStarts2);
    }

    private ReadLikelihoods<Allele> createTestLikelihoods(final Optional<Integer> shiftModulus) {
        final int numChromosomes = 2;
        final int startingChromosome = 1;
        final int chromosomeSize = 1000;
        final SAMFileHeader SAM_HEADER = ArtificialReadUtils.createArtificialSamHeader(numChromosomes, startingChromosome, chromosomeSize);

        final int alignmentStart = 100;
        final int readLength = 11;
        final byte baseq = 30;

        final int numRefReads = 10;
        final int numReads = numAltReads + numRefReads;

        final Map<String, List<GATKRead>> readMap = new LinkedHashMap<>();
        final List<GATKRead> reads = new ArrayList<>(numReads);
        final byte[] quals = new byte[readLength];
        Arrays.fill(quals, baseq);

        // add alt reads, with the start position shifted by i mod shiftModulus
        for (int i = 0; i < numAltReads; i++) {
            final int startPosition = shiftModulus.isPresent() ? alignmentStart + (i % shiftModulus.get()) : alignmentStart;
            final GATKRead read = ArtificialReadUtils.createArtificialRead(SAM_HEADER, "Read" + i, chromosomeIndex, startPosition,
                    "CCCCCACCCCC".getBytes(), quals, "11M");
            reads.add(read);
        }

        // add ref reads
        for (int j = numAltReads; j < numReads; j++) {
            final GATKRead read = ArtificialReadUtils.createArtificialRead(SAM_HEADER, "Read" + j, chromosomeIndex, alignmentStart,
                    "CCCCCCCCCCC".getBytes(), quals, "11M");
            reads.add(read);
        }

        readMap.put(sampleName, reads);

        final ReadLikelihoods<Allele> likelihoods = new ReadLikelihoods<>(sampleList, alleleList, readMap);
        final int sampleIndex = 0;
        final LikelihoodMatrix<Allele> matrix = likelihoods.sampleMatrix(sampleIndex);

        final double logLikelihoodOfBestAllele = 10.0;
        final int refAlleleIndex = 0;
        final int altAlleleIndex = 1;

        // Log likelihoods are initialized to 0 for all alleles. For each read_i, set its log likelihood for ALT allele to a positive value.
        // This makes read_i an ALT read
        for (int i = 0; i < numAltReads; i++) {
            matrix.set(altAlleleIndex, i, logLikelihoodOfBestAllele);
        }

        // Analogously, make read_j a REF read
        for (int j = numAltReads; j < numReads; j++) {
            matrix.set(refAlleleIndex, j, logLikelihoodOfBestAllele);
        }

        return likelihoods;
    }
}