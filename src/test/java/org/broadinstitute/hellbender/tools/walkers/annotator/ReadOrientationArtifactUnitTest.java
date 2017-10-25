package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.variant.variantcontext.*;
import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.RandomGeneratorFactory;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceFileSource;
import org.broadinstitute.hellbender.utils.GATKProtectedVariantContextUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.genotyper.*;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.*;

import static org.testng.Assert.*;

/**
 * Created by tsato on 10/23/17.
 */
public class ReadOrientationArtifactUnitTest extends BaseTest {
    final int chromosomeIndex = 1;
    final int variantSite = 100_000;
    final String sampleName = "someTumor";

    @DataProvider(name = "readOrientationArtifact")
    public Object[][] readOrientationArtifactData() {
        return new Object[][]{
                {0.3, 0.5, 0.0, 0.0}, // allele fraction 0.3 variant with no bias
                {0.3, 0.99, 1.0, 0.0}, // allele fraction 0.3 variant with heavy F1R2 bias
                {0.3, 0.94, 1.0, 0.0}, // allele fraction 0.2 variant with moderate F1R2 bias
                {0.2, 1e-2, 0.0, 1.0}}; // allele fraction 0.2 variant with a heavy F2R1 bias
    }

    @Test(dataProvider="readOrientationArtifact")
    public void test(final double alleleFraction, final double altF1R2Fraction,
                     final double expectedF1R2Prob, final double expectedF2R1Prob){
        final List<Allele> alleles = Arrays.asList(Allele.create((byte) 'C', true), Allele.create((byte) 'A', false));
        final VariantContext vc = new VariantContextBuilder("source", Integer.toString(chromosomeIndex), variantSite, variantSite, alleles)
                .make();

        final ReferenceContext ref = new ReferenceContext(new ReferenceFileSource(new File(hg19_chr1_1M_Reference)),
                new SimpleInterval("1", variantSite, variantSite));

        ReadOrientationArtifact annotation = new ReadOrientationArtifact();

        final GenotypeBuilder genotypeBuilder = new GenotypeBuilder(sampleName);
        final ReadLikelihoods<Allele> readLikelihoods = createReadLikelihoods(alleleFraction, altF1R2Fraction,
                alleles, genotypeBuilder);
        genotypeBuilder.alleles(alleles);

        File hyperparmaeterTable = new File(publicTestDir, "HSCX1127T-hyperparameters.tsv");
        annotation.setHyperparameterTable(hyperparmaeterTable);

        annotation.annotate(ref, vc, genotypeBuilder.make(), genotypeBuilder, readLikelihoods);

        final Genotype genotype = genotypeBuilder.make();
        final double[] posteriorProbabilities =
                GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(genotype, GATKVCFConstants.READ_ORIENTATION_POSTERIOR_KEY, () -> null, -1);
        final double posteriorProbOfF1R2Artifact = posteriorProbabilities[0];
        final double posteriorProbOfF2R1Artifact = posteriorProbabilities[1];

        final double epsilon = 1e-3;
        Assert.assertEquals(posteriorProbOfF1R2Artifact, expectedF1R2Prob, epsilon);
        Assert.assertEquals(posteriorProbOfF2R1Artifact, expectedF2R1Prob, epsilon);
    }

    /** Update the genotype builder by side effect, too **/
    private ReadLikelihoods<Allele> createReadLikelihoods(final double alleleFraction, final double altF1R2Fraction,
                                                          final List<Allele> alleles, final GenotypeBuilder genotypeBuilder){
        final RandomGenerator randomGenerator = RandomGeneratorFactory.createRandomGenerator(new Random(111));

        final int depth = 100;
        final BinomialDistribution altReadDistribution = new BinomialDistribution(randomGenerator, depth, alleleFraction);
        final int numAltReads = altReadDistribution.sample();

        final BinomialDistribution altF1R2Distribution = new BinomialDistribution(randomGenerator, numAltReads, altF1R2Fraction);
        final int numAltF1R2Reads = altF1R2Distribution.sample();
        final BinomialDistribution refF1R2Distribution = new BinomialDistribution(randomGenerator, numAltReads, 0.5);
        final int numRefF1R2Reads = refF1R2Distribution.sample();

        final byte baseq = 20;
        final int mapq = 30;

        // create a sam header
        final int numChromosomes = 1;
        final int startingChromosome = 0;
        final int chromosomeSize = 1_000_000;

        final SAMFileHeader samHeader = ArtificialReadUtils.createArtificialSamHeader(numChromosomes, startingChromosome, chromosomeSize);

        final SampleList sampleList = new IndexedSampleList(genotypeBuilder.make().getSampleName());
        final AlleleList<Allele> alleleList = new IndexedAlleleList<>(alleles);

        final int chromosomeIndex = 0;
        final int alignmentStart = 100;
        final int readLength = 9;

        final Map<String, List<GATKRead>> readMap = new LinkedHashMap<>();
        final List<GATKRead> reads = new ArrayList<>(depth);
        final byte[] quals = new byte[readLength];
        Arrays.fill(quals, baseq);

        // add alt reads
        for (int i = 0; i < numAltReads; i++) {
            final GATKRead read = ArtificialReadUtils.createArtificialRead(samHeader, "Read" + i, chromosomeIndex, alignmentStart,
                    "CCCCACCCC".getBytes(), quals, "9M");
            read.setMappingQuality(mapq);
            if (i < numAltF1R2Reads){
                read.setIsFirstOfPair(); // set F1R2
            } else {
                read.setIsSecondOfPair(); // set F2R1
            }
            reads.add(read);
        }

        // REF reads have good read orientation balance
        for (int i = numAltReads; i < depth; i++) {
            final GATKRead read = ArtificialReadUtils.createArtificialRead(samHeader, "Read" + i, chromosomeIndex, alignmentStart,
                    "CCCCCCCCC".getBytes(), quals, "9M");
            read.setMappingQuality(mapq);

            if (i < numAltReads + numRefF1R2Reads){
                read.setIsFirstOfPair(); // set F1R2
            } else {
                read.setIsSecondOfPair(); // set F2R1
            }
            reads.add(read);
        }

        readMap.put(sampleName, reads);

        final ReadLikelihoods<Allele> readLikelihoods = new ReadLikelihoods<>(sampleList, alleleList, readMap);

        final LikelihoodMatrix<Allele> matrix = readLikelihoods.sampleMatrix(0);

        final double logLikelihoodOfBestAllele = 10.0;
        final int refAlleleIndex = 0;
        final int altAlleleIndex = 1;

        // Log likelihoods are initialized to 0 for all alleles.
        // For each read_i, set its likelihood for the alt allele to +10. Although this is not a valid likelihood, it's ok
        // as long as we get the relative values correct - the sample matrix will shift/normalize the values to valid
        // log likelihoods (e.g. < 0)
        for (int i = 0; i < numAltReads; i++) {
            matrix.set(altAlleleIndex, i, logLikelihoodOfBestAllele);
        }

        // Analogously, make read_j a ref read
        for (int j = numAltReads; j < depth; j++) {
            matrix.set(refAlleleIndex, j, logLikelihoodOfBestAllele);
        }

        genotypeBuilder.AD(new int[]{depth - numAltReads, numAltReads});

        return readLikelihoods;
    }
}