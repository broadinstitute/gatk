package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.variant.variantcontext.*;
import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.RandomGeneratorFactory;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceFileSource;
import org.broadinstitute.hellbender.tools.walkers.readorientation.*;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.genotyper.*;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.*;

public class ReadOrientationArtifactUnitTest extends GATKBaseTest {
    final int chromosomeIndex = 1;
    final int variantPosition = 100_000;
    final String sampleName = "someTumor";
    final int readLength = 9;

    @DataProvider(name = "readOrientationArtifact")
    public Object[][] readOrientationArtifactData() {
        return new Object[][]{
                {150, 0.5, 0.5, 0.0, 1e-3, Optional.empty()}, // a likely germline het with no bias
                {150, 0.3, 0.5, 0.0, 1e-3, Optional.empty()}, // allele fraction 0.3 variant with no bias
                {150, 0.3, 0.99, 1.0, 1e-3, Optional.of(ReadOrientation.F1R2)}, // allele fraction 0.3 variant with a heavy F1R2 bias
                {150, 0.3, 0.06, 1.0, 1e-2, Optional.of(ReadOrientation.F2R1)}, // allele fraction 0.2 variant with a moderate F2R1 bias
                {150, 0.2, 0.01, 1.0, 1e-2, Optional.of(ReadOrientation.F2R1)}, // allele fraction 0.2 variant with a heavy F2R1 bias
                {150, 0.03, 0.0, 1.0, 1e-1, Optional.of(ReadOrientation.F2R1)}, // low allele fraction with a heavy F2R1 bias
                {300, 0.02, 1.0, 1.0, 1e-1, Optional.of(ReadOrientation.F1R2)}, // low allele fraction, high depth
                {200, 0.02, 1.0, 1.0, 1e-1, Optional.of(ReadOrientation.F1R2)}, // low allele fraction, high depth - NOTE: to get this case right, *maybe* we should learn f.
                {700, 0.025, 1.0, 1.0, 1e-1, Optional.of(ReadOrientation.F1R2)}, // high depth hom ref site with some error
                {50, 0.01, 0.8, 1.0, 1e-1, Optional.of(ReadOrientation.F1R2)}, // medium depth hom ref site with some error
                {160, 0.025, 1.0, 1.0, 1e-1, Optional.of(ReadOrientation.F1R2)},
                {150, 0.025, 1.0, 1.0, 1e-1, Optional.of(ReadOrientation.F1R2)},
                {700, 0.01, 1.0, 1.0, 1e-1, Optional.of(ReadOrientation.F1R2)}};
    }

    @DataProvider(name = "failingCases")
    public Object[][] failingCasesData() {
        return new Object[][]{
                {72, 0.041, 1.0, 1.0, 1e-1, Optional.of(ReadOrientation.F1R2)}};
    }

    @Test(dataProvider="readOrientationArtifact")
    public void test(final int depth, final double alleleFraction, final double altF1R2Fraction,
                     final double expectedArtifactProb, final double epsilon,
                     final Optional<ReadOrientation> expectedArtifactType) throws IOException {
        final ReferenceContext reference = new ReferenceContext(new ReferenceFileSource(IOUtils.getPath(hg19_chr1_1M_Reference)),
                new SimpleInterval("1", variantPosition -5, variantPosition +5));
        final String refContext = reference.getKmerAround(variantPosition, F1R2FilterConstants.REF_CONTEXT_PADDING);
        final Allele refAllele = Allele.create(F1R2FilterUtils.getMiddleBase(refContext).toString(), true);
        final Nucleotide altBase = Nucleotide.A;
        final Allele altAllele = Allele.create(altBase.toString(), false); // C -> A transition
        final List<Allele> alleles = Arrays.asList(refAllele, altAllele);
        final VariantContext vc  = new VariantContextBuilder("source", Integer.toString(chromosomeIndex), variantPosition, variantPosition, alleles)
                .attribute(GATKVCFConstants.TUMOR_LOD_KEY, new double[]{ 7.0 })
                .make();

        /**
         * The prior distribution of the states predicts that most of the time a site is hom ref, and the remaining
         * probability mass is distributed equally among the remaining states
         */
        double[] pi = createReasonablePrior(altBase);

        // Arbitrarily choose the number of ref and alt examples observed to compute the priors. They do not affect inference.
        int numRefExamples = 1_000_000;
        int numAltExamples = 1000;

        final ArtifactPriorCollection priors = new ArtifactPriorCollection();
        priors.set(new ArtifactPrior(refContext, pi, numRefExamples, numAltExamples));
        final File table = File.createTempFile("prior", "table");
        priors.writeArtifactPriors(table);

        final ReadOrientationArtifact annotation = new ReadOrientationArtifact(table);

        final GenotypeBuilder genotypeBuilder = new GenotypeBuilder(sampleName);
        final ReadLikelihoods<Allele> readLikelihoods = createReadLikelihoods(depth, alleleFraction, altF1R2Fraction,
                alleles, genotypeBuilder);
        genotypeBuilder.alleles(alleles);
        annotation.annotate(reference, vc, genotypeBuilder.make(), genotypeBuilder, readLikelihoods);

        final Genotype genotype = genotypeBuilder.make();
        final double posterior = GATKProtectedVariantContextUtils.getAttributeAsDouble(genotype, GATKVCFConstants.ROF_POSTERIOR_KEY, -1.0);

        if (expectedArtifactType.isPresent()){
            final ReadOrientation artifactType = ReadOrientation.valueOf(
                    GATKProtectedVariantContextUtils.getAttributeAsString(genotype, GATKVCFConstants.ROF_TYPE_KEY, null));

            Assert.assertEquals(posterior, expectedArtifactProb, epsilon);
            Assert.assertEquals(artifactType, expectedArtifactType.get());
        } else {
            Assert.assertTrue(posterior < 0.01);
        }
    }

    private double[] createReasonablePrior(final Nucleotide altAllele) {
        double[] pi = new double[ArtifactState.values().length];
        pi[ArtifactState.HOM_REF.ordinal()] = 1000;
        pi[ArtifactState.GERMLINE_HET.ordinal()] = 10;
        pi[ArtifactState.SOMATIC_HET.ordinal()] = 1;
        pi[ArtifactState.HOM_VAR.ordinal()] = 1;

        pi[ArtifactState.getF1R2StateForAlt(altAllele).ordinal()] = 1;
        pi[ArtifactState.getF2R1StateForAlt(altAllele).ordinal()] = 1;

        return MathUtils.normalizeFromRealSpace(pi);
    }

    /** Update the genotype builder by side effect, too **/
    private ReadLikelihoods<Allele> createReadLikelihoods(final int depth, final double alleleFraction, final double altF1R2Fraction,
                                                          final List<Allele> alleles, final GenotypeBuilder genotypeBuilder){
        final RandomGenerator randomGenerator = RandomGeneratorFactory.createRandomGenerator(new Random(111));

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
        final int alignmentStart = variantPosition - readLength/2;

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