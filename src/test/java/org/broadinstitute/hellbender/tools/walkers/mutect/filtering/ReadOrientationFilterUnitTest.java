package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceFileSource;
import org.broadinstitute.hellbender.tools.walkers.readorientation.*;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Optional;

public class ReadOrientationFilterUnitTest extends GATKBaseTest {
    final int chromosomeIndex = 1;
    final int variantPosition = 100_000;
    final String sampleName = "someTumor";

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
                {100, 0.01, 0.8, 1.0, 1e-1, Optional.of(ReadOrientation.F1R2)}, // medium depth hom ref site with some error
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

        final int altCount = (int) (depth * alleleFraction);
        final int refCount = depth - altCount;
        final int altF1R2 = (int) (altCount * altF1R2Fraction);
        final int[] f1r2 = new int[] { refCount/2, altF1R2};
        final int[] f2r1 = new int[] { refCount - refCount/2, altCount - altF1R2};

        // Arbitrarily choose the number of ref and alt examples observed to compute the priors. They do not affect inference.
        int numRefExamples = 1_000_000;
        int numAltExamples = 1000;


        final Nucleotide[] altBases = new Nucleotide[] {Nucleotide.A, Nucleotide.C, Nucleotide.A};
        final StringBuilder altBasesStringBuilder = new StringBuilder();
        for (final Nucleotide nuc : altBases) {
            altBasesStringBuilder.append((char) nuc.encodeAsByte());
        }
        final String altBasesString = altBasesStringBuilder.toString();

        final ArtifactPriorCollection priors = new ArtifactPriorCollection(sampleName);

        for (int n = 0; n < altBases.length; n++) {
            // most of the time a site is hom ref, and the remaining probability mass is distributed equally among the remaining states
            double[] pi = createReasonablePrior(altBases[n]);
            priors.set(new ArtifactPrior(reference.getKmerAround(variantPosition + n, F1R2FilterConstants.REF_CONTEXT_PADDING), pi, numRefExamples, numAltExamples));
        }
        final File table = IOUtils.createTempFile("prior", "table");
        priors.writeArtifactPriors(table);

        final ReadOrientationFilter filter = new ReadOrientationFilter(Collections.singletonList(table));

        for (int mnpLength = 1; mnpLength <= altBases.length; mnpLength++) {
            final Allele refAllele = Allele.create(reference.getBases(new SimpleInterval("1", variantPosition, variantPosition + mnpLength - 1)), true);
            final Allele altAllele = Allele.create(altBasesString.substring(0, mnpLength), false); // C -> A transition
            final List<Allele> alleles = Arrays.asList(refAllele, altAllele);
            final VariantContext vc = new VariantContextBuilder("source", Integer.toString(chromosomeIndex), variantPosition, variantPosition + mnpLength - 1, alleles)
                    .attribute(GATKVCFConstants.TUMOR_LOG_10_ODDS_KEY, new double[]{7.0})
                    .make();

            final GenotypeBuilder genotypeBuilder = new GenotypeBuilder(sampleName)
                    .attribute(GATKVCFConstants.F1R2_KEY, f1r2)
                    .attribute(GATKVCFConstants.F2R1_KEY, f2r1);
            genotypeBuilder.alleles(alleles);
            final double posterior = filter.artifactProbability(reference, vc, genotypeBuilder.make());

            if (expectedArtifactType.isPresent()) {
                Assert.assertEquals(posterior, expectedArtifactProb, epsilon);
            } else {
                Assert.assertTrue(posterior < 0.01);
            }
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

        return MathUtils.normalizeSumToOne(pi);
    }

}