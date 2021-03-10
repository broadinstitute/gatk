package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.*;
import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.mutect.SomaticLikelihoodsEngine;
import org.broadinstitute.hellbender.tools.walkers.mutect.SubsettedLikelihoodMatrix;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import picard.util.MathUtil;

import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

@DocumentedFeature(groupName= HelpConstants.DOC_CAT_ANNOTATORS, groupSummary=HelpConstants.DOC_CAT_ANNOTATORS_SUMMARY, summary="Total depth of coverage per sample and over all samples (DP)")
public final class DirichletAlleleDepthAndFraction extends GenotypeAnnotation {

    private static DecimalFormat DEPTH_FORMAT = new DecimalFormat("#.##");
    private static DecimalFormat FRACTION_FORMAT = new DecimalFormat("#.####");

    public static final String DEPTH_KEY = "DD";
    public static final String FRACTION_KEY = "DF";

    public static final VCFFormatHeaderLine DEPTH_HEADER_LINE =
            new VCFFormatHeaderLine(DEPTH_KEY, VCFHeaderLineCount.R, VCFHeaderLineType.Float, "Allele depth based on Dirchlet posterior pseudo-counts");
    public static final VCFFormatHeaderLine FRACTION_HEADER_LINE =
            new VCFFormatHeaderLine(FRACTION_KEY, VCFHeaderLineCount.R, VCFHeaderLineType.Float, "Allele Fraction based on Dirchilet postrior pseudo-counts");

    @Argument(fullName = "dirichlet-prior-pseudo-count",
              doc = "Pseudo-count use as prior for all alleles. The default is 1.0 resulting in a flat prior")
    public double priorPseudoCount = 0.01;

    public boolean weightPseudoCounts = true;

    private static final List<String> KEYS = Arrays.asList(DEPTH_KEY, FRACTION_KEY);
    private static final List<VCFFormatHeaderLine> HEADER_LINES = Arrays.asList(DEPTH_HEADER_LINE, FRACTION_HEADER_LINE);

    private double[][] priorPseudoCounts = new double[10][]; // possibly won't ever go beyond 10.

    @Override
    public List<VCFFormatHeaderLine> getDescriptions() {
        return HEADER_LINES;
    }

    @Override
    public List<String> getKeyNames() {
        return KEYS;
    }

    @Override
    public void annotate(final ReferenceContext ref, final VariantContext vc,
                         final Genotype g, final GenotypeBuilder gb,
                         final AlleleLikelihoods<GATKRead, Allele> likelihoods) {
        if (likelihoods == null)
            return;
        final List<Allele> allelesToEmit = vc.getAlleles();
        final double[] prior = composePriorPseudoCounts(allelesToEmit.size());
        final LikelihoodMatrix<GATKRead, Allele> sampleMatrix = likelihoods.sampleMatrix(likelihoods.indexOfSample(g.getSampleName()));
        final LikelihoodMatrix<GATKRead, Allele> sampleMatrixForAlleles =
                allelesToEmit.size() == likelihoods.numberOfAlleles()
                        ? sampleMatrix : new SubsettedLikelihoodMatrix<>(sampleMatrix, allelesToEmit);
        final RealMatrix lkMatrixLog10 = sampleMatrixForAlleles.asRealMatrix();
        final RealMatrix lkMatrix = lkMatrixLog10.copy();
        for (int i = 0; i < lkMatrix.getRowDimension(); i++) {
            for (int j = 0; j < lkMatrix.getColumnDimension(); j++) {
                lkMatrix.setEntry(i, j, lkMatrix.getEntry(i, j) * Math.log(10));
            }
        }

        final double[] weights = calculateWeights(lkMatrix);
        final double[] posteriors =  sampleMatrix.evidenceCount() == 0 ? prior :
               SomaticLikelihoodsEngine.alleleFractionsPosterior(lkMatrix, prior, weights);
        final double[] frequencies = MathUtils.normalizeSumToOne(posteriors);
        gb.attribute(DEPTH_KEY, Arrays.stream(posteriors)
                .mapToObj(DEPTH_FORMAT::format)
                .collect(Collectors.joining(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR)));
        gb.attribute(FRACTION_KEY,Arrays.stream(frequencies)
                .mapToObj(FRACTION_FORMAT::format)
                .collect(Collectors.joining(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR)));
    }

    private double[] calculateWeights(RealMatrix lkMatrix) {
        if (!weightPseudoCounts) {
            return null;
        } else {
            final double[] result = new double[lkMatrix.getColumnDimension()];
            for (int i = 0; i < result.length; i++) {
                double best = lkMatrix.getEntry(0, i);
                double secondBest = Double.NEGATIVE_INFINITY;
                for (int j = 1; j < lkMatrix.getRowDimension(); j++) {
                    final double val = lkMatrix.getEntry(j, i);
                    if (val > best) {
                        secondBest = best;
                        best = val;
                    } else if (val > secondBest) {
                        secondBest = val;
                    }
                }
                final double secondAdjusted = secondBest - best;
                result[i] = 1 - Math.pow(10, secondAdjusted);
            }
            return result;
        }
    }

    /**
     * Composes prior pseudo-count arrays on demand.
     *
     * @param numberOfAlleles
     * @return never {@code null.}
     */
    private double[] composePriorPseudoCounts(final int numberOfAlleles) {
            if (priorPseudoCounts.length < numberOfAlleles) {
                priorPseudoCounts = Arrays.copyOf(priorPseudoCounts, numberOfAlleles << 1);
            }
            double[] result = priorPseudoCounts[numberOfAlleles];
            if (result == null) {
                Arrays.fill(result = priorPseudoCounts[numberOfAlleles] = new double[numberOfAlleles], priorPseudoCount);
            }
            return result;
    }
}
