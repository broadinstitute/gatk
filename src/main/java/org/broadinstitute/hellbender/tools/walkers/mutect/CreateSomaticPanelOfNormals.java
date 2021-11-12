package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.tools.walkers.readorientation.BetaDistributionShape;
import org.broadinstitute.hellbender.tools.walkers.validation.basicshortmutpileup.BetaBinomialDistribution;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.OptimizationUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import picard.cmdline.programgroups.VariantFilteringProgramGroup;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.function.DoubleUnaryOperator;
import java.util.stream.Collectors;

/**
 * Create a panel of normals (PoN) containing germline and artifactual sites for use with Mutect2.
 *
 * <p>
 *     The tool takes multiple normal sample callsets produced by {@link Mutect2}'s tumor-only mode and collates sites present in multiple samples
 *     (two by default, set by the --min-sample-count argument) into a sites-only VCF. The PoN captures common artifacts.  Mutect2 then
 *     uses the PoN to filter variants at the site-level.
 *
 *     The --max-germline-probability argument sets the threshold for possible germline variants to be included in the PoN.  By default this
 *     is set to 0.5, so that likely germline events are excluded.  This is usually the correct behavior as germline variants are best handled
 *     by probabilistic modeling via Mutect2's --germline-resource argument.  A germline resource, such as gnomAD in the case of humans, is a much
 *     more refined tool for germline filtering than any PoN could be.
 * </p>
 * <p>
 *     This tool is featured in the Somatic Short Mutation calling Best Practice Workflow.
 *     See <a href="https://software.broadinstitute.org/gatk/documentation/article?id=11136">Tutorial#11136</a> for a
 *     step-by-step description of the workflow and <a href="https://software.broadinstitute.org/gatk/documentation/article?id=11127">Article#11127</a>
 *     for an overview of what traditional somatic calling entails. For the latest pipeline scripts, see the
 *     <a href="https://github.com/broadinstitute/gatk/tree/master/scripts/mutect2_wdl">Mutect2 WDL scripts directory</a>.
 * </p>
 * <h3>Example workflow</h3>
 *
 * <h4>Step 1. Run Mutect2 in tumor-only mode for each normal sample.</h4>
 * <p>
 *   Note that as of May, 2019 -max-mnp-distance must be set to zero to avoid a bug in GenomicsDBImport.
 * </p>
 * <pre>
 * gatk Mutect2 -R reference.fasta -I normal1.bam -max-mnp-distance 0 -O normal1.vcf.gz
 * </pre>
 *
 * <h4>Step 2. Create a GenomicsDB from the normal Mutect2 calls.</h4>
 *
 *  <pre>
 *    gatk GenomicsDBImport -R reference.fasta -L intervals.interval_list \
 *       --genomicsdb-workspace-path pon_db \
 *       -V normal1.vcf.gz \
 *       -V normal2.vcf.gz \
 *       -V normal3.vcf.gz
 *  </pre>
 *
 * <h4>Step 3. Combine the normal calls using CreateSomaticPanelOfNormals.</h4>
 *
 * <pre>
 * gatk CreateSomaticPanelOfNormals -R reference.fasta -V gendb://pon_db -O pon.vcf.gz
 * </pre>
 */
@CommandLineProgramProperties(
        summary = "Make a panel of normals (PoN) for use with Mutect2",
        oneLineSummary = "Make a panel of normals for use with Mutect2",
        programGroup = VariantFilteringProgramGroup.class
)
@DocumentedFeature
@BetaFeature
public class CreateSomaticPanelOfNormals extends VariantWalker {

    public static final String MIN_SAMPLE_COUNT_LONG_NAME = "min-sample-count";
    public static final int DEFAULT_MIN_SAMPLE_COUNT = 2;

    public static final String MAX_GERMLINE_PROBABILITY_LONG_NAME = "max-germline-probability";
    public static final double DEFAULT_MAX_GERMLINE_PROBABILITY = 0.5;

    public static final String FRACTION_INFO_FIELD = "FRACTION";
    public static final String BETA_SHAPE_INFO_FIELD = "BETA";

    // the prior probability that any given site has an artifact
    private static final double ARTIFACT_PRIOR = 0.001;

    // beta distribution for artifact allele fractions
    private static final double ARTIFACT_ALPHA = 1;
    private static final double ARTIFACT_BETA = 7;

    private static final double NEGLIGIBLE_ALLELE_FREQUENCY = 1.0e-8;

    @Argument(fullName = MIN_SAMPLE_COUNT_LONG_NAME, doc="Number of samples containing a variant site required to include it in the panel of normals.", optional = true)
    private int minSampleCount = DEFAULT_MIN_SAMPLE_COUNT;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc="Output vcf")
    private String outputVcf;

    /**
     * A resource, such as gnomAD, containing population allele frequencies of common and rare variants.  We use this to remove germline variants from the panel
     * of normals, keeping only technical artifacts
     */
    @Argument(fullName= M2ArgumentCollection.GERMLINE_RESOURCE_LONG_NAME, doc="Population vcf of germline sequencing containing allele fractions.", optional = true)
    public FeatureInput<VariantContext> germlineResource;

    @Argument(fullName= MAX_GERMLINE_PROBABILITY_LONG_NAME, doc="Skip genotypes with germline probability greater than this value", optional = true)
    public double maxGermlineProbability = DEFAULT_MAX_GERMLINE_PROBABILITY;

    private VariantContextWriter vcfWriter;

    private int numSamples;

    @Override
    public void onTraversalStart() {
        final Set<VCFHeaderLine> headerInfo = new HashSet<>(getDefaultToolVCFHeaderLines());
        headerInfo.add(VCFHeader.makeHeaderVersionLine(VCFHeader.DEFAULT_VCF_VERSION));
        headerInfo.add(new VCFInfoHeaderLine(FRACTION_INFO_FIELD, 1, VCFHeaderLineType.Float, "Fraction of samples exhibiting artifact"));
        headerInfo.add(new VCFInfoHeaderLine(BETA_SHAPE_INFO_FIELD, 2, VCFHeaderLineType.Float, "Beta distribution parameters to fit artifact allele fractions"));

        getHeaderForVariants().getGenotypeSamples()
                .forEach(sample -> headerInfo.add(new VCFHeaderLine(Mutect2Engine.NORMAL_SAMPLE_KEY_IN_VCF_HEADER, sample)));

        vcfWriter = createVCFWriter(IOUtils.getPath(outputVcf));
        final VCFHeader outputHeader = new VCFHeader(headerInfo);
        outputHeader.setSequenceDictionary(getHeaderForVariants().getSequenceDictionary());
        vcfWriter.writeHeader(outputHeader);

        numSamples = getHeaderForVariants().getNGenotypeSamples();
    }

    @Override
    public void apply(final VariantContext vc, final ReadsContext rc, final ReferenceContext ref, final FeatureContext fc) {

        // if there are no non-span-del alt alleles, there's no artifact
        if (vc.getNAlleles() == 1) {
            return;
        } else if (vc.getNAlleles() == 2 && vc.getAlternateAllele(0).basesMatch(Allele.SPAN_DEL)) {
            return;
        }

        final List<VariantContext> germline = fc.getValues(germlineResource, new SimpleInterval(vc));
        final double germlineAF = germline.isEmpty() ? 0 : Mutect2Engine.getAttributeAsDoubleList(germline.get(0), VCFConstants.ALLELE_FREQUENCY_KEY, 0.0)
                .stream().mapToDouble(af -> af).sum();

        // note: if at this site some input vcfs had a variant and some had only a spanning deletion from an upstream event,
        // GenomicsDBImport removes the spanning deletion from the ADs and so the altCount logic works and counts only
        // real variants, not spanning deletions.

        // TODO: replace this logic where multiallelic sites are excluded from germline analysis with something principled
        // TODO that calculates germline probabilities by allele.  See GATK issue #6744 https://github.com/broadinstitute/gatk/issues/6744
        // TODO for some analysis as to why the AD field may not be available barring changes to GenomicsDB
        final List<Genotype> variantGenotypes = vc.getNAlleles() > 2 ? vc.getGenotypes() : vc.getGenotypes().stream()
                .filter(g -> hasArtifact(g, germlineAF)).collect(Collectors.toList());

        if (variantGenotypes.size() < minSampleCount) {
            return;
        }

        final double fraction = (double) variantGenotypes.size() / numSamples;

        final List<int[]> altAndRefCounts = variantGenotypes.stream()
                .filter(Genotype::hasAD)
                .map(g -> new int[] {altCount(g), g.getAD()[0]})
                .collect(Collectors.toList());

        final BetaDistributionShape betaDistributionShape = fitBeta(altAndRefCounts);

        final VariantContext outputVc = new VariantContextBuilder(vc.getSource(), vc.getContig(), vc.getStart(), vc.getEnd(), vc.getAlleles())
                .attribute(FRACTION_INFO_FIELD, fraction)
                .attribute(BETA_SHAPE_INFO_FIELD, new double[] { betaDistributionShape.getAlpha(), betaDistributionShape.getBeta()})
                .make();

        vcfWriter.add(outputVc);
    }

    private final boolean hasArtifact(final Genotype g, final double populationAlleleFrequency) {
        final int altCount = altCount(g);
        if (altCount == 0) {
            return false;
        }
        final int totalCount = (int) MathUtils.sum(g.getAD());

        return germlineProbability(populationAlleleFrequency, altCount, totalCount) < maxGermlineProbability;
    }

    private static final int altCount(final Genotype g) {
        return g.hasAD() ? (int) MathUtils.sum(g.getAD()) - g.getAD()[0] : 0;
    }

    private static final double germlineProbability(final double alleleFrequency, final int altCount, final int totalCount) {
        if (alleleFrequency < NEGLIGIBLE_ALLELE_FREQUENCY || alleleFrequency > 1) {
            return 0;
        }

        final double hetPrior = alleleFrequency * (1 - alleleFrequency) * 2;
        final double homPrior = MathUtils.square(alleleFrequency);
        final double hetLikelihood = MathUtils.binomialProbability(totalCount, altCount, 0.5);
        final double homLikelihood = MathUtils.binomialProbability(totalCount, altCount, 0.98);

        final double artifactLikelihood = new BetaBinomialDistribution(null, ARTIFACT_ALPHA, ARTIFACT_BETA, totalCount).probability(altCount);

        // the indices of this array are 0 -- germline het, 1 -- germline hom alt, 2 -- not germline
        final double[] relativeProbsOfHetHomArtifact = {hetPrior * hetLikelihood + homPrior * homLikelihood, ARTIFACT_PRIOR * artifactLikelihood};

        // check for invalid probabilities just in case of finite precision error
        return MathUtils.sum(relativeProbsOfHetHomArtifact) < 0 ? 0 : MathUtils.normalizeSumToOne(relativeProbsOfHetHomArtifact)[0];
    }

    private BetaDistributionShape fitBeta(final List<int[]> altAndRefCounts) {
        final int totalAltCount = altAndRefCounts.stream().mapToInt(pair -> pair[0]).sum();
        final int totalRefCount = altAndRefCounts.stream().mapToInt(pair -> pair[1]).sum();
        final int min = Math.min(totalAltCount, totalRefCount);

        // keeping the ratio of alpha and beta equal to the ratio of baseAlpha and baseBeta gives the empirical mean
        final double baseAlpha = (totalAltCount + 1.0) / (min + 1);
        final double baseBeta = (totalRefCount + 1.0) / (min + 1);

        final DoubleUnaryOperator logLikelihood = s -> {
            final double alpha = baseAlpha * s;
            final double beta = baseBeta * s;

            return altAndRefCounts.stream().mapToDouble(pair -> {
                final int n = pair[0] + pair[1];
                final int k = pair[0];
                return new BetaBinomialDistribution(null, alpha, beta, n).logProbability(k);
            }).sum();
        };

        final double scale = OptimizationUtils.max(logLikelihood, 0.01, 100, 1, 0.01, 0.1, 100).getPoint();

        return new BetaDistributionShape(baseAlpha * scale, baseBeta * scale);
    }

    @Override
    public Object onTraversalSuccess() {
        return "SUCCESS";
    }

    @Override
    public void closeTool() {
        if ( vcfWriter != null ) {
            vcfWriter.close();
        }
    }
}
