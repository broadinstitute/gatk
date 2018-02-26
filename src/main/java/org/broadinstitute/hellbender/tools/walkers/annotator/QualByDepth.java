package org.broadinstitute.hellbender.tools.walkers.annotator;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;

import java.util.Collections;
import java.util.List;
import java.util.Map;

/**
 * Variant confidence normalized by unfiltered depth of variant samples
 *
 * <p>This annotation puts the variant confidence QUAL score into perspective by normalizing for the amount of coverage available. Because each read contributes a little to the QUAL score, variants in regions with deep coverage can have artificially inflated QUAL scores, giving the impression that the call is supported by more evidence than it really is. To compensate for this, we normalize the variant confidence by depth, which gives us a more objective picture of how well supported the call is.</p>
 *
 * <h3>Statistical notes</h3>
 * <p>The QD is the QUAL score normalized by allele depth (AD) for a variant. For a single sample, the HaplotypeCaller calculates the QD by taking QUAL/AD. For multiple samples, HaplotypeCaller and GenotypeGVCFs calculate the QD by taking QUAL/AD of samples with a non hom-ref genotype call. The reason we leave out the samples with a hom-ref call is to not penalize the QUAL for the other samples with the variant call.</p>
 * <p>Here is a single sample example:</p>
 * <pre>2	37629	.	C	G	1063.77	.	AC=2;AF=1.00;AN=2;DP=31;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=58.50;QD=34.32;SOR=2.376	GT:AD:DP:GQ:PL:QSS	1/1:0,31:31:93:1092,93,0:0,960</pre>
   <p>QUAL/AD = 1063.77/31 = 34.32 = QD</p>
 * <p>Here is a multi-sample example:</p>
 * <pre>10	8046	.	C	T	4107.13	.	AC=1;AF=0.167;AN=6;BaseQRankSum=-3.717;DP=1063;FS=1.616;MLEAC=1;MLEAF=0.167;QD=11.54
   GT:AD:DP:GQ:PL:QSS	0/0:369,4:373:99:0,1007,12207:10548,98	    0/0:331,1:332:99:0,967,11125:9576,27	    0/1:192,164:356:99:4138,0,5291:5501,4505</pre>
 * <p>QUAL/AD = 4107.13/356 = 11.54 = QD</p>
 *
 * <h3>Caveat</h3>
 * <p>This annotation can only be calculated for sites for which at least one sample was genotyped as carrying a variant allele.</p>
 *
 * <h3>Related annotations</h3>
 * <ul>
 *     <li><b><a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_annotator_Coverage.php">Coverage</a></b> gives the filtered depth of coverage for each sample and the unfiltered depth across all samples.</li>
 *     <li><b><a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_annotator_DepthPerAlleleBySample.php">DepthPerAlleleBySample</a></b> calculates depth of coverage for each allele per sample (AD).</li>
 * </ul>
 */
@DocumentedFeature(groupName=HelpConstants.DOC_CAT_ANNOTATORS, groupSummary=HelpConstants.DOC_CAT_ANNOTATORS_SUMMARY, summary="Variant confidence normalized by unfiltered depth of variant samples (QD)")
public final class QualByDepth extends InfoFieldAnnotation implements StandardAnnotation {

    @VisibleForTesting
    static final double MAX_QD_BEFORE_FIXING = 35;

    @VisibleForTesting
    static final double IDEAL_HIGH_QD = 30;
    private static final double JITTER_SIGMA = 3;

    @Override
    public Map<String, Object> annotate(final ReferenceContext ref,
                                        final VariantContext vc,
                                        final ReadLikelihoods<Allele> likelihoods) {
        Utils.nonNull(vc);
        if ( !vc.hasLog10PError() ) {
            return Collections.emptyMap();
        }

        final GenotypesContext genotypes = vc.getGenotypes();
        if ( genotypes == null || genotypes.isEmpty() ) {
            return Collections.emptyMap();
        }

        int depth = 0;
        int ADrestrictedDepth = 0;

        for ( final Genotype genotype : genotypes ) {
            // we care only about variant calls with likelihoods
            if ( !genotype.isHet() && !genotype.isHomVar() ) {
                continue;
            }

            // if we have the AD values for this sample, let's make sure that the variant depth is greater than 1!
            if ( genotype.hasAD() ) {
                final int[] AD = genotype.getAD();
                final int totalADdepth = (int) MathUtils.sum(AD);
                if ( totalADdepth != 0 ) {
                    if (totalADdepth - AD[0] > 1) {
                        ADrestrictedDepth += totalADdepth;
                    }
                    depth += totalADdepth;
                    continue;
                }
            }
            // if there is no AD value or it is a dummy value, we want to look to other means to get the depth
            if (likelihoods != null) {
                depth += likelihoods.sampleReadCount(likelihoods.indexOfSample(genotype.getSampleName()));
            } else if ( genotype.hasDP() ) {
                depth += genotype.getDP();
            }
        }

        // if the AD-restricted depth is a usable value (i.e. not zero), then we should use that one going forward
        if ( ADrestrictedDepth > 0 ) {
            depth = ADrestrictedDepth;
        }

        if ( depth == 0 ) {
            return Collections.emptyMap();
        }

        final double qual = -10.0 * vc.getLog10PError();
        double QD = qual / depth;

        // Hack: see note in the fixTooHighQD method below
        QD = fixTooHighQD(QD);

        return Collections.singletonMap(getKeyNames().get(0), String.format("%.2f", QD));
    }

    /**
     * The haplotype caller generates very high quality scores when multiple events are on the
     * same haplotype.  This causes some very good variants to have unusually high QD values,
     * and VQSR will filter these out.  This code looks at the QD value, and if it is above
     * threshold we map it down to the mean high QD value, with some jittering
     *
     * @param QD the raw QD score
     * @return a QD value
     */
    public static double fixTooHighQD(final double QD) {
        if ( QD < MAX_QD_BEFORE_FIXING ) {
            return QD;
        } else {
            return IDEAL_HIGH_QD + Utils.getRandomGenerator().nextGaussian() * JITTER_SIGMA;
        }
    }

    @Override
    public List<String> getKeyNames() { return Collections.singletonList(GATKVCFConstants.QUAL_BY_DEPTH_KEY); }
}
