package org.broadinstitute.hellbender.tools.walkers.annotator;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.help.HelpConstants;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.IntStream;


/**
 * Count of all reads with MAPQ = 0 across all samples
 *
 * <p>This anotation gives you the count of all reads that have MAPQ = 0 across all samples. The count of reads with MAPQ0 can be used for quality control; high counts typically indicate regions where it is difficult to make confident calls.</p>
 *
 * <h3>Caveat</h3>
 * <p>It is not useful to apply this annotation with HaplotypeCaller because HC filters out all reads with MQ0 upfront, so the annotation will always return a value of 0.</p>
 *
 * <h3>Related annotations</h3>
 * <ul>
 *     <li><b><a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_annotator_MappingQualityZeroBySample.php">MappingQualityZeroBySample</a></b> gives the count of reads with MAPQ=0 for each individual sample.</li>
 *     <li><b><a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_annotator_LowMQ.php">LowMQ</a></b> gives the proportion of reads with low mapping quality (MAPQ below 10, including 0).</li>
 * </ul>
 *
 */
@DocumentedFeature(groupName=HelpConstants.DOC_CAT_ANNOTATORS, groupSummary=HelpConstants.DOC_CAT_ANNOTATORS_SUMMARY, summary="Count of all reads with MAPQ = 0 across all samples (MQ0)")
public final class MappingQualityZero extends InfoFieldAnnotation {

    @Override
    public Map<String, Object> annotate(final ReferenceContext ref,
                                        final VariantContext vc,
                                        final ReadLikelihoods<Allele> likelihoods) {
        Utils.nonNull(vc);
        if (!vc.isVariant() || likelihoods == null){
            return Collections.emptyMap();
        }
        //NOTE: unlike other annotations, this one returns 0 if likelihoods are empty
        final long mq0 = IntStream.range(0, likelihoods.numberOfSamples()).boxed()
                .flatMap(s -> likelihoods.sampleReads(s).stream())
                .filter(r -> r.getMappingQuality() == 0)
                .count();

        return Collections.singletonMap(getKeyNames().get(0), formattedValue(mq0));
    }

    @VisibleForTesting
    static String formattedValue(long mq0) {
        return String.format("%d", mq0);
    }

    @Override
    public List<String> getKeyNames() { return Collections.singletonList(VCFConstants.MAPPING_QUALITY_ZERO_KEY); }

    @Override
    public List<VCFInfoHeaderLine> getDescriptions() {
        return Collections.singletonList(VCFStandardHeaderLines.getInfoLine(getKeyNames().get(0)));
    }
}
