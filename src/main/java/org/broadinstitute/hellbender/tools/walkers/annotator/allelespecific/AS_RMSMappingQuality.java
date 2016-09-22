package org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.apache.log4j.Logger;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.annotator.InfoFieldAnnotation;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;

import java.util.*;


/**
 * Allele-specific Root Mean Square of the mapping quality of reads across all samples.
 *
 * <p>This annotation provides an estimation of the mapping quality of reads supporting each alternate allele in a variant call. Depending on the tool it is called from, it produces either raw data (sum of squared MQs) or the calculated root mean square.</p>
 *
 * The raw data is used to accurately calculate the root mean square when combining more than one sample.
 *
 * <h3>Statistical notes</h3>
 * <p>The root mean square is equivalent to the mean of the mapping qualities plus the standard deviation of the mapping qualities.</p>
 *
 * <h3>Related annotations</h3>
 * <ul>
 *     <li><b><a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_annotator_MappingQualityRankSumTest.php">MappingQualityRankSumTest</a></b> compares the mapping quality of reads supporting the REF and ALT alleles.</li>
 * </ul>
 *
 * <h3>Caveat</h3>
 * <p>Uninformative reads are not used in this annotation.</p>
 *
 * <h3>Related annotations</h3>
 * <ul>
 *     <li><b><a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_annotator_RMSMappingQuality.php">RMSMappingQuality</a></b> outputs a version of this annotation that includes all alternate alleles in a single calculation.</li>
 *     <li><b><a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_annotator_MappingQualityRankSumTest.php">MappingQualityRankSumTest</a></b> compares the mapping quality of reads supporting the REF and ALT alleles.</li>
 * </ul>
 */
public final class AS_RMSMappingQuality extends InfoFieldAnnotation implements AS_StandardAnnotation, ReducibleAnnotation {

    private final String printFormat = "%.2f";

    private static final Logger logger = Logger.getLogger(AS_RMSMappingQuality.class);
    public static final String SPLIT_DELIM = "\\|"; //String.split takes a regex, so we need to escape the pipe
    public static final String PRINT_DELIM = "|";

    @Override
    public List<VCFInfoHeaderLine> getDescriptions() {
        //TODO I'm just going to ignore this for now
//        if (AnnotationUtils.walkerRequiresRawData(callingWalker))
            return Arrays.asList(GATKVCFHeaderLines.getInfoLine(getRawKeyName()));
//        else
//            return Arrays.asList(GATKVCFHeaderLines.getInfoLine(getKeyNames().get(0)));
    }

    @Override
    public Map<String, Object> annotate(final ReferenceContext ref,
                                        final VariantContext vc,
                                        final ReadLikelihoods<Allele> likelihoods) {
        return annotateRawData(ref, vc, likelihoods);
    }

    @Override
    public Map<String, Object> annotateRawData(final ReferenceContext ref,
                                               final VariantContext vc,
                                               final ReadLikelihoods<Allele> likelihoods ) {
        Utils.nonNull(vc);
        if ( likelihoods == null) {
            return Collections.emptyMap();
        }

        final Map<String, Object> annotations = new LinkedHashMap<>();
        final ReducibleAnnotationData<Number> myData = new ReducibleAnnotationData<>(null);
        calculateRawData(vc, likelihoods, myData);
        final String annotationString = makeRawAnnotationString(vc.getAlleles(), myData.getAttributeMap());
        annotations.put(getRawKeyName(), annotationString);
        return annotations;
    }

    @SuppressWarnings({"unchecked", "rawtypes"})//FIXME
    @Override
    public void calculateRawData(final VariantContext vc,
                                 final ReadLikelihoods<Allele> likelihoods,
                                 final ReducibleAnnotationData myData){
        //For the raw data here, we're only keeping track of the sum of the squares of our values
        //When we go to reduce, we'll use the AD info to get the number of reads

        //must use likelihoods for allele-specific annotations
        if (likelihoods == null) {
            return;
        }
        getRMSDataFromLikelihoods(likelihoods, myData);
    }

    @Override
    public List<String> getKeyNames() { return Arrays.asList(GATKVCFConstants.AS_RMS_MAPPING_QUALITY_KEY); }

    @Override
    public String getRawKeyName() { return GATKVCFConstants.AS_RAW_RMS_MAPPING_QUALITY_KEY; }

    private void getRMSDataFromLikelihoods(final ReadLikelihoods<Allele> likelihoods, ReducibleAnnotationData<Number> myData) {
        for ( final ReadLikelihoods<Allele>.BestAllele bestAllele : likelihoods.bestAlleles() ) {
            if (bestAllele.isInformative()) {
                final int mq = bestAllele.read.getMappingQuality();
                if ( mq != QualityUtils.MAPPING_QUALITY_UNAVAILABLE ) {
                    final double currSquareSum = myData.hasAttribute(bestAllele.allele) ? (double) myData.getAttribute(bestAllele.allele) : 0;
                    myData.putAttribute(bestAllele.allele, currSquareSum + mq * mq);
                }
            }
        }
    }

    private String makeRawAnnotationString(final List<Allele> vcAlleles, final Map<Allele, Number> perAlleleValues) {
        String annotationString = "";
        for (final Allele current : vcAlleles) {
            if (!annotationString.isEmpty()) {
                annotationString += PRINT_DELIM;
            }
            if(perAlleleValues.get(current) != null) {
                annotationString += String.format(printFormat, perAlleleValues.get(current));
            } else {
                annotationString += String.format(printFormat, 0.0);
            }
        }
        return annotationString;
    }
}
