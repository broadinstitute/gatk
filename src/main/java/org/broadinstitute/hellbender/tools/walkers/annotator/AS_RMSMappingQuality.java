package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.apache.log4j.Logger;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;


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
                                        final Map<String, PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap) {
        return annotateRawData(ref, vc, perReadAlleleLikelihoodMap);
    }

    @Override
    public Map<String, Object> annotateRawData(final ReferenceContext ref,
                                               final VariantContext vc,
                                               final Map<String, PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap ) {
        Utils.nonNull(vc);
        if ( perReadAlleleLikelihoodMap == null) {
            return null;
        }

        final Map<String, Object> annotations = new HashMap<>();
        final ReducibleAnnotationData<Number> myData = new ReducibleAnnotationData<>(null);
        calculateRawData(vc, perReadAlleleLikelihoodMap, myData);
        final String annotationString = makeRawAnnotationString(vc.getAlleles(), myData.getAttributeMap());
        annotations.put(getRawKeyName(), annotationString);
        return annotations;
    }

    @SuppressWarnings({"unchecked", "rawtypes"})//FIXME
    @Override
    public void calculateRawData(final VariantContext vc,
                                 final Map<String, PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap,
                                 final ReducibleAnnotationData myData){
        //For the raw data here, we're only keeping track of the sum of the squares of our values
        //When we go to reduce, we'll use the AD info to get the number of reads

        //must use perReadAlleleLikelihoodMap for allele-specific annotations
        if (perReadAlleleLikelihoodMap == null || perReadAlleleLikelihoodMap.isEmpty()) {
            return;
        }
        getRMSDataFromPRALM(perReadAlleleLikelihoodMap, myData);
    }

    @Override
    public List<String> getKeyNames() { return Arrays.asList(GATKVCFConstants.AS_RMS_MAPPING_QUALITY_KEY); }

    @Override
    public String getRawKeyName() { return GATKVCFConstants.AS_RAW_RMS_MAPPING_QUALITY_KEY; }

    public void getRMSDataFromPRALM(Map<String, PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap, ReducibleAnnotationData<Number> myData) {
       //over all the samples in the Map...
        for ( final PerReadAlleleLikelihoodMap perReadLikelihoods : perReadAlleleLikelihoodMap.values() ) {
            //for each read...
            for ( final Map.Entry<GATKRead, Map<Allele,Double>> readLikelihoods : perReadLikelihoods.getLikelihoodReadMap().entrySet() ) {
                final int mq = readLikelihoods.getKey().getMappingQuality();
                if ( mq != QualityUtils.MAPPING_QUALITY_UNAVAILABLE ) {
                    if (!PerReadAlleleLikelihoodMap.getMostLikelyAllele(readLikelihoods.getValue()).isInformative()) {
                        continue;
                    }
                    final Allele bestAllele = PerReadAlleleLikelihoodMap.getMostLikelyAllele(readLikelihoods.getValue()).getMostLikelyAllele();
                    double currSquareSum = 0;
                    if (myData.hasAttribute(bestAllele)) {
                        currSquareSum += (double) myData.getAttribute(bestAllele);
                    }
                    myData.putAttribute(bestAllele, currSquareSum + mq * mq);
                }
            }
        }
    }

    public String makeRawAnnotationString(final List<Allele> vcAlleles, final Map<Allele, Number> perAlleleValues) {
        String annotationString = "";
        for (final Allele current : vcAlleles) {
            if (!annotationString.isEmpty())
                annotationString += PRINT_DELIM;
            if(perAlleleValues.get(current) != null)
                annotationString += String.format(printFormat,perAlleleValues.get(current));
            else
                annotationString += String.format(printFormat, 0.0);
        }
        return annotationString;
    }
}