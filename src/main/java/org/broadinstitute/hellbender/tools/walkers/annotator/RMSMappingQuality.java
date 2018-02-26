package org.broadinstitute.hellbender.tools.walkers.annotator;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.ReducibleAnnotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.ReducibleAnnotationData;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;

import java.util.*;
import java.util.stream.IntStream;


/**
 * Root Mean Square of the mapping quality of reads across all samples.
 *
 * <p>This annotation provides an estimation of the overall mapping quality of reads supporting a variant call, averaged over all samples in a cohort.</p>
 *
 * <h3>Statistical notes</h3>
 * <p>The root mean square is equivalent to the mean of the mapping qualities plus the standard deviation of the mapping qualities.</p>
 *
 * <h3>Related annotations</h3>
 * <ul>
 *     <li><b><a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_annotator_MappingQualityRankSumTest.php">MappingQualityRankSumTest</a></b> compares the mapping quality of reads supporting the REF and ALT alleles.</li>
 * </ul>
 *
 */
@DocumentedFeature(groupName=HelpConstants.DOC_CAT_ANNOTATORS, groupSummary=HelpConstants.DOC_CAT_ANNOTATORS_SUMMARY, summary="Root mean square of the mapping quality of reads across all samples (MQ)")
public final class RMSMappingQuality extends InfoFieldAnnotation implements StandardAnnotation, ReducibleAnnotation {
    private static final RMSMappingQuality instance = new RMSMappingQuality();

    @Override
    public String getRawKeyName() { return GATKVCFConstants.RAW_RMS_MAPPING_QUALITY_KEY;}

    @Override
    public List<String> getKeyNames() {
        return Arrays.asList(VCFConstants.RMS_MAPPING_QUALITY_KEY, getRawKeyName());
    }

    @Override
    public List<VCFInfoHeaderLine> getDescriptions() {
        return Arrays.asList(VCFStandardHeaderLines.getInfoLine(getKeyNames().get(0)), GATKVCFHeaderLines.getInfoLine(getRawKeyName()));
    }

    @Override
    public List<VCFInfoHeaderLine> getRawDescriptions() {
        return getDescriptions();
    }

    /**
     * Generate the raw data necessary to calculate the annotation. Raw data is the final endpoint for gVCFs.
     */
    @Override
    public Map<String, Object> annotateRawData(final ReferenceContext ref,
                                               final VariantContext vc,
                                               final ReadLikelihoods<Allele> likelihoods){
        Utils.nonNull(vc);
        if (likelihoods == null || likelihoods.readCount() == 0) {
            return Collections.emptyMap();
        }

        final Map<String, Object> annotations = new HashMap<>();
        final ReducibleAnnotationData<Number> myData = new ReducibleAnnotationData<>(null);
        calculateRawData(vc, likelihoods, myData);
        final String annotationString = formattedValue((double) myData.getAttributeMap().get(Allele.NO_CALL));
        annotations.put(getRawKeyName(), annotationString);
        return annotations;
    }

    @Override
    @SuppressWarnings({"unchecked", "rawtypes"})//FIXME generics here blow up
    public Map<String, Object> combineRawData(final List<Allele> vcAlleles, final List<ReducibleAnnotationData<?>>  annotationList) {
        //VC already contains merged alleles from ReferenceConfidenceVariantContextMerger
        ReducibleAnnotationData combinedData = new ReducibleAnnotationData(null);

        for (final ReducibleAnnotationData currentValue : annotationList) {
            parseRawDataString(currentValue);
            combineAttributeMap(currentValue, combinedData);

        }
        final Map<String, Object> annotations = new HashMap<>();
        String annotationString = makeRawAnnotationString(vcAlleles, combinedData.getAttributeMap());
        annotations.put(getRawKeyName(), annotationString);
        return annotations;
    }

    public String makeRawAnnotationString(final List<Allele> vcAlleles, final Map<Allele, Number> perAlleleData) {
        return String.format("%.2f", perAlleleData.get(Allele.NO_CALL));
    }

    @Override
    @SuppressWarnings({"unchecked", "rawtypes"})//FIXME generics here blow up
    public Map<String, Object> finalizeRawData(final VariantContext vc, final VariantContext originalVC) {
        if (!vc.hasAttribute(getRawKeyName()))
            return new HashMap<>();
        String rawMQdata = vc.getAttributeAsString(getRawKeyName(),null);
        if (rawMQdata == null)
            return new HashMap<>();

        ReducibleAnnotationData myData = new ReducibleAnnotationData(rawMQdata);
        parseRawDataString(myData);

        String annotationString = makeFinalizedAnnotationString(getNumOfReads(vc, null), myData.getAttributeMap());
        return Collections.singletonMap(getKeyNames().get(0), annotationString);
    }

    public String makeFinalizedAnnotationString(final int numOfReads, final Map<Allele, Double> perAlleleData) {
        return String.format("%.2f", Math.sqrt(perAlleleData.get(Allele.NO_CALL)/numOfReads));
    }


    public void combineAttributeMap(ReducibleAnnotationData<Double> toAdd, ReducibleAnnotationData<Double> combined) {
        if (combined.getAttribute(Allele.NO_CALL) != null) {
            combined.putAttribute(Allele.NO_CALL, combined.getAttribute(Allele.NO_CALL) + toAdd.getAttribute(Allele.NO_CALL));
        } else {
            combined.putAttribute(Allele.NO_CALL, toAdd.getAttribute(Allele.NO_CALL));
        }
    }

    @SuppressWarnings({"unchecked", "rawtypes"})//FIXME
    public void calculateRawData(final VariantContext vc,
                                 final ReadLikelihoods<Allele> likelihoods,
                                 final ReducibleAnnotationData rawAnnotations){
        //put this as a double, like GATK3.5
        final double squareSum = IntStream.range(0, likelihoods.numberOfSamples()).boxed()
                .flatMap(s -> likelihoods.sampleReads(s).stream())
                .map(GATKRead::getMappingQuality)
                .filter(mq -> mq != QualityUtils.MAPPING_QUALITY_UNAVAILABLE)
                .mapToDouble(mq -> mq * mq).sum();

        rawAnnotations.putAttribute(Allele.NO_CALL, squareSum);
    }

    @Override
    public Map<String, Object> annotate(final ReferenceContext ref,
                                        final VariantContext vc,
                                        final ReadLikelihoods<Allele> likelihoods) {
        Utils.nonNull(vc);
        if (likelihoods == null || likelihoods.readCount() < 1 ) {
            return new HashMap<>();
        }

        final Map<String, Object> annotations = new HashMap<>();
        final ReducibleAnnotationData<Double> myData = new ReducibleAnnotationData<>(null);
        calculateRawData(vc, likelihoods, myData);
        final String annotationString = makeFinalizedAnnotationString(getNumOfReads(vc, likelihoods), myData.getAttributeMap());
        annotations.put(getKeyNames().get(0), annotationString);
        return annotations;
    }

    @VisibleForTesting
    static String formattedValue(double rms) {
        return String.format("%.2f", rms);
    }

    /**
     * converts {@link GATKVCFConstants#RAW_RMS_MAPPING_QUALITY_KEY} into  {@link VCFConstants#RMS_MAPPING_QUALITY_KEY}  annotation if present
     * @param vc which potentially contains rawMQ
     * @return if vc contained {@link GATKVCFConstants#RAW_RMS_MAPPING_QUALITY_KEY} it will be replaced with {@link VCFConstants#RMS_MAPPING_QUALITY_KEY}
     * otherwise return the original vc
     */
    public VariantContext finalizeRawMQ(final VariantContext vc) {
        final String rawMQdata = vc.getAttributeAsString(getRawKeyName(), null);
        if (rawMQdata == null) {
            return vc;
        } else {
            final double squareSum = parseRawDataString(rawMQdata);
            final int numOfReads = getNumOfReads(vc, null);
            final double rms = Math.sqrt(squareSum / (double)numOfReads);
            final String finalizedRMSMAppingQuality = formattedValue(rms);
            return new VariantContextBuilder(vc)
                    .rmAttribute(getRawKeyName())
                    .attribute(getKeyNames().get(0), finalizedRMSMAppingQuality)
                    .make();
        }
    }

    protected void parseRawDataString(ReducibleAnnotationData<Number> myData) {
        final String rawDataString = myData.getRawData();
        String[] rawMQdataAsStringVector;
        rawMQdataAsStringVector = rawDataString.split(",");
        double squareSum = Double.parseDouble(rawMQdataAsStringVector[0]);
        myData.putAttribute(Allele.NO_CALL, squareSum);
    }

    //TODO once the AS annotations have been added genotype gvcfs this can be removed for a more generic approach
    private static double parseRawDataString(String rawDataString) {
        try {
            /*
             * TODO: this is copied from gatk3 where it ignored all but the first value, we should figure out if this is
             * the right thing to do or if it should just convert the string without trying to split it and fail if
             * there is more than one value
             */
            final double squareSum = Double.parseDouble(rawDataString.split(",")[0]);
            return squareSum;
        } catch (final NumberFormatException e) {
            throw new UserException.BadInput("malformed " + GATKVCFConstants.RAW_RMS_MAPPING_QUALITY_KEY + " annotation: " + rawDataString);
        }
    }


    /**
     *
     * @return the number of reads at the given site, calculated as InfoField {@link VCFConstants#DEPTH_KEY} minus the
     * format field {@link GATKVCFConstants#MIN_DP_FORMAT_KEY} or DP of each of the HomRef genotypes at that site.
     * Will fall back to calculating the reads from the stratifiedContexts then AlleleLikelyhoods data if provided.
     * @throws UserException.BadInput if the {@link VCFConstants#DEPTH_KEY} is missing or if the calculated depth is <= 0
     */
    @VisibleForTesting
    static int getNumOfReads(final VariantContext vc,
                             final ReadLikelihoods<Allele> likelihoods) {
        int numOfReads = 0;
        if (vc.hasAttribute(VCFConstants.DEPTH_KEY)) {
            numOfReads = vc.getAttributeAsInt(VCFConstants.DEPTH_KEY, -1);
            if(vc.hasGenotypes()) {
                for(final Genotype gt : vc.getGenotypes()) {
                    if(gt.isHomRef()) {
                        //site-level DP contribution will come from MIN_DP for gVCF-called reference variants or DP for BP resolution
                        if (gt.hasExtendedAttribute(GATKVCFConstants.MIN_DP_FORMAT_KEY)) {
                            numOfReads -= Integer.parseInt(gt.getExtendedAttribute(GATKVCFConstants.MIN_DP_FORMAT_KEY).toString());
                        } else if (gt.hasDP()) {
                            numOfReads -= gt.getDP();
                        }
                    }
                }
            }

        // If there is no depth key, try to compute from the likelihoods
        } else if (likelihoods != null && likelihoods.numberOfAlleles() != 0) {
            for (int i = 0; i < likelihoods.numberOfSamples(); i++) {
                for (GATKRead read : likelihoods.sampleReads(i)) {
                    if (read.getMappingQuality() != QualityUtils.MAPPING_QUALITY_UNAVAILABLE) {
                        numOfReads++;
                    }
                }
            }
        }
        if (numOfReads <= 0) {
            numOfReads = -1;  //return -1 to result in a NaN
        }
        return numOfReads;
    }

    public static RMSMappingQuality getInstance() {
        return instance;
    }
}
