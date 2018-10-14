package org.broadinstitute.hellbender.tools.walkers.annotator;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.ReducibleAnnotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.ReducibleAnnotationData;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.logging.OneShotLogger;
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
 * <p>The raw data format for this annotation consists of a list of two entries: the sum of the squared mapping qualities and the number of reads across variant (not homRef) genotypes</p>
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
    private static final OneShotLogger logger = new OneShotLogger(RMSMappingQuality.class);
    private static final RMSMappingQuality instance = new RMSMappingQuality();
    public static final int NUM_LIST_ENTRIES = 2;
    public static final int SUM_OF_SQUARES_INDEX = 0;
    public static final int TOTAL_DEPTH_INDEX = 1;

    @Override
    public String getRawKeyName() { return GATKVCFConstants.RAW_MAPPING_QUALITY_WITH_DEPTH_KEY;}   //new key for the two-value MQ data to prevent version mismatch catastrophes

    public static String getDeprecatedRawKeyName() { return GATKVCFConstants.RAW_RMS_MAPPING_QUALITY_KEY;}   //new key for the two-value MQ data to prevent version mismatch catastrophes

    @Override
    public List<String> getKeyNames() {
        return Arrays.asList(VCFConstants.RMS_MAPPING_QUALITY_KEY, getRawKeyName());
    }

    @Override
    public List<VCFInfoHeaderLine> getDescriptions() {
        return Arrays.asList(VCFStandardHeaderLines.getInfoLine(getKeyNames().get(0)));
    }

    @Override
    public List<VCFInfoHeaderLine> getRawDescriptions() {
        return Arrays.asList(GATKVCFHeaderLines.getInfoLine(getRawKeyName()));
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
        final ReducibleAnnotationData<List<Integer>> myData = new ReducibleAnnotationData<>(null);
        calculateRawData(vc, likelihoods, myData);
        final String annotationString = makeRawAnnotationString(vc.getAlleles(), myData.getAttributeMap());
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

    public String makeRawAnnotationString(final List<Allele> vcAlleles, final Map<Allele, List<Integer>> perAlleleData) {
        return String.format("%d,%d", perAlleleData.get(Allele.NO_CALL).get(0), perAlleleData.get(Allele.NO_CALL).get(1));
    }

    @Override
    @SuppressWarnings({"unchecked", "rawtypes"})//FIXME generics here blow up
    public Map<String, Object> finalizeRawData(final VariantContext vc, final VariantContext originalVC) {
        String rawMQdata;
        if (vc.hasAttribute(getRawKeyName())) {
            rawMQdata = vc.getAttributeAsString(getRawKeyName(), null);
        }
        else if (vc.hasAttribute(getDeprecatedRawKeyName())) {
            rawMQdata = vc.getAttributeAsString(getDeprecatedRawKeyName(), null);
            //the original version of ReblockGVCF produces a different MQ format -- try to handle that gracefully here just in case those files go through GenotypeGVCFs
            if (vc.hasAttribute("MQ_DP")) {
                logger.warn("Presence of MQ_DP key indicates that this tool may be running on an older output of ReblockGVCF " +
                        "that may not have compatible annotations with this GATK version. Attempting to reformat MQ data.");
                final String rawMQdepth = vc.getAttributeAsString("MQ_DP",null);
                if (rawMQdepth == null) {
                    throw new UserException.BadInput("MQ annotation data is not properly formatted. This version expects an " +
                            "int tuple of sum of squared MQ values and total reads over variant genotypes.");
                }
                rawMQdata = Math.round(Double.parseDouble(rawMQdata)) + "," + rawMQdepth;  //deprecated format was double so it needs to be converted to int
            }
            else {
                logger.warn("MQ annotation data is not properly formatted. This GATK version expects key "
                        + getRawKeyName() + " with an int tuple of sum of squared MQ values and total reads over variant "
                        + "genotypes as the value. Attempting to use deprecated MQ calculation.");
                final int numOfReads = getNumOfReads(vc, null);
                rawMQdata = Math.round(Double.parseDouble(rawMQdata)) + "," + numOfReads;   //deprecated format was double so it needs to be converted to int
            }
        }
        else {
            return new HashMap<>();
        }
        if (rawMQdata == null) {
            return new HashMap<>();
        }

        ReducibleAnnotationData<List<Integer>> myData = new ReducibleAnnotationData(rawMQdata);
        parseRawDataString(myData);

        final String annotationString = makeFinalizedAnnotationString(myData.getAttribute(Allele.NO_CALL).get(TOTAL_DEPTH_INDEX), myData.getAttribute(Allele.NO_CALL).get(SUM_OF_SQUARES_INDEX));
        return Collections.singletonMap(getKeyNames().get(0), annotationString);
    }

    public String makeFinalizedAnnotationString(final int numOfReads, final int sumOfSquaredMQs) {
        return String.format("%.2f", Math.sqrt(sumOfSquaredMQs/(double)numOfReads));
    }


    /**
     * Combine the int tuples: sum of squared mapping quality and read counts over variant genotypes
     * Since this is not an allele-specific annotation, store the result with the NO_CALL allele key.
     * @param toAdd new data
     * @param combined passed in as previous data, modified to store the combined result
     */
    public void combineAttributeMap(ReducibleAnnotationData<List<Integer>> toAdd, ReducibleAnnotationData<List<Integer>> combined) {
        if (combined.getAttribute(Allele.NO_CALL) != null) {
            combined.putAttribute(Allele.NO_CALL, Arrays.asList(combined.getAttribute(Allele.NO_CALL).get(0) + toAdd.getAttribute(Allele.NO_CALL).get(0),
                    combined.getAttribute(Allele.NO_CALL).get(1) + toAdd.getAttribute(Allele.NO_CALL).get(1)));
        } else {
            combined.putAttribute(Allele.NO_CALL, toAdd.getAttribute(Allele.NO_CALL));
        }
    }

    @SuppressWarnings({"unchecked", "rawtypes"})//FIXME
    public void calculateRawData(final VariantContext vc,
                                 final ReadLikelihoods<Allele> likelihoods,
                                 final ReducibleAnnotationData rawAnnotations){
        //GATK3.5 had a double, but change this to an int for the tuple representation (square sum, read count)
        int squareSum = 0;
        int numReadsUsed = 0;
        for (int i = 0; i < likelihoods.numberOfSamples(); i++) {
            for (final GATKRead read : likelihoods.sampleReads(i)) {
                int mq = read.getMappingQuality();
                if (mq != QualityUtils.MAPPING_QUALITY_UNAVAILABLE) {
                    squareSum += mq * mq;
                    numReadsUsed++;
                }
            }
        }
        rawAnnotations.putAttribute(Allele.NO_CALL, Arrays.asList(squareSum, numReadsUsed));
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
        final ReducibleAnnotationData<List<Integer>> myData = new ReducibleAnnotationData<>(null);
        calculateRawData(vc, likelihoods, myData);
        final String annotationString = makeFinalizedAnnotationString(myData.getAttribute(Allele.NO_CALL).get(TOTAL_DEPTH_INDEX), myData.getAttribute(Allele.NO_CALL).get(SUM_OF_SQUARES_INDEX));
        annotations.put(getKeyNames().get(0), annotationString);
        return annotations;
    }

    @VisibleForTesting
    static String formattedValue(double rms) {
        return String.format("%.2f", rms);
    }

    /**
     * converts {@link GATKVCFConstants#RAW_MAPPING_QUALITY_WITH_DEPTH_KEY} into  {@link VCFConstants#RMS_MAPPING_QUALITY_KEY}  annotation if present
     * NOTE: this is currently only used by HaplotypeCaller in VCF mode
     * @param vc which potentially contains rawMQ
     * @return if vc contained {@link GATKVCFConstants#RAW_MAPPING_QUALITY_WITH_DEPTH_KEY} it will be replaced with {@link VCFConstants#RMS_MAPPING_QUALITY_KEY}
     * otherwise return the original vc
     */
    public VariantContext finalizeRawMQ(final VariantContext vc) {
        final String rawMQdata = vc.getAttributeAsString(getRawKeyName(), null);
        if (rawMQdata == null) {
            return vc;
        } else {
            final List<Integer> SSQMQandDP = parseRawDataString(rawMQdata);
            final double rms = Math.sqrt(SSQMQandDP.get(SUM_OF_SQUARES_INDEX) / (double)SSQMQandDP.get(TOTAL_DEPTH_INDEX));
            final String finalizedRMSMAppingQuality = formattedValue(rms);
            return new VariantContextBuilder(vc)
                    .rmAttribute(getRawKeyName())
                    .attribute(getKeyNames().get(0), finalizedRMSMAppingQuality)
                    .make();
        }
    }

    protected void parseRawDataString(ReducibleAnnotationData<List<Integer>> myData) {
        myData.putAttribute(Allele.NO_CALL, parseRawDataString(myData.getRawData()));
    }

    //TODO once the AS annotations have been added genotype gvcfs this can be removed for a more generic approach
    private static List<Integer> parseRawDataString(String rawDataString) {
        try {
            final String[] parsed = rawDataString.split(",");
            if (parsed.length != NUM_LIST_ENTRIES) {
                throw new UserException.BadInput("Raw value for annotation has " + parsed.length + " values, expected " + NUM_LIST_ENTRIES);
            }
            final int squareSum = Integer.parseInt(parsed[SUM_OF_SQUARES_INDEX]);
            final int totalDP = Integer.parseInt(parsed[TOTAL_DEPTH_INDEX]);
            return Arrays.asList(squareSum,totalDP);
        } catch (final NumberFormatException e) {
            throw new UserException.BadInput("malformed " + GATKVCFConstants.RAW_RMS_MAPPING_QUALITY_KEY + " annotation: " + rawDataString);
        }
    }


    /**
     *
     * @return the number of reads at the given site, trying first {@Link GATKVCFConstants.RAW_MAPPING_QUALITY_WITH_DEPTH_KEY},
     * falling back to calculating the value as InfoField {@link VCFConstants#DEPTH_KEY} minus the
     * format field {@link GATKVCFConstants#MIN_DP_FORMAT_KEY} or DP of each of the HomRef genotypes at that site.
     * If neither of those is possible, will fall back to calculating the reads from the ReadLikelihoods data if provided.
     * @throws UserException.BadInput if the {@link VCFConstants#DEPTH_KEY} is missing or if the calculated depth is <= 0
     */
    @VisibleForTesting
    static int getNumOfReads(final VariantContext vc,
                             final ReadLikelihoods<Allele> likelihoods) {
        if(vc.hasAttribute(GATKVCFConstants.RAW_MAPPING_QUALITY_WITH_DEPTH_KEY)) {
            List<Integer> MQtuple = vc.getAttributeAsIntList(GATKVCFConstants.RAW_MAPPING_QUALITY_WITH_DEPTH_KEY,0);
            if (MQtuple.get(TOTAL_DEPTH_INDEX) > 0) {
                return MQtuple.get(TOTAL_DEPTH_INDEX);
            }
        }

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
