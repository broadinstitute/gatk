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
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static org.broadinstitute.hellbender.utils.GATKProtectedVariantContextUtils.getAttributeAsLong;
import static org.broadinstitute.hellbender.utils.GATKProtectedVariantContextUtils.getAttributeAsLongList;


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
    private static final int NUM_LIST_ENTRIES = 2;
    private static final int SUM_OF_SQUARES_INDEX = 0;
    private static final int TOTAL_DEPTH_INDEX = 1;

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
        final ReducibleAnnotationData<List<Long>> myData = new ReducibleAnnotationData<>(null);
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

    private String makeRawAnnotationString(final List<Allele> vcAlleles, final Map<Allele, List<Long>> perAlleleData) {
        return String.format("%d,%d", perAlleleData.get(Allele.NO_CALL).get(SUM_OF_SQUARES_INDEX), perAlleleData.get(Allele.NO_CALL).get(TOTAL_DEPTH_INDEX));
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
                    throw new UserException.BadInput("MQ annotation data is not properly formatted. This version expects a " +
                            "long tuple of sum of squared MQ values and total reads over variant genotypes.");
                }
                rawMQdata = Math.round(Double.parseDouble(rawMQdata)) + "," + rawMQdepth;  //deprecated format was double so it needs to be converted to long
            }
            else {
                logger.warn("MQ annotation data is not properly formatted. This GATK version expects key "
                        + getRawKeyName() + " with an long tuple of sum of squared MQ values and total reads over variant "
                        + "genotypes as the value. Attempting to use deprecated MQ calculation.");
                final long numOfReads = getNumOfReads(vc, null);
                rawMQdata = Math.round(Double.parseDouble(rawMQdata)) + "," + numOfReads;   //deprecated format was double so it needs to be converted to long
            }
        }
        else {
            return new HashMap<>();
        }
        if (rawMQdata == null) {
            return new HashMap<>();
        }

        ReducibleAnnotationData<List<Long>> myData = new ReducibleAnnotationData(rawMQdata);
        parseRawDataString(myData);

        final String annotationString = makeFinalizedAnnotationString(myData.getAttribute(Allele.NO_CALL).get(TOTAL_DEPTH_INDEX), myData.getAttribute(Allele.NO_CALL).get(SUM_OF_SQUARES_INDEX));
        return Collections.singletonMap(getKeyNames().get(0), annotationString);
    }

    private String makeFinalizedAnnotationString(final long numOfReads, final long sumOfSquaredMQs) {
        return String.format("%.2f", Math.sqrt(sumOfSquaredMQs/(double)numOfReads));
    }


    /**
     * Combine the long tuples: sum of squared mapping quality and read counts over variant genotypes
     * Since this is not an allele-specific annotation, store the result with the NO_CALL allele key.
     * @param toAdd new data
     * @param combined passed in as previous data, modified to store the combined result
     */
    private void combineAttributeMap(ReducibleAnnotationData<List<Long>> toAdd, ReducibleAnnotationData<List<Long>> combined) {
        if (combined.getAttribute(Allele.NO_CALL) != null) {
            combined.putAttribute(Allele.NO_CALL, Arrays.asList(combined.getAttribute(Allele.NO_CALL).get(SUM_OF_SQUARES_INDEX) + toAdd.getAttribute(Allele.NO_CALL).get(SUM_OF_SQUARES_INDEX),
                    combined.getAttribute(Allele.NO_CALL).get(TOTAL_DEPTH_INDEX) + toAdd.getAttribute(Allele.NO_CALL).get(TOTAL_DEPTH_INDEX)));
        } else {
            combined.putAttribute(Allele.NO_CALL, toAdd.getAttribute(Allele.NO_CALL));
        }
    }

    @SuppressWarnings({"unchecked", "rawtypes"})//FIXME
    private void calculateRawData(final VariantContext vc,
                                 final ReadLikelihoods<Allele> likelihoods,
                                 final ReducibleAnnotationData rawAnnotations){
        //GATK3.5 had a double, but change this to an long for the tuple representation (square sum, read count)
        long squareSum = 0;
        long numReadsUsed = 0;
        for (int i = 0; i < likelihoods.numberOfSamples(); i++) {
            for (final GATKRead read : likelihoods.sampleReads(i)) {
                long mq = read.getMappingQuality();
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
        final ReducibleAnnotationData<List<Long>> myData = new ReducibleAnnotationData<>(null);
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
     * NOTE: this is currently only used by HaplotypeCaller in VCF mode and GnarlyGenotyper
     * @param vc which potentially contains rawMQ
     * @return if vc contained {@link GATKVCFConstants#RAW_MAPPING_QUALITY_WITH_DEPTH_KEY} it will be replaced with {@link VCFConstants#RMS_MAPPING_QUALITY_KEY}
     * otherwise return the original vc
     */
    public VariantContext finalizeRawMQ(final VariantContext vc) {
        final String rawMQdata = vc.getAttributeAsString(getRawKeyName(), null);
        if (rawMQdata == null) {
            if (!vc.hasAttribute(GATKVCFConstants.MAPPING_QUALITY_DEPTH)) {
                return vc;
            }
            if (vc.hasAttribute(GATKVCFConstants.MAPPING_QUALITY_DEPTH)) {
                final int numOfReads = vc.getAttributeAsInt(GATKVCFConstants.MAPPING_QUALITY_DEPTH, getNumOfReads(vc));  //MQ_DP is an undocumented hack for the Gnarly Pipeline -- improved version uses RAW_MQ_and_DP tuple format (see #4969)
                final String deprecatedRawMQdata = vc.getAttributeAsString(getDeprecatedRawKeyName(), null);
                final double squareSum = parseDeprecatedRawDataString(deprecatedRawMQdata);
                final double rms = Math.sqrt(squareSum / (double)numOfReads);
                final String finalizedRMSMAppingQuality = formattedValue(rms);
                return new VariantContextBuilder(vc)
                        .rmAttribute(getDeprecatedRawKeyName())
                        .attribute(getKeyNames().get(0), finalizedRMSMAppingQuality)
                        .make();
            }

        } else {
            final List<Long> SSQMQandDP = parseRawDataString(rawMQdata);
            final double rms = Math.sqrt(SSQMQandDP.get(SUM_OF_SQUARES_INDEX) / (double)SSQMQandDP.get(TOTAL_DEPTH_INDEX));
            final String finalizedRMSMAppingQuality = formattedValue(rms);
            return new VariantContextBuilder(vc)
                    .rmAttribute(getRawKeyName())
                    .attribute(getKeyNames().get(0), finalizedRMSMAppingQuality)
                    .make();
        }
        return vc;
    }

    private void parseRawDataString(ReducibleAnnotationData<List<Long>> myData) {
        myData.putAttribute(Allele.NO_CALL, parseRawDataString(myData.getRawData()));
    }

    //TODO once the AS annotations have been added genotype gvcfs this can be removed for a more generic approach
    private static List<Long> parseRawDataString(String rawDataString) {
        try {
            final String[] parsed = rawDataString.trim().replaceAll("\\[|\\]", "").split(", *");
            if (parsed.length != NUM_LIST_ENTRIES) {
                throw new UserException.BadInput("Raw value for annotation has " + parsed.length + " values, expected " + NUM_LIST_ENTRIES);
            }
            final long squareSum = Long.parseLong(parsed[SUM_OF_SQUARES_INDEX]);
            final long totalDP = Long.parseLong(parsed[TOTAL_DEPTH_INDEX]);
            return Arrays.asList(squareSum,totalDP);
        } catch (final NumberFormatException e) {
            throw new UserException.BadInput("malformed " + GATKVCFConstants.RAW_MAPPING_QUALITY_WITH_DEPTH_KEY + " annotation: " + rawDataString);
        }
    }

    //Maintain some semblance of backward compatability by keeping the ability to use the old annotation key and format
    private static double parseDeprecatedRawDataString(String rawDataString) {
        try {
            /*
             * TODO: this is copied from gatk3 where it ignored all but the first value, we should figure out if this is
             * the right thing to do or if it should just convert the string without trying to split it and fail if
             * there is more than one value
             */
            final double squareSum = Double.parseDouble(rawDataString.split(",")[0]);
            return squareSum;
        } catch (final NumberFormatException e) {
            throw new UserException.BadInput("malformed " + getDeprecatedRawKeyName() + " annotation: " + rawDataString);
        }
    }

    /**
     *
     * @return the number of reads at the given site, calculated as InfoField {@link VCFConstants#DEPTH_KEY} minus the
     * format field {@link GATKVCFConstants#MIN_DP_FORMAT_KEY} or DP of each of the HomRef genotypes at that site
     * @throws UserException.BadInput if the {@link VCFConstants#DEPTH_KEY} is missing or if the calculated depth is <= 0
     */
    @VisibleForTesting
    private static int getNumOfReads(final VariantContext vc) {
        if(vc.hasAttribute(GATKVCFConstants.MAPPING_QUALITY_DEPTH)) {
            int mqDP = vc.getAttributeAsInt(GATKVCFConstants.MAPPING_QUALITY_DEPTH, 0);
            if (mqDP > 0) {
                return mqDP;
            }
        }

        //don't use the full depth because we don't calculate MQ for reference blocks
        //don't count spanning deletion calls towards number of reads
        int numOfReads = vc.getAttributeAsInt(VCFConstants.DEPTH_KEY, -1);
        if(vc.hasGenotypes()) {
            for(final Genotype gt : vc.getGenotypes()) {
               if(hasReferenceDepth(gt)) {
                    //site-level DP contribution will come from MIN_DP for gVCF-called reference variants or DP for BP resolution
                    if (gt.hasExtendedAttribute(GATKVCFConstants.MIN_DP_FORMAT_KEY)) {
                        numOfReads -= Integer.parseInt(gt.getExtendedAttribute(GATKVCFConstants.MIN_DP_FORMAT_KEY).toString());
                    } else if (gt.hasDP()) {
                        numOfReads -= gt.getDP();
                    }
                }
                else if(hasSpanningDeletionAllele(gt)) {
                    //site-level DP contribution will come from MIN_DP for gVCF-called reference variants or DP for BP resolution
                    if (gt.hasExtendedAttribute(GATKVCFConstants.MIN_DP_FORMAT_KEY)) {
                        numOfReads -= Integer.parseInt(gt.getExtendedAttribute(GATKVCFConstants.MIN_DP_FORMAT_KEY).toString());
                    } else if (gt.hasDP()) {
                        numOfReads -= gt.getDP();
                    }
                }
            }
        }
        if (numOfReads <= 0){
            numOfReads = -1;  //return -1 to result in a NaN
        }
        return numOfReads;
    }

    //In the new reducible framework only samples that get annotated at the GVCF level contribute to MQ
    //The problem is that DP includes those samples plus the min_DP of the homRef blocks, which don't contribute MQ
    //The fix is to pull out reference blocks, whether or not they have a called GT, but don't subtract depth from PL=[0,0,0] sites because they're still "variant"
    //This is still inaccurate if there's an annotated homRef in the GVCF, which does happen for really low evidence alleles, but we won't know after the samples are merged
    private static boolean hasReferenceDepth(Genotype gt) {
        return gt.isHomRef() || (gt.isNoCall() && gt.hasPL() && gt.getPL()[0] == 0 && gt.getPL()[1] != 0);
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
    static long getNumOfReads(final VariantContext vc,
                             final ReadLikelihoods<Allele> likelihoods) {
        if(vc.hasAttribute(GATKVCFConstants.RAW_MAPPING_QUALITY_WITH_DEPTH_KEY)) {
            List<Long> MQtuple = getAttributeAsLongList(vc, GATKVCFConstants.RAW_MAPPING_QUALITY_WITH_DEPTH_KEY,0L);
            if (MQtuple.get(TOTAL_DEPTH_INDEX) > 0) {
                return MQtuple.get(TOTAL_DEPTH_INDEX);
            }
        }

        long numOfReads = 0;
        if (vc.hasAttribute(VCFConstants.DEPTH_KEY)) {
            numOfReads = getAttributeAsLong(vc, VCFConstants.DEPTH_KEY, -1L);
            if(vc.hasGenotypes()) {
                for(final Genotype gt : vc.getGenotypes()) {
                    if(gt.isHomRef()) {
                        //site-level DP contribution will come from MIN_DP for gVCF-called reference variants or DP for BP resolution
                        if (gt.hasExtendedAttribute(GATKVCFConstants.MIN_DP_FORMAT_KEY)) {
                            numOfReads -= Long.parseLong(gt.getExtendedAttribute(GATKVCFConstants.MIN_DP_FORMAT_KEY).toString());
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

    private static boolean hasSpanningDeletionAllele(final Genotype gt) {
        for(final Allele a : gt.getAlleles()) {
            boolean hasSpanningDeletion = GATKVCFConstants.isSpanningDeletion(a);
            if(hasSpanningDeletion) {
                return true;
            }
        }
        return false;
    }

    public static RMSMappingQuality getInstance() {
        return instance;
    }
}
