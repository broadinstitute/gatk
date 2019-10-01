package org.broadinstitute.hellbender.tools.walkers.annotator;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.ReducibleAnnotation;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.logging.OneShotLogger;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;

import java.util.*;

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


public final class RMSMappingQuality extends InfoFieldAnnotation implements StandardAnnotation, ReducibleAnnotation<RMSMappingQuality.Raw> {


    static class Raw extends AbstractList<Long> {
        private long sumOfSquares;
        private long depth;

        private Raw(final long sumOfSquares, final long depth) {
            this.sumOfSquares = sumOfSquares;
            this.depth = depth;
        }

        private void add(final Raw raw) {
            sumOfSquares += raw.sumOfSquares;
            depth += raw.depth;
        }

        @Override
        public int size() {
            return 2;
        }

        public String toString() {
            return sumOfSquares + "," + depth;
        }

        public static Raw parse(final String str) {
            final int commaIndex = str.indexOf(',');
            if (commaIndex == -1) {
                throw new IllegalArgumentException("bad format");
            } else {
                final long sumOfSquares = parseStrictPositiveLong(str, 0, commaIndex);
                final long depth = parseStrictPositiveLong(str, commaIndex + 1, str.length());
                return new Raw(sumOfSquares, depth);
            }
        }

        private static long parseStrictPositiveLong(final String str, final int start, final int end) {
            int idx;
            char ch = 0;
            for(idx = start; idx < end; idx++) {
                ch = str.charAt(idx);
                if (ch != ' ') {
                    break;
                }
            }
            if (idx == end) {
                throw new IllegalArgumentException("bad long: " + str.substring(start, end));
            } else {
                long total = 0;
                while (ch >= '0' && ch <= '9' && idx < end) {
                    total = total * 10 + (ch - '0');
                    if (++idx == end) {
                        break;
                    } else {
                        ch = str.charAt(idx);
                    }
                }
                while (idx < end) {
                    if (ch != ' ') {
                        throw new IllegalArgumentException("bad long: " + str.substring(start, end));
                    } else if (++idx == end) {
                        break;
                    } else {
                        ch = str.charAt(idx);
                    }
                }
                return total;
            }
        }

        @Override
        public Long get(final int index) {
            if (index == 0) {
                return sumOfSquares;
            } else if (index == 1) {
                return depth;
            } else {
                throw new IndexOutOfBoundsException("" + index);
            }
        }

        public double toRMS() {
            return sumOfSquares / (double) depth;
        }
    }

    private static final OneShotLogger logger = new OneShotLogger(RMSMappingQuality.class);
    private static final RMSMappingQuality INSTANCE = new RMSMappingQuality();
    @SuppressWarnings("unused")
    private static final int SUM_OF_SQUARES_INDEX = 0;
    private static final int TOTAL_DEPTH_INDEX = 1;
    private static final String OUTPUT_PRECISION = "%.2f";
    public static final String RMS_MAPPING_QUALITY_OLD_BEHAVIOR_OVERRIDE_ARGUMENT = "allow-old-rms-mapping-quality-annotation-data";

    private RMSMappingQuality() {};

    @Argument(fullName = RMS_MAPPING_QUALITY_OLD_BEHAVIOR_OVERRIDE_ARGUMENT, doc="Override to allow old RMSMappingQuality annotated VCFs to function", optional=true)
    public boolean allowOlderRawKeyValues = false;

    @Override
    public String getKeyName() {
        return VCFConstants.RMS_MAPPING_QUALITY_KEY;
    }

    @Override
    public String getRawKeyName() { return GATKVCFConstants.RAW_MAPPING_QUALITY_WITH_DEPTH_KEY;}   //new key for the two-value MQ data to prevent version mismatch catastrophes

    @Override
    public Raw computeRaw(final VariantContextBuilder builder, final ReferenceContext ref, final VariantContext vc, final ReadLikelihoods<Allele> likelihoods) {
        Utils.nonNull(vc);
        if (likelihoods == null || likelihoods.readCount() < 1) {
            return null;
        } else {
            final Raw result = new Raw(0, 0);
            for (int i = 0; i < likelihoods.numberOfSamples(); i++) {
                for (final GATKRead read : likelihoods.sampleReads(i)) {
                    long mq = read.getMappingQuality();
                    if (mq != QualityUtils.MAPPING_QUALITY_UNAVAILABLE) {
                        result.sumOfSquares += mq * mq;
                        result.depth++;
                    }
                }
            }
            return result;
        }
    }

    @Override
    public Raw reduceRaws(final VariantContextBuilder builder, final List<Raw> values) {
        final int length = values.size();
        if (length == 1) {
            return values.get(0);
        } else if (length == 0) {
            return null;
        } else {
            final Raw result = values.get(0);
            for (int i = 1; i < length; i++) {
                result.add(values.get(i));
            }
            return result;
        }
    }

    @Override
    public void setRaw(final VariantContextBuilder builder, final Raw value) {
        builder.attribute(getRawKeyName(), String.format("%d,%d", value.sumOfSquares, value.depth));
    }

    private Raw parseRaw(final String asString) {
        if (asString == null) {
            return null;
        } else {
            final int commaIndex = asString.indexOf(',');
            if (commaIndex == -1) {
                throw new GATKException("bad annotation value: " + asString);
            } else {
                try {
                    final Raw result = new Raw(0, 0);
                    result.sumOfSquares = Long.parseLong(asString.substring(0, commaIndex).trim());
                    result.depth = Long.parseLong(asString.substring(commaIndex + 1).trim());
                    return result;
                } catch (final NumberFormatException ex) {
                    throw new GATKException("bad annotation value: " + asString, ex);
                }
            }
        }
    }

    public Raw getRaw(final VariantContext context) {
        if (context.hasAttribute(getRawKeyName())) {
            return parseRaw(context.getAttributeAsString(getRawKeyName(), null));
        } else if (context.hasAttribute(getDeprecatedRawKeyName())) {
            if (!allowOlderRawKeyValues) {
                throw new UserException.BadInput("Presence of '-"+getDeprecatedRawKeyName()+"' annotation is detected. This GATK version expects key "
                        + getRawKeyName() + " with a tuple of sum of squared MQ values and total reads over variant "
                        + "genotypes as the value. This could indicate that the provided input was produced with an older version of GATK. " +
                        "Use the argument '--"+RMS_MAPPING_QUALITY_OLD_BEHAVIOR_OVERRIDE_ARGUMENT+"' to override and attempt the deprecated MQ calculation. There " +
                        "may be differences in how newer GATK versions calculate DP and MQ that may result in worse MQ results. Use at your own risk.");
            } else if (context.hasAttribute(GATKVCFConstants.MAPPING_QUALITY_DEPTH_DEPRECATED)) {
                logger.warn("Presence of " + GATKVCFConstants.MAPPING_QUALITY_DEPTH_DEPRECATED + " key indicates that this tool may be running on an older output of ReblockGVCF " +
                        "that may not have compatible annotations with this GATK version. Attempting to reformat MQ data.");
                final String rawMQdepth = context.getAttributeAsString(GATKVCFConstants.MAPPING_QUALITY_DEPTH_DEPRECATED,null);
                if (rawMQdepth == null) {
                    throw new UserException.BadInput("MQ annotation data is not properly formatted. This version expects a " +
                            "long tuple of sum of squared MQ values and total reads over variant genotypes.");
                }

                return parseRaw(Math.round(context.getAttributeAsDouble(getDeprecatedRawKeyName(), Double.NaN)) + "," + rawMQdepth);
            } else {
                logger.warn("MQ annotation data is not properly formatted. This GATK version expects key "
                        + getRawKeyName() + " with a tuple of sum of squared MQ values and total reads over variant "
                        + "genotypes as the value. Attempting to use deprecated MQ calculation.");
                final long numOfReads = getNumOfReads(context, null);
                return parseRaw(Math.round(context.getAttributeAsDouble(getDeprecatedRawKeyName(), Double.NaN)) + "," + numOfReads);
            }
        } else {
            return null;
        }
    }

    @Override
    public Object toFinalizedAnnotation(final Raw raw) {
        return raw == null ? null : raw.toRMS();
    }


    public static String getDeprecatedRawKeyName() { return GATKVCFConstants.RAW_RMS_MAPPING_QUALITY_KEY;}   //new key for the two-value MQ data to prevent version mismatch catastrophes

    @Override
    public List<String> getKeyNames() {
        return Arrays.asList(getKeyName(), getRawKeyName());
    }

    @Override
    public List<VCFInfoHeaderLine> getDescriptions() {
        return Arrays.asList(VCFStandardHeaderLines.getInfoLine(getKeyName()));
    }

    @Override
    public List<VCFInfoHeaderLine> getRawDescriptions() {
        return Arrays.asList(GATKVCFHeaderLines.getInfoLine(getRawKeyName()));
    }

    public void annotateUsingRaw(final VariantContextBuilder builder, final Raw raw) {
        builder.attribute(getKeyName(), raw.toRMS());
    }

    public void annotateUsingRawIfNeeded(final VariantContextBuilder builder, final VariantContext context) {
        if (context.hasAttribute(getKeyName())) {
            builder.attribute(getKeyName(), context.getAttribute(getKeyName()));
        } else {
            final Raw raw = getRaw(context);
            if (raw != null) {
                builder.attribute(getKeyName(), raw.toRMS());
            }
        }
    }


    @Override
    public Map<String, Object> annotate(final ReferenceContext ref,
                                        final VariantContext vc,
                                        final ReadLikelihoods<Allele> likelihoods) {
        Utils.nonNull(vc);
        if (likelihoods == null || likelihoods.readCount() < 1 ) {
            return new HashMap<>();
        }

        final Raw raw = this.computeRaw(null, ref, vc, likelihoods);
        return Collections.singletonMap(getKeyName(), raw.toRMS());
    }

    @VisibleForTesting
    static String formattedValue(double rms) {
        return String.format(OUTPUT_PRECISION, rms);
    }

    /**
     * converts {@link GATKVCFConstants#RAW_MAPPING_QUALITY_WITH_DEPTH_KEY} into  {@link VCFConstants#RMS_MAPPING_QUALITY_KEY}  annotation if present
     * NOTE: this is currently only used by HaplotypeCaller in VCF mode and GnarlyGenotyper
     * @param vc which potentially contains rawMQ
     * @return if vc contained {@link GATKVCFConstants#RAW_MAPPING_QUALITY_WITH_DEPTH_KEY} it will be replaced with {@link VCFConstants#RMS_MAPPING_QUALITY_KEY}
     * otherwise return the original vc
     */
    public void finalizeAnnotation(final VariantContextBuilder builder, final VariantContext vc) {
        final Raw raw = getRaw(vc);
        if (raw != null) {
            builder.attribute(getKeyName(), raw.toRMS());
            builder.rmAttribute(getRawKeyName());
            builder.rmAttribute(getDeprecatedRawKeyName());
        }
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
        if (vc.hasAttribute(GATKVCFConstants.RAW_MAPPING_QUALITY_WITH_DEPTH_KEY)) {
            final List<Long> mqTuple = getAttributeAsLongList(vc, GATKVCFConstants.RAW_MAPPING_QUALITY_WITH_DEPTH_KEY, 0L);
            if (mqTuple.get(TOTAL_DEPTH_INDEX) > 0) {
                return mqTuple.get(TOTAL_DEPTH_INDEX);
            }
        }
        if (vc.hasAttribute(GATKVCFConstants.MAPPING_QUALITY_DEPTH_DEPRECATED)) {
            final int mqDP = vc.getAttributeAsInt(GATKVCFConstants.MAPPING_QUALITY_DEPTH_DEPRECATED, 0);
            if (mqDP > 0) {
                return mqDP;
            }
        }
        if (vc.hasAttribute(VCFConstants.DEPTH_KEY)) {
            long numOfReads = getAttributeAsLong(vc, VCFConstants.DEPTH_KEY, -1L);
            if (vc.hasGenotypes()) {
                for (final Genotype gt : vc.getGenotypes()) {
                    if (gt.isHomRef()) {
                        //site-level DP contribution will come from MIN_DP for gVCF-called reference variants or DP for BP resolution
                        if (gt.hasExtendedAttribute(GATKVCFConstants.MIN_DP_FORMAT_KEY)) {
                            numOfReads -= Long.parseLong(gt.getExtendedAttribute(GATKVCFConstants.MIN_DP_FORMAT_KEY).toString());
                        } else if (gt.hasDP()) {
                            numOfReads -= gt.getDP();
                        }
                    }
                }
            }
            return numOfReads <= 0 ? -1 : numOfReads;
        } else if (likelihoods != null && likelihoods.numberOfAlleles() != 0) {
            long numOfReads = 0;
            for (int i = 0; i < likelihoods.numberOfSamples(); i++) {
                for (GATKRead read : likelihoods.sampleReads(i)) {
                    if (read.getMappingQuality() != QualityUtils.MAPPING_QUALITY_UNAVAILABLE) {
                        numOfReads++;
                    }
                }
            }
            return numOfReads <= 0 ? -1 : numOfReads;
        } else {
            return -1;
        }
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
        return INSTANCE;
    }
}
