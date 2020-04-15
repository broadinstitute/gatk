package org.broadinstitute.hellbender.tools.copynumber.gcnv;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.FeatureReader;
import htsjdk.tribble.TribbleException;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.apache.commons.math3.util.FastMath;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.copynumber.GermlineCNVCaller;
import org.broadinstitute.hellbender.tools.copynumber.PostprocessGermlineCNVCalls;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.IntegerCopyNumberSegment;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFHeaderLines;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.reference.ReferenceUtils;
import org.broadinstitute.hellbender.utils.variant.GATKSVVariantContextUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Helper class for {@link PostprocessGermlineCNVCalls} for single-sample postprocessing of segmented
 * {@link GermlineCNVCaller} calls.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class GermlineCNVSegmentVariantComposer extends GermlineCNVVariantComposer<IntegerCopyNumberSegment> {

    /* VCF FORMAT header keys */

    /**
     * Number of points in the segment
     */
    public static final String NP = "NP";

    /**
     * Quality metric (some points called)
     */
    public static final String QS = "QS";

    /**
     * Quality metric (all points called)
     */
    public static final String QA = "QA";

    /**
     * Quality metric (segment start)
     */
    public static final String QSS = "QSS";

    /**
     * Quality metric (segment end)
     */
    public static final String QSE = "QSE";

    public static final String SCHEMA_HEADER_KEY = "gcnvVcfSchemaVersion";
    public static final String CURRENT_SCHEMA_VERSION = "2.0";

    protected final Logger logger = LogManager.getLogger(this.getClass());

    private final IntegerCopyNumberState refAutosomalCopyNumberState;
    private final Set<String> allosomalContigSet;
    private final ReferenceSequenceFile reference;
    private final int dupeQSThreshold;
    private final int hetDelQSThreshold;
    private final int homDelQSThreshold;
    private final double siteFreqThreshold;
    private final File clusteredCohortVcf;
    private final FeatureReader<VariantContext> clusteredVCFReader;

    /**
     * Constructor.
     *
     * @param outputWriter variant context writer
     * @param sampleName sample name
     * @param refAutosomalCopyNumberState ref copy-number state on autosomal contigs
     * @param allosomalContigSet set of allosomal contigs (ref copy-number allele be chosen according to
     *                           given contig baseline copy-number states)
     * @param reference may be null
     * @param dupeQSThreshold QS filtering threshold for duplications
     * @param hetDelQSThreshold QS filtering threshold for heteozygous deletions
     * @param homDelQSThreshold QS filtering threshold for homozygous deletions
     * @param siteFreqThreshold site frequency to eliminate common events
     * @param clusteredCohortVcf VCF describing unified segments across the whole cohort
     *
     */
    public GermlineCNVSegmentVariantComposer(final VariantContextWriter outputWriter,
                                             final String sampleName,
                                             final IntegerCopyNumberState refAutosomalCopyNumberState,
                                             final Set<String> allosomalContigSet,
                                             final ReferenceSequenceFile reference,
                                             final int dupeQSThreshold,
                                             final int hetDelQSThreshold,
                                             final int homDelQSThreshold,
                                             final double siteFreqThreshold,
                                             final File clusteredCohortVcf) {
        super(outputWriter, sampleName);
        this.refAutosomalCopyNumberState = Utils.nonNull(refAutosomalCopyNumberState);
        this.allosomalContigSet = Utils.nonNull(allosomalContigSet);
        this.reference = reference;
        this.dupeQSThreshold = dupeQSThreshold;
        this.hetDelQSThreshold = hetDelQSThreshold;
        this.homDelQSThreshold = homDelQSThreshold;
        this.siteFreqThreshold = siteFreqThreshold;
        this.clusteredCohortVcf = clusteredCohortVcf;
        if (clusteredCohortVcf != null) {
            try {
                clusteredVCFReader = AbstractFeatureReader.getFeatureReader(clusteredCohortVcf.getAbsolutePath(), new VCFCodec());
            } catch ( final TribbleException ex ) {
                throw new GATKException("Error - IO problem with file " + clusteredCohortVcf, ex);
            }
        } else {
            clusteredVCFReader = null;
        }
    }

    public GermlineCNVSegmentVariantComposer(final VariantContextWriter outputWriter,
                                             final String sampleName,
                                             final IntegerCopyNumberState refAutosomalCopyNumberState,
                                             final Set<String> allosomalContigSet,
                                             final ReferenceSequenceFile reference) {
        this(outputWriter, sampleName, refAutosomalCopyNumberState, allosomalContigSet, reference, 0, 0, 0, 0, null);
    }

    @Override
    public void composeVariantContextHeader(final SAMSequenceDictionary sequenceDictionary,
                                            final Set<VCFHeaderLine> vcfDefaultToolHeaderLines) {
        final VCFHeader result = new VCFHeader(Collections.emptySet(), Collections.singletonList(sampleName));

        /* add VCF version */
        result.addMetaDataLine(new VCFHeaderLine(VCFHeaderVersion.VCF4_2.getFormatString(),
                VCFHeaderVersion.VCF4_2.getVersionString()));

        /* add default tool header lines */
        vcfDefaultToolHeaderLines.forEach(result::addMetaDataLine);

        // add tool output schema version
        result.addMetaDataLine(new VCFHeaderLine(SCHEMA_HEADER_KEY, CURRENT_SCHEMA_VERSION));

        result.setSequenceDictionary(sequenceDictionary);

        /* header lines for annotations copied from cohort VCF */
        result.addMetaDataLine(GATKSVVCFHeaderLines.getInfoLine(GATKSVVCFConstants.SVTYPE));
        result.addMetaDataLine(GATKSVVCFHeaderLines.getInfoLine(GATKSVVCFConstants.SVLEN));
        result.addMetaDataLine(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.ORIGINAL_AC_KEY));
        result.addMetaDataLine(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.ORIGINAL_AF_KEY));
        result.addMetaDataLine(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.ORIGINAL_AN_KEY));


        /* header lines related to genotype formatting */
        result.addMetaDataLine(new VCFFormatHeaderLine(VCFConstants.GENOTYPE_KEY, 1,
                VCFHeaderLineType.String, "Segment genotype"));
        result.addMetaDataLine(new VCFFormatHeaderLine(GATKSVVCFConstants.COPY_NUMBER_FORMAT, 1,
                VCFHeaderLineType.Integer, "Segment most-likely copy-number call"));
        result.addMetaDataLine(new VCFFormatHeaderLine(NP, 1,
                VCFHeaderLineType.Integer, "Number of points (i.e. targets or bins) in the segment"));
        result.addMetaDataLine(new VCFFormatHeaderLine(QS, 1,
                VCFHeaderLineType.Integer, "Complementary Phred-scaled probability that at least one point " +
                "(i.e. target or bin) in the segment agrees with the segment copy-number call"));
        result.addMetaDataLine(new VCFFormatHeaderLine(QA, 1,
                VCFHeaderLineType.Integer, "Complementary Phred-scaled probability that all points " +
                "(i.e. targets or bins) in the segment agree with the segment copy-number call"));
        result.addMetaDataLine(new VCFFormatHeaderLine(QSS, 1,
                VCFHeaderLineType.Integer, "Complementary Phred-scaled probability that the segment start " +
                "position is a genuine copy-number changepoint"));
        result.addMetaDataLine(new VCFFormatHeaderLine(QSE, 1,
                VCFHeaderLineType.Integer, "Complementary Phred-scaled probability that the segment end " +
                "position is a genuine copy-number changepoint"));

        /* INFO header lines */
        result.addMetaDataLine(new VCFInfoHeaderLine(VCFConstants.END_KEY, 1,
                VCFHeaderLineType.Integer, "End coordinate of the variant"));

        /*FILTER header lines */
        result.addMetaDataLine(GATKSVVCFHeaderLines.getFilterLine(GATKSVVCFConstants.LOW_QS_SCORE_FILTER_KEY));
        result.addMetaDataLine(GATKSVVCFHeaderLines.getFilterLine(GATKSVVCFConstants.FREQUENCY_FILTER_KEY));

        outputWriter.writeHeader(result);
    }

    /**
     * Compose a variant context from a given {@link IntegerCopyNumberSegment}
     *
     * @param segment an instance of {@link IntegerCopyNumberSegment}
     * @return composed variant context
     */
    @VisibleForTesting
    VariantContext composeVariantContext(final IntegerCopyNumberSegment segment) {
        final String contig = segment.getContig();
        final int start = segment.getStart();
        final int end = segment.getEnd();
        final int copyNumberCall = segment.getCallIntegerCopyNumberState().getCopyNumber();
        final Allele refAllele = reference == null ? REF_ALLELE : Allele.create(ReferenceUtils.getRefBaseAtPosition(reference, contig, start), true);

        final VariantContextBuilder variantContextBuilder = new VariantContextBuilder();
        variantContextBuilder.chr(contig);
        variantContextBuilder.start(start);
        variantContextBuilder.stop(end);
        variantContextBuilder.id(String.format(VARIANT_PREFIX + "_%s_%d_%d", contig, start, end));

        final GenotypeBuilder genotypeBuilder = new GenotypeBuilder(sampleName);
        final IntegerCopyNumberState refCopyNumber = allosomalContigSet.contains(contig)
                ? segment.getBaselineIntegerCopyNumberState()
                : refAutosomalCopyNumberState;
        genotypeBuilder.alleles(GATKSVVariantContextUtils.makeGenotypeAllelesFromCopyNumber(copyNumberCall, refCopyNumber.getCopyNumber(), refAllele));
        genotypeBuilder.attribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT, copyNumberCall);
        genotypeBuilder.attribute(NP, segment.getNumPoints());
        genotypeBuilder.attribute(QS, FastMath.round(segment.getQualitySomeCalled()));
        genotypeBuilder.attribute(QA, FastMath.round(segment.getQualityAllCalled()));
        genotypeBuilder.attribute(QSS, FastMath.round(segment.getQualityStart()));
        genotypeBuilder.attribute(QSE, FastMath.round(segment.getQualityEnd()));
        //don't build yet because genotype might get updated from input cohort VCF

        final Set<Allele> uniquifiedAlleles = new HashSet<>();
        uniquifiedAlleles.add(refAllele);
        if (copyNumberCall > refCopyNumber.getCopyNumber()) {
            uniquifiedAlleles.add(GATKSVVCFConstants.DUP_ALLELE);  //dupes need additional ALTs since their genotypes are no-call
        } else if (copyNumberCall < refCopyNumber.getCopyNumber()) {
            uniquifiedAlleles.add(GATKSVVCFConstants.DEL_ALLELE);  //dels may be no-call, in which case we need to add their ALT
        }
        variantContextBuilder.alleles(uniquifiedAlleles);
        variantContextBuilder.attribute(VCFConstants.END_KEY, end);

        //copy over allele frequency etc.
        List<Double> alleleFrequency = new ArrayList<>();
        if (clusteredVCFReader != null) {
            try {
                final List<VariantContext> matches = clusteredVCFReader.query(
                        segment.getContig(), segment.getStart(), segment.getEnd()).stream().filter(vc -> vc.getStart() == segment.getStart() && vc.getEnd() == segment.getEnd()).collect(Collectors.toList());
                if (!matches.isEmpty() && !(matches.get(0) == null)) {
                    final VariantContext cohortVC = matches.get(0);
                    copyAnnotationIfPresent(cohortVC, variantContextBuilder, VCFConstants.ALLELE_COUNT_KEY, GATKVCFConstants.ORIGINAL_AC_KEY);
                    copyAnnotationIfPresent(cohortVC, variantContextBuilder, VCFConstants.ALLELE_FREQUENCY_KEY, GATKVCFConstants.ORIGINAL_AF_KEY);
                    alleleFrequency = cohortVC.getAttributeAsDoubleList(VCFConstants.ALLELE_FREQUENCY_KEY, 0.0);
                    copyAnnotationIfPresent(cohortVC, variantContextBuilder, VCFConstants.ALLELE_NUMBER_KEY, GATKVCFConstants.ORIGINAL_AN_KEY);
                    copyAnnotationIfPresent(cohortVC, variantContextBuilder, GATKSVVCFConstants.SVTYPE, GATKSVVCFConstants.SVTYPE);
                    copyAnnotationIfPresent(cohortVC, variantContextBuilder, GATKSVVCFConstants.SVLEN, GATKSVVCFConstants.SVLEN);
                    //the joint segmentation goes through a lot of trouble to correct genotypes for overlapping events, so prefer those
                    if (cohortVC.hasGenotype(sampleName)) {
                        final Genotype genotype = cohortVC.getGenotype(sampleName);
                        final List<Allele> cohortVCAlleles = genotype.getAlleles();
                        genotypeBuilder.alleles(cohortVCAlleles);
                        //if ref is still N, update it
                        if (!uniquifiedAlleles.contains(cohortVC.getReference())) {
                            uniquifiedAlleles.remove(REF_ALLELE);
                            uniquifiedAlleles.add(cohortVC.getReference());
                        }
                        if (!genotype.isNoCall() && !uniquifiedAlleles.containsAll(cohortVCAlleles)) {
                            uniquifiedAlleles.addAll(cohortVCAlleles);
                        }
                        variantContextBuilder.alleles(uniquifiedAlleles);
                    }
                } else {
                    logger.warn("No matching cohort VC at " + segment.getContig() + ":" + segment.getStart());
                }
            } catch (final IOException e) {
                throw new GATKException("Error querying file " + clusteredCohortVcf + " over interval " +
                        new SimpleInterval(segment.getContig(), segment.getStart(), segment.getEnd()), e);
            }
        }
        final Genotype genotype = genotypeBuilder.make();
        variantContextBuilder.genotypes(genotype);
        variantContextBuilder.log10PError(segment.getQualitySomeCalled()/-10.0);  //use QS as site-level QUAL
        //apply filters if we're running against the cohort VCF
        if (clusteredCohortVcf != null) {
            variantContextBuilder.filters(getInfoFilters(segment, alleleFrequency));
        }
        return variantContextBuilder.make();
    }

    /**
     * Check if an annotation exists in the source VC and add it to the builder if found
     * @param cohortVC variant context source of annotations
     * @param variantContextBuilder is modified to add attribute (if found)
     * @param annotationKeyToQuery look for this key in the source
     * @param annotationKeyToWrite write this key back to VC builder
     */
    private void copyAnnotationIfPresent(final VariantContext cohortVC, final VariantContextBuilder variantContextBuilder, final String annotationKeyToQuery, final String annotationKeyToWrite) {
        if (cohortVC.hasAttribute(annotationKeyToQuery)) {
            final Object cohortValue = cohortVC.getAttribute(annotationKeyToQuery);
            if (cohortValue != null) {
                variantContextBuilder.attribute(annotationKeyToWrite, cohortValue);
            }
        }
    }

    /**
     * Determine quality and site frequency filters for a copy number event
     * @param segment copy number event with its quality metrics
     * @param alleleFrequency threshold above which common events are filtered out
     * @return
     */
    private Set<String> getInfoFilters(final IntegerCopyNumberSegment segment, final List<Double> alleleFrequency) {
        final Set<String> returnFilters = new LinkedHashSet<>();
        final int qsThreshold;
        //determine event type and applicable threshold
        if (segment.getCallIntegerCopyNumberState().getCopyNumber() == 0) {
            qsThreshold = homDelQSThreshold;
        } else if (segment.getCallIntegerCopyNumberState().getCopyNumber() > segment.getBaselineIntegerCopyNumberState().getCopyNumber()) {
            qsThreshold = dupeQSThreshold;
        } else if (segment.getCallIntegerCopyNumberState().getCopyNumber() < segment.getBaselineIntegerCopyNumberState().getCopyNumber()) {
            qsThreshold = hetDelQSThreshold;
        } else {
            qsThreshold = 0;
        }
        if (segment.getQualitySomeCalled() < qsThreshold) {
            returnFilters.add(GATKSVVCFConstants.LOW_QS_SCORE_FILTER_KEY);
        }
        if (alleleFrequency.stream().allMatch(af -> af >= siteFreqThreshold)) {
            returnFilters.add(GATKSVVCFConstants.FREQUENCY_FILTER_KEY);
        }

        return returnFilters;
    }
}
