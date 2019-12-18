package org.broadinstitute.hellbender.tools.copynumber.gcnv;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.hellbender.tools.copynumber.GermlineCNVCaller;
import org.broadinstitute.hellbender.tools.copynumber.PostprocessGermlineCNVCalls;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.IntegerCopyNumberSegment;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;

/**
 * Helper class for {@link PostprocessGermlineCNVCalls} for single-sample postprocessing of segmented
 * {@link GermlineCNVCaller} calls.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class GermlineCNVSegmentVariantComposer extends GermlineCNVVariantComposer<IntegerCopyNumberSegment> {

    /* VCF FORMAT header keys */

    /**
     * Segment copy-number call
     */
    public static final String CN = "CN";

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

    private final IntegerCopyNumberState refAutosomalCopyNumberState;
    private final Set<String> allosomalContigSet;

    /**
     * Constructor.
     *
     * @param outputWriter variant context writer
     * @param sampleName sample name
     * @param refAutosomalCopyNumberState ref copy-number state on autosomal contigs
     * @param allosomalContigSet set of allosomal contigs (ref copy-number allele be chosen according to
     *                           given contig baseline copy-number states)
     */
    public GermlineCNVSegmentVariantComposer(final VariantContextWriter outputWriter,
                                             final String sampleName,
                                             final IntegerCopyNumberState refAutosomalCopyNumberState,
                                             final Set<String> allosomalContigSet) {
        super(outputWriter, sampleName);
        this.refAutosomalCopyNumberState = Utils.nonNull(refAutosomalCopyNumberState);
        this.allosomalContigSet = Utils.nonNull(allosomalContigSet);
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

        result.setSequenceDictionary(sequenceDictionary);

        /* header lines related to genotype formatting */
        result.addMetaDataLine(new VCFFormatHeaderLine(VCFConstants.GENOTYPE_KEY, 1,
                VCFHeaderLineType.Integer, "Segment genotype"));
        result.addMetaDataLine(new VCFFormatHeaderLine(CN, 1,
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
        int copyNumberCall = segment.getCallIntegerCopyNumberState().getCopyNumber();

        final VariantContextBuilder variantContextBuilder = new VariantContextBuilder();
        variantContextBuilder.chr(contig);
        variantContextBuilder.start(start);
        variantContextBuilder.stop(end);
        variantContextBuilder.id(String.format(VARIANT_PREFIX + "_%s_%d_%d", contig, start, end));

        final GenotypeBuilder genotypeBuilder = new GenotypeBuilder(sampleName);
        final IntegerCopyNumberState refCopyNumber = allosomalContigSet.contains(contig)
                ? segment.getBaselineIntegerCopyNumberState()
                : refAutosomalCopyNumberState;
        final Allele allele;
        if (copyNumberCall > refCopyNumber.getCopyNumber()) {
            allele = DUP_ALLELE;
        } else if (copyNumberCall < refCopyNumber.getCopyNumber()) {
            allele = DEL_ALLELE;
        } else {
            allele = REF_ALLELE;
        }
        genotypeBuilder.alleles(Collections.singletonList(allele));
        genotypeBuilder.attribute(CN, copyNumberCall);
        genotypeBuilder.attribute(NP, segment.getNumPoints());
        genotypeBuilder.attribute(QS, FastMath.round(segment.getQualitySomeCalled()));
        genotypeBuilder.attribute(QA, FastMath.round(segment.getQualityAllCalled()));
        genotypeBuilder.attribute(QSS, FastMath.round(segment.getQualityStart()));
        genotypeBuilder.attribute(QSE, FastMath.round(segment.getQualityEnd()));
        final Genotype genotype = genotypeBuilder.make();

        final List<Allele> vcAlleles = new ArrayList<>(Collections.singletonList(REF_ALLELE));
        if (!allele.equals(REF_ALLELE)) {
            vcAlleles.add(allele);
        }
        variantContextBuilder.alleles(vcAlleles);
        variantContextBuilder.attribute(VCFConstants.END_KEY, end);
        variantContextBuilder.genotypes(genotype);
        variantContextBuilder.log10PError(segment.getQualitySomeCalled()/-10.0);
        return variantContextBuilder.make();
    }
}
