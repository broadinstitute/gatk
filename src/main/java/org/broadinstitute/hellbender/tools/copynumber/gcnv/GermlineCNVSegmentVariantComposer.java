package org.broadinstitute.hellbender.tools.copynumber.gcnv;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.hellbender.tools.copynumber.GermlineCNVCaller;
import org.broadinstitute.hellbender.tools.copynumber.PostprocessGermlineCNVCalls;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.IntegerCopyNumberSegmentCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.IntegerCopyNumberSegment;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Set;

/**
 * Helper class for {@link PostprocessGermlineCNVCalls} for single-sample postprocessing of segmented
 * {@link GermlineCNVCaller} calls.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class GermlineCNVSegmentVariantComposer {

    private static final String VARIANT_PREFIX = "CNV";

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

    private final String sampleName;
    private final VariantContextWriter outputWriter;
    private final IntegerCopyNumberState refAutosomalCopyNumberState;
    private final Set<String> allosomalContigSet;

    public static final Allele REF_ALLELE = Allele.create("N", true);
    public static final Allele DEL_ALLELE = Allele.create("<DEL>", false);
    public static final Allele DUP_ALLELE = Allele.create("<DUP>", false);
    public static final List<Allele> ALL_ALLELES = Arrays.asList(REF_ALLELE, DEL_ALLELE, DUP_ALLELE);

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
        this.outputWriter = Utils.nonNull(outputWriter);
        this.sampleName = Utils.nonEmpty(sampleName);
        this.refAutosomalCopyNumberState = Utils.nonNull(refAutosomalCopyNumberState);
        this.allosomalContigSet = Utils.nonNull(allosomalContigSet);
    }

    /**
     * Compose the header of the variant context
     *
     * @param vcfDefaultToolHeaderLines default header lines of the VCF generation tool
     */
    public void composeVariantContextHeader(final Set<VCFHeaderLine> vcfDefaultToolHeaderLines) {
        final VCFHeader result = new VCFHeader(Collections.emptySet(), Collections.singletonList(sampleName));

        /* add VCF version */
        result.addMetaDataLine(new VCFHeaderLine(VCFHeaderVersion.VCF4_2.getFormatString(),
                VCFHeaderVersion.VCF4_2.getVersionString()));

        /* add default tool header lines */
        vcfDefaultToolHeaderLines.forEach(line -> result.addMetaDataLine(line));

        /* header lines related to genotype formatting */
        result.addMetaDataLine(new VCFFormatHeaderLine(VCFConstants.GENOTYPE_KEY, 1,
                VCFHeaderLineType.Integer, "Segment genotype"));
        result.addMetaDataLine(new VCFFormatHeaderLine(CN, 1,
                VCFHeaderLineType.Integer, "Segment most-likely copy-number call"));
        result.addMetaDataLine(new VCFFormatHeaderLine(NP, 1,
                VCFHeaderLineType.Integer, "Number of points (i.e. targets or bins) in the segment"));
        result.addMetaDataLine(new VCFFormatHeaderLine(QS, 1,
                VCFHeaderLineType.Integer, "Phred-scaled quality of one or more of the points in the segment " +
                "agreeing with the most-likely copy-number call (normalized by the number of points)"));
        result.addMetaDataLine(new VCFFormatHeaderLine(QA, 1,
                VCFHeaderLineType.Integer, "Phred-scaled quality of all of the points in the segment " +
                "agreeing with the most-likely copy-number call"));
        result.addMetaDataLine(new VCFFormatHeaderLine(QSS, 1,
                VCFHeaderLineType.Integer, "Phred-scaled quality of determining the segment start position"));
        result.addMetaDataLine(new VCFFormatHeaderLine(QSE, 1,
                VCFHeaderLineType.Integer, "Phred-scaled quality of determining the segment end position"));

        /* INFO header lines */
        result.addMetaDataLine(new VCFInfoHeaderLine(VCFConstants.END_KEY, 1,
                VCFHeaderLineType.Integer, "End coordinate of the variant"));
        outputWriter.writeHeader(result);
    }

    /**
     * Write all variants.
     */
    public void writeAll(final IntegerCopyNumberSegmentCollection collection) {
        Utils.nonNull(collection);
        Utils.validateArg(collection.getMetadata().getSampleName().equals(sampleName),
                String.format("The sample name in the given collection differs from the expected " +
                        "sample name (collection: %s, expected: %s).",
                        collection.getMetadata().getSampleName(), sampleName));
        for (final IntegerCopyNumberSegment segment : collection.getRecords()) {
            final VariantContext variant = composeSegmentVariantContext(segment);
            outputWriter.add(variant);
        }
    }

    /**
     * Compose a variant context from a given {@link IntegerCopyNumberSegment}
     *
     * @param segment an instance of {@link IntegerCopyNumberSegment}
     * @return composed variant context
     */
    @VisibleForTesting
    VariantContext composeSegmentVariantContext(final IntegerCopyNumberSegment segment) {
        final String contig = segment.getContig();
        final int start = segment.getStart();
        final int end = segment.getEnd();
        int copyNumberCall = segment.getCallIntegerCopyNumberState().getCopyNumber();

        final VariantContextBuilder variantContextBuilder = new VariantContextBuilder();
        variantContextBuilder.alleles(ALL_ALLELES);
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

        variantContextBuilder.attribute(VCFConstants.END_KEY, end);
        variantContextBuilder.genotypes(genotype);
        return variantContextBuilder.make();
    }
}
