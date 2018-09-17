package org.broadinstitute.hellbender.tools.spark.sv.utils;

import htsjdk.samtools.*;
import htsjdk.variant.variantcontext.StructuralVariantType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.engine.spark.datasources.ReferenceMultiSparkSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.ReadMetadata;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;

import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * Variant context with additional method to mine the structural variant specific information from structural
 * variant records.
 */
public final class SVContext extends VariantContext {

    private static final long serialVersionUID = 1L;

    /**
     * Indicates that the 'end' has not been calculated yet.
     */
    private static final int MISSING_END = -2;

    /**
     * Indicates that the length has not been calculated yet.
     */
    private static final int MISSING_LENGTH = -2;

    /**
     * Indicates that the variant does not have a length or this could not be determined, e.g. for BND records.
     */
    public static final int NO_LENGTH = -1;

    private int end = MISSING_END;
    private int length = MISSING_LENGTH;

    /**
     * Returns an instance given a plain {@link VariantContext} instance.
     * <p>
     *     This method will fail with a {@link IllegalArgumentException} if the input variant context is not
     *     structural. A structural variant context must be biallelic and properly annotated with
     *     with its type using the INFO field {@value VCFConstants#SVTYPE}.
     *
     * </p>
     * @param vc the input variant context.
     * @throws IllegalArgumentException if {@code vc} is {@code null} or does not seem to be
     * a structural variant context based on its annotations.
     * @return never {@code null}, potentially the same as the input if it happens to be an instance of this class.
     */
    public static SVContext of(final VariantContext vc) {
        if (vc instanceof SVContext) {
            return (SVContext) vc;
        }
        assertIsStructuralVariantContext(vc);
        return new SVContext(vc);
    }

    private static void assertIsStructuralVariantContext(final VariantContext vc) {
        Utils.nonNull(vc, "the input variant context must not be null");
        Utils.nonNull(vc.getStructuralVariantType(), "the input variant-context is not structural; the SVTYPE annotation is missing");
        if (vc.getNAlleles() != 2) {
            throw new IllegalArgumentException("structural variant context must be biallelic");
        } else if (vc.getAlleles().get(0).isNonReference()) {
            throw new IllegalArgumentException("the allele with index 0 must be reference");
        } else if (vc.getAlleles().get(1).isReference()) {
            throw new IllegalArgumentException("the allele with index 1 must be non-reference");
        }
    }

    private SVContext(final VariantContext other) {
        super(other);
    }

    /**
     * The end position of this variant context.
     * This end position is determined in this order:
     * <ul>
     *     <li>the value of the annotation {@value VCFConstants#END_KEY} if present, otherwise</li>
     *     <li>the value that results from applying the length annotated under {@value  GATKSVVCFConstants#SVLEN} based on the structural type if such annotation is present, otherwise</li>
     *     <li>the variant context start position + the length of the reference allele - 1 as per {@link VariantContext#getEnd()}.</li>
     * </ul>
     * @return whatever is returned by {@link #getStart()} or greater.
     */
    @Override
    public int getEnd() {
        if (end == MISSING_END) {
            end = getAttributeAsInt(VCFConstants.END_KEY, MISSING_END);
            if (end == MISSING_END) { // need to test a second time as "END=." would result in still a missing-end.
                if (getStructuralVariantType() == StructuralVariantType.DEL) {
                    end = getStart() + getStructuralVariantLength();
                } else {
                    end = super.getEnd();
                }
            }
        }
        return end;
    }

    /**
     * Returns the ids of the assembled contig that support this context's structural variant.
     * <p>
     *     The list returned is an immutable list.
     * </p>
     * @throws IllegalStateException if the {@link GATKSVVCFConstants#CONTIG_NAMES} annotation contains
     * undefined contig names (.)
     *
     * @return never {@code null}, an empty list if no structural variant is specified.
     */
    public List<String> getSupportingContigIds() {
        if (!hasAttribute(GATKSVVCFConstants.CONTIG_NAMES)) {
            return Collections.emptyList();
        } else {
            final List<String> result = getAttributeAsStringList(GATKSVVCFConstants.CONTIG_NAMES, null);
            if (result.contains(null)) {
                throw new IllegalStateException("the contig names annotation contains undefined values");
            }
            return result;
        }
    }

    /**
     * Composes the haplotype that corresponds to an allele based on the reference sequence.
     * <p>
     *     The resulting haplotype will be annotated with the appropriate Cigar and genome location.
     * </p>
     * @param index the index of the target allele.
     * @param paddingSize extra bases from the reference sequence to be added on either side.
     * @param reference the reference to use as source.
     * @return never {@code null}.
     */
    public Haplotype composeHaplotypeBasedOnReference(final int index, final int paddingSize, final ReferenceMultiSparkSource reference)  {
        Utils.nonNull(reference, "the input reference cannot be null");
        ParamUtils.isPositiveOrZero(paddingSize, "the input padding must be 0 or greater");
        ParamUtils.inRange(index, 0, 1, "the input allele index must be 0 or 1");

        final SAMSequenceDictionary dictionary = reference.getReferenceSequenceDictionary(null);
        Utils.nonNull(dictionary, "the input reference does not have a dictionary");

        final SAMSequenceRecord contigRecord = dictionary.getSequence(getContig());
        Utils.nonNull(contigRecord, "the input reference does not have a contig named: " + getContig());

        final int contigLength = contigRecord.getSequenceLength();
        if (contigLength < getEnd()) {
            throw new IllegalArgumentException(String.format("this variant goes beyond the end of "
                    + "the containing contig based on the input reference: "
                    + "contig %s length is %d but variant end is %d", getContig(), contigLength, getEnd()));
        }

        final SimpleInterval referenceInterval = composePaddedInterval(getContig(), contigLength, getStart(), getEnd(), paddingSize);

        final ReferenceBases bases;
        try {
            bases = reference.getReferenceBases(referenceInterval);
        } catch (final IOException ex) {
            throw new GATKException("could not read reference file");
        }
        if (index == 0) {
            final Haplotype result = new Haplotype(bases.getBases(), true);
            result.setCigar(new Cigar(Collections.singletonList(new CigarElement(referenceInterval.size(), CigarOperator.M))));
            result.setGenomeLocation(referenceInterval);
            return result;
        } else { //index == 1
            switch (getStructuralVariantType()) {
                case INS:
                    return composeInsertionHaplotype(bases);
                case DEL:
                    return composeDeletionHaplotype(bases);
                default: // not jet supported. Please add more types as needed.
                    throw new UnsupportedOperationException("not supported yet");
            }
        }
    }

    private Haplotype composeDeletionHaplotype(final ReferenceBases referenceBases) {
        final int deletionSize = getStructuralVariantLength();
        final byte[] resultBases = new byte[referenceBases.getInterval().size() - deletionSize];
        final int leftPaddingSize = getStart() - referenceBases.getInterval().getStart() + 1;
        final int rightPaddingSize = referenceBases.getInterval().getEnd() - getStart() - deletionSize;
        final byte[] referenceBaseBytes = referenceBases.getBases();
        System.arraycopy(referenceBaseBytes, 0, resultBases, 0, leftPaddingSize);
        System.arraycopy(referenceBaseBytes, leftPaddingSize + deletionSize, resultBases, leftPaddingSize, rightPaddingSize);
        final Cigar cigar = new Cigar(Arrays.asList(new CigarElement(leftPaddingSize, CigarOperator.M),
                new CigarElement(deletionSize, CigarOperator.D),
                new CigarElement(rightPaddingSize, CigarOperator.M)));
        final Haplotype result = new Haplotype(resultBases, false);
        result.setCigar(cigar);
        result.setGenomeLocation(referenceBases.getInterval());
        return result;
    }

    private Haplotype composeInsertionHaplotype(final ReferenceBases referenceBases) {
        final byte[] insertedSequence = getInsertedSequence();
        final byte[] referenceBaseBytes = referenceBases.getBases();
        final byte[] resultBases = new byte[referenceBases.getInterval().size() + insertedSequence.length];
        final int leftPaddingSize = getStart() - referenceBases.getInterval().getStart() + 1;
        final int rightPaddingSize = referenceBases.getInterval().getEnd() - getStart();
        System.arraycopy(referenceBaseBytes, 0, resultBases, 0, leftPaddingSize);
        System.arraycopy(insertedSequence, 0, resultBases, leftPaddingSize, insertedSequence.length);
        System.arraycopy(referenceBaseBytes, leftPaddingSize    , resultBases, leftPaddingSize + insertedSequence.length, rightPaddingSize);
        final Cigar cigar = new Cigar(Arrays.asList(new CigarElement(leftPaddingSize, CigarOperator.M),
                new CigarElement(insertedSequence.length, CigarOperator.I),
                new CigarElement(rightPaddingSize, CigarOperator.M)));
        final Haplotype result = new Haplotype(resultBases, false);
        result.setCigar(cigar);
        result.setGenomeLocation(referenceBases.getInterval());
        return result;
    }

    /**
     * Returns the inserted sequence.
     * <p>
     *     Currently this method relies on the annotation {@value GATKSVVCFConstants#INSERTED_SEQUENCE}.
     * </p>
     * @return {@code null} if there is no inserted sequence.
     */
    public byte[] getInsertedSequence() {
        if (hasAttribute(GATKSVVCFConstants.INSERTED_SEQUENCE)) {
            final String asString = getAttributeAsString(GATKSVVCFConstants.INSERTED_SEQUENCE, null);
            return asString == null ? null : asString.getBytes();
        } else {
            return null;
        }
    }

    /**
     * Returns the absolute structural variant length as recorded in the SVLEN annotation.
     * <p>
     *     When SVLEN is absent this method returns {@link #NO_LENGTH} ({@value #NO_LENGTH}).
     * <p>
     *     Otherwise this method will return 0 or positive value, so for example for a deletion of 10bp
     *     it will return {@code 10} despite that the SVLEN annotation should be {@code -10}.
     * </p>
     * @return {@link #NO_LENGTH} if there is no SVLEN annotation and the length could not be inferred, 0 or greater otherwise.
     */
    public int getStructuralVariantLength() {
        if (length == MISSING_LENGTH) {
            length = Math.abs(getAttributeAsInt(GATKSVVCFConstants.SVLEN, NO_LENGTH));
        }
        return length;
    }

    /**
     * Returns break point intervals for this structural variant.
     * <p>
     *     Typically the break point interval would be located between two point locations in the reference genome.
     *     In that case this method will return {@code padding} bases up- and down-stream from that
     *     inter-base position.
     * </p>
     * <p>
     *     In case the padding is set to zero, since the 0-length interval is not valid, this method would
     *     return an interval including the base just before the break-point.
     * </p>
     * <p>
     *      If padForHomology is set, the breakpoint interval will include the region specified as homologous
     *      sequence in the HOMOLOGY_LENGTH attribute of VariantContext, in addition to the
     *      normal padding specified in the first parameter.
     * </p>
     *
     * @param padding the padding around the exact location of the break point to be included in the padded interval.
     * @param dictionary reference meta-data.
     * @param padForHomology Add homologous sequence around the breakpoint to the interval
     * @return never {@code null}, potentially 0-length but typically at least one element.
     */
    public List<SimpleInterval> getBreakPointIntervals(final int padding,
                                                       final SAMSequenceDictionary dictionary,
                                                       final boolean padForHomology) {
        ParamUtils.isPositiveOrZero(padding, "the input padding must be 0 or greater");
        Utils.nonNull(dictionary, "the input dictionary cannot be null");
        final String contigName = getContig();
        final int contigLength = dictionary.getSequence(contigName).getSequenceLength();
        final StructuralVariantType type = getStructuralVariantType();
        final int start = getStart();
        if (type == StructuralVariantType.INS) {
            return Collections.singletonList(
                    composePaddedInterval(contigName, contigLength, start, start, padding));
        } else if (type == StructuralVariantType.DEL) {
            final int end = getEnd();
            final int homologyPadding = (padForHomology ? getAttributeAsInt(GATKSVVCFConstants.HOMOLOGY_LENGTH, 0) : 0);
            return Arrays.asList(
                    composePaddedInterval(contigName, contigLength, start + 1,
                            start + 1 + homologyPadding, padding),
                    composePaddedInterval(contigName, contigLength, end,
                            end + homologyPadding, padding));
        } else {
            // Please, add more types as needed!
            throw new UnsupportedOperationException("currently only supported for INS and DELs");
        }
    }

    private static SimpleInterval composePaddedInterval(final String contig, final int contigSize, final int start,
                                                        final int end, final int padding) {
        return new SimpleInterval(contig,
                Math.max(1, padding > 0 ? start - padding + 1 : start),
                Math.min(contigSize, end + padding));

    }

    public PairedStrandedIntervals getPairedStrandedIntervals(final ReadMetadata metadata, final SAMSequenceDictionary referenceSequenceDictionary, final int padding) {
        final StructuralVariantType type = getStructuralVariantType();
        if (type == StructuralVariantType.DEL) {
            final List<SimpleInterval> breakPointIntervals = getBreakPointIntervals(padding, referenceSequenceDictionary, true);
            final SimpleInterval leftBreakpointSimpleInterval = breakPointIntervals.get(0);
            final SVInterval leftBreakpointInterval = new SVInterval(
                    metadata.getContigID(leftBreakpointSimpleInterval.getContig()),
                    leftBreakpointSimpleInterval.getStart(),
                    leftBreakpointSimpleInterval.getEnd() + 1);
            final SimpleInterval rightBreakpointSimpleInterval = breakPointIntervals.get(1);
            final SVInterval rightBreakpointInterval = new SVInterval(
                    metadata.getContigID(rightBreakpointSimpleInterval.getContig()),
                    rightBreakpointSimpleInterval.getStart(),
                    rightBreakpointSimpleInterval.getEnd() + 1);

            return new PairedStrandedIntervals(new StrandedInterval(leftBreakpointInterval, true),
                    new StrandedInterval(rightBreakpointInterval, false));
        } else {
            throw new UnsupportedOperationException("currently only supported for DELs");
        }
    }
}