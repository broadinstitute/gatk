package org.broadinstitute.hellbender.tools.spark.sv;

import htsjdk.samtools.*;
import htsjdk.variant.variantcontext.StructuralVariantType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFHeaderLines;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.function.Supplier;

/**
 * Variant context with additional method to mind the structural variant specific information from structural
 * variant records.
 */
public class StructuralVariantContext extends VariantContext {

    private static final long serialVersionUID = 1L;

    private int end = -1;
    private int length = -1;

    /**
     * Returns a instance given a plain {@link VariantContext} instance.
     * <p>
     *     This method will fail with a {@link IllegalArgumentException} if the input variant context is not
     *     structural. Currently we only check that it has been properly annotated with the {@value VCFConstants#SVTYPE}
     *     annotation.
     * </p>
     * @param vc the input variant context.
     * @throws IllegalArgumentException if {@code vc} is {@code null} or does not seem to be
     * a structural variant context based on its annotations.
     * @return never {@code null}, potentially the same as the input if it happens to be an instance of this class.
     */
    public static StructuralVariantContext create(final VariantContext vc) {
        if (vc instanceof StructuralVariantContext) {
            return (StructuralVariantContext) vc;
        } else {
            checkIsStructuralVariantContext(vc);
            return new StructuralVariantContext(vc);
        }
    }

    private static void checkIsStructuralVariantContext(final VariantContext vc) {
        Utils.nonNull(vc, "the input variant context must not be null");
        Utils.nonNull(vc.getAttributeAsString(VCFConstants.SVTYPE, null), "the SVTYPE annotation is missing");
    }

    private StructuralVariantContext(final VariantContext other) {
        super(other);
    }

    /**
     * The end position of this variant context.
     * This end position is determined with this priority:
     * <ul>
     *     <li>the value of the annotation {@value VCFConstants#END_KEY} if present, otherwise</li>
     *     <li>the value that results from applying the length annotated under {@value  GATKSVVCFConstants#SVLEN} based on the structural type if such annoation is present, otherwise</li>
     *     <li>the variant context start position + the length of the reference allele - 1.</li>
     * </ul>
     * @return whatever is returned by {@link #getStart()} or greater.
     */
    @Override
    public int getEnd() {
        if (end <= 0) {
            if (hasAttribute(VCFConstants.END_KEY)) {
                end = getAttributeAsInt(VCFConstants.END_KEY, -1);
            } else if (getStructuralVariantType() == StructuralVariantType.DEL) {
                end = getStart() + getStructuralVariantLength();
            } else {
                end = super.getEnd();
            }
        }
        return end;
    }

    /**
     * Returns the assembly ids for this context's structural variant.
     * <p>
     *     The list returned is an immutable list.
     * </p>
     * @throws IllegalStateException if the {@link GATKSVVCFConstants#CONTIG_NAMES} annotation contains
     * undefined contig names (.)
     *
     * @return never {@code null}, an empty list if no structural variant is specified.
     */
    public List<String> getContigNames() {
        if (!hasAttribute(GATKSVVCFConstants.CONTIG_NAMES)) {
            return Collections.emptyList();
        } else {
            final List<String> contigNames = getAttributeAsStringList(GATKSVVCFConstants.CONTIG_NAMES, null);
            if (contigNames.contains(null)) {
                throw new IllegalStateException("the contig names annotation contains undefined values");
            }
            return contigNames;
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
    public Haplotype composeHaplotypeBasedOnReference(final int index, final int paddingSize, final ReferenceMultiSource reference)  {
        Utils.nonNull(reference, "the input reference cannot be null");
        ParamUtils.isPositiveOrZero(paddingSize, "the input padding must be 0 or greater");
        ParamUtils.inRange(index, 0, 1, "the input allele index must be 0 or 1");

        final SAMSequenceDictionary dictionary = reference.getReferenceSequenceDictionary(null);
        if (dictionary == null) {
            throw new IllegalArgumentException("the input reference does not have a dictionary");
        }
        final SAMSequenceRecord contigRecord = dictionary.getSequence(getContig());
        if (contigRecord == null) {
            throw new IllegalArgumentException("the input reference does not have a contig named: " + getContig());
        }
        final int contigLength = contigRecord.getSequenceLength();
        if (contigLength < getEnd()) {
            throw new IllegalArgumentException(String.format("this variant goes beyond the end of "
                    + "the containing contig based on the input reference: "
                    + "contig %s length is %d but variant end is %d", getContig(), contigLength, getEnd()));
        }
        final SimpleInterval referenceInterval = new SimpleInterval(getContig(),
                Math.max(1, getStart() - paddingSize + 1),
                Math.min(contigLength, getEnd() + paddingSize));

        final ReferenceBases bases;
        try {
            bases = reference.getReferenceBases(null, referenceInterval);
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
                default: // not jet supported.
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
            if (asString == null) {
                return null;
            } else {
                return asString.getBytes();
            }
        } else {
            return null;
        }
    }

    /**
     * Returns the absolute structural variant length as recorded in the SVLEN annotation.
     * <p>
     *     This method always return a 0 or positive value, so for example for a deletion of 10bp
     *     it will return 10 despite that the SVLEN annotation should be -10.
     * </p>
     * @return -1 if there is no SVLEN annotation, 0 or greater otherwise.
     */
    public int getStructuralVariantLength() {
        if (length < 0) {
            if (hasAttribute(GATKSVVCFConstants.SVLEN)) {
                length = Math.abs(getAttributeAsInt(GATKSVVCFConstants.SVLEN, 0));
            }
        }
        return length;
    }

    /**
     * Returns break point intervals for this structural variant.
     * <p>
     *     Typically the break point would be located between to bases in the reference genome.
     *     In that case this method will return {@code padding} bases up- and down-stream from that
     *     inter-base position.
     * </p>
     * <p>
     *     In case the padding is set two zero, since the 0-length interval is not valid, this method would
     *     return an interval included the base just before the break-point.
     * </p>
     *
     * @param padding the padding around the exact location of the break point to be included in the padding interval.
     * @param dictionary reference meta-data.
     * @return never {@code null}, potentially 0-length but typically at least one element.
     */
    public List<SimpleInterval> getBreakPointIntervals(final int padding, final SAMSequenceDictionary dictionary) {
        Utils.nonNull(dictionary, "the input dictionary cannot be null");
        final String contigName = getContig();
        final int contigLength = dictionary.getSequence(contigName).getSequenceLength();
        final StructuralVariantType type = getStructuralVariantType();
        final int start = getStart();
        if (type == StructuralVariantType.INS) {
            return Collections.singletonList(new SimpleInterval(contigName,
                    Math.max(1, padding > 0 ? (start - padding + 1) : start),
                    Math.min(contigLength, start + padding)));
        } else if (type == StructuralVariantType.DEL) { // must be <DEL>
            final int length = getStructuralVariantLength();
            return Arrays.asList(
                    new SimpleInterval(contigName,
                            Math.max(1, padding > 0 ? (start - padding + 1) : start) ,
                            Math.min(contigLength, start + padding)), // left-end.
                    new SimpleInterval(contigName,
                            Math.max(1, length + (padding > 0 ? (start - padding + 1) : start)),
                            Math.min(contigLength, start + length + padding)));
        } else {
            throw new UnsupportedOperationException("currently only supported for INS and DELs");
        }
    }
}