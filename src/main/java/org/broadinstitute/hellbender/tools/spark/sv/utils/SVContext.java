package org.broadinstitute.hellbender.tools.spark.sv.utils;

import htsjdk.samtools.*;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.StructuralVariantType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.hadoop.io.MD5Hash;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignedContig;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.ReadMetadata;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

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

    /**
     * Holds the unique id for this structural variant. A {@link null} indicates that the ID has not yet been
     * calculated.
     */
    public String uid = null;

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
     * @param paddingSize extra bases from the reference sequence to be added on either side.
     * @param reference the reference to use as source.
     * @return never {@code null}.
     */
    public List<SVHaplotype> composeHaplotypesBasedOnReference(final int paddingSize,
                                                        final ReferenceMultiSource reference,
                                                        final List<AlignedContig> contigs)  {
        Utils.nonNull(reference, "the input reference cannot be null");
        ParamUtils.isPositiveOrZero(paddingSize, "the input padding must be 0 or greater");

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

        final List<ReferenceBases> referenceBases = calculateReferenceIntervals(paddingSize, dictionary, contigs).stream()
                .map(pi -> { try { return reference.getReferenceBases(null, pi); } catch (final IOException ex) { throw new UncheckedIOException(ex); }})
                .collect(Collectors.toList());
        final List<SVHaplotype> result = new ArrayList<>(2);

        final SVHaplotype referenceHaplotype = composeReferenceHaplotype(referenceBases);

        result.add(referenceHaplotype);
        switch (getStructuralVariantType()) {
            case INS:
                    result.add(composeInsertionHaplotype(referenceBases));
                    break;
            case DEL:
                    result.add(composeDeletionHaplotype(referenceBases));
                    break;
            default: // not jet supported. Please add more types as needed.
                    throw new UnsupportedOperationException("not supported yet");
        }
        return result;
    }

    private SVHaplotype composeReferenceHaplotype(final List<ReferenceBases> referenceBases) {
        final byte[] bases = Utils.concat(referenceBases.stream().map(ReferenceBases::getBases).toArray(byte[][]::new));
        final List<CigarElement> cigarElements = new ArrayList<>();
        SimpleInterval previousInterval;
        cigarElements.add(new CigarElement((previousInterval = referenceBases.get(0).getInterval()).getLengthOnReference(), CigarOperator.M));
        for (int i = 1; i < referenceBases.size(); i++) {
            final ReferenceBases nextBases = referenceBases.get(i);
            final SimpleInterval nextInterval = nextBases.getInterval();
            cigarElements.add(new CigarElement(nextInterval.getStart() - previousInterval.getEnd() - 1, CigarOperator.N));
            cigarElements.add(new CigarElement(nextInterval.getLengthOnReference(), CigarOperator.M));
            previousInterval = nextInterval;
        }
        final Cigar cigar = new Cigar(cigarElements);
        final SimpleInterval startInterval = referenceBases.get(0).getInterval();
        return new ArraySVHaplotype(SVHaplotype.REF_HAPLOTYPE_NAME,
                Collections.singletonList(new AlignmentInterval(startInterval.getContig(), startInterval.getStart(), true, cigar,
                        SAMRecord.UNKNOWN_MAPPING_QUALITY, 0, AlignmentInterval.NO_AS)), bases, this.getUniqueID(), getStartPositionInterval(), false);
    }

    /**
     * Composes a list of the relevant interval in the reference when it comes to compose the haplotypes.
     *
     * <p>
     *     This list results from applying the requested padding to all break points extending to include
     *     any overlapping contig.
     * </p>
     *
     * <p>
     *     Thus, contig or contigs alignment intervals that fall far from the break points will be ignored,
     *     whereas those that are close enough to the break points wil extend the padding so that these are
     *     fully included.
     * </p>
     *
     * <p>
     *     The resulting list is sorted by the reference dictionary and contains strictly non-overlapping
     *     intervals; if they overlap the would have been merged into a single element.
     * </p>
     *
     * @param padding
     * @param dictionary
     * @param contigs
     * @return
     */
    private List<SimpleInterval> calculateReferenceIntervals(final int padding, final SAMSequenceDictionary dictionary,
                                                             final List<AlignedContig> contigs) {
        final List<SimpleInterval> breakPoints = getBreakPointIntervals(padding, dictionary, true);
        SVIntervalLocator locator = new SVIntervalLocator(dictionary);
        final SVIntervalTree<SVInterval> tree = new SVIntervalTree<>();
        Stream.concat(breakPoints.stream().map(si -> locator.toSVInterval(si)),
                contigs.stream().flatMap(ctg -> ctg.alignmentIntervals.stream())
                        .map(ai -> locator.toSVInterval(ai.referenceSpan, padding + 1)))
                .forEach(svIntervalPlus -> {
                    final Iterator<SVIntervalTree.Entry<SVInterval>> overlappers = tree.overlappers(svIntervalPlus);
                    int minStart = svIntervalPlus.getStart() + 1; // + 1 and -1 to remove the aboutter padding.
                    int maxEnd = svIntervalPlus.getEnd() - 1;
                    while (overlappers.hasNext()) {
                        final SVInterval overlapper = overlappers.next().getInterval();
                        if (overlapper.getStart() < minStart) { minStart = overlapper.getStart(); }
                        if (overlapper.getEnd() > maxEnd) { maxEnd = overlapper.getEnd(); }
                        overlappers.remove();
                    }
                    final SVInterval newInterval = new SVInterval(svIntervalPlus.getContig(), minStart, maxEnd);
                    tree.put(newInterval, newInterval);
                });
        return breakPoints.stream()
                .map(bp -> locator.toSVInterval(bp, 1))
                .map(svi -> tree.overlappers(svi).next().getInterval())
                .map(locator::toSimpleInterval)
                .map(ivl -> {
                    final int seqLength = dictionary.getSequence(ivl.getContig()).getSequenceLength();
                    return seqLength >= ivl.getEnd() ? ivl : new SimpleInterval(ivl.getContig(), ivl.getStart(), seqLength);
                })
                .distinct().collect(Collectors.toList());
    }

    private SVHaplotype composeDeletionHaplotype(final List<ReferenceBases> allReferenceBases) {
        final ReferenceBases leftReferenceBases = allReferenceBases.stream()
                .filter(rb -> rb.getInterval().contains(getStartPositionInterval()))
                .findFirst()
                .orElseThrow(() -> new IllegalStateException("input reference bases does not contain the variant position"));
        final ReferenceBases rightReferenceBases = allReferenceBases.stream()
                .filter(rb -> rb.getInterval().contains(getStopPositionInterval()))
                .findFirst()
                .orElseThrow(() -> new IllegalStateException("input reference bases does not contain the stop position"));

        final int deletionSize = getStructuralVariantLength();
        final byte[] leftBases = leftReferenceBases.getSubsetBases(leftReferenceBases.getInterval().getStart(), getStart());
        final byte[] rightBases = rightReferenceBases.getSubsetBases(getStart() + getStructuralVariantLength() + 1, rightReferenceBases.getInterval().getEnd());
        final byte[] resultBases = Utils.concat(leftBases, rightBases);
        final Cigar cigar = new Cigar(Arrays.asList(new CigarElement(leftBases.length, CigarOperator.M),
                new CigarElement(deletionSize, CigarOperator.D),
                new CigarElement(rightBases.length, CigarOperator.M)));
        return new ArraySVHaplotype(SVHaplotype.ALT_HAPLOTYPE_NAME,
                Collections.singletonList(new AlignmentInterval(leftReferenceBases.getInterval().getContig(),
                        leftReferenceBases.getInterval().getStart(), true, cigar, SAMRecord.UNKNOWN_MAPPING_QUALITY, 0,
                        AlignmentInterval.NO_AS)),
                resultBases, getUniqueID(), getStartPositionInterval(), false);
    }

    private SVHaplotype composeInsertionHaplotype(final List<ReferenceBases> allReferenceBases) {
        final ReferenceBases referenceBases = allReferenceBases.stream()
                .filter(rb -> rb.getInterval().contains(getStartPositionInterval()))
                .findFirst()
                .orElseThrow(() -> new IllegalStateException("input reference bases does not contain the variant position"));

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
        final SimpleInterval referenceBasesInterval = referenceBases.getInterval();
        return new ArraySVHaplotype(SVHaplotype.ALT_HAPLOTYPE_NAME,
                Collections.singletonList(new AlignmentInterval(referenceBasesInterval.getContig(), referenceBasesInterval.getStart(), true, cigar, SAMRecord.UNKNOWN_MAPPING_QUALITY, 0, AlignmentInterval.NO_AS)),
                resultBases, getUniqueID(), getStartPositionInterval(), false);
    }

    /**
     * Returns the inserted sequence.
     * <p>
     *     Currently this method relies on the annotation {@value GATKSVVCFConstants#INSERTED_SEQUENCE}.
     * </p>
     * @return {@code null} if there is no inserted sequence.
     */
    public byte[] getInsertedSequence() {
        if (hasAttribute("INSERTED_SEQUENCE")) {
            final String asString = getAttributeAsString("INSERTED_SEQUENCE", null);
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

    public static int getStructuralVariantLength(final VariantContext vc) {
        Utils.nonNull(vc);
        if (vc.hasAttribute(GATKSVVCFConstants.SVLEN)) {
            return Math.abs(vc.getAttributeAsInt(GATKSVVCFConstants.SVLEN, NO_LENGTH));
        } else {
            return NO_LENGTH;
        }
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

    /**
     * Returns a unique id for this structural variant based on its location, type and
     * attributes.
     * @return never {@code null}.
     */
    public String getUniqueID() {
        if (uid == null) {
            uid = composeUniqueID();
        }
        return uid;
    }

    private String composeUniqueID() {
        final StringBuilder builder = new StringBuilder(50);
        builder.append("sv/");
        if (hasID()) {
            builder.append(getID()).append('/');
        }

        builder.append(getStructuralVariantType().name().toLowerCase()).append('/');
        final String attributeString = getAttributes().entrySet().stream()
                .map(entry -> entry.getKey() + "=" + String.valueOf(entry.getValue()))
                .sorted()
                .collect(Collectors.joining(";"));
        builder.append("dg:").append(Integer.toHexString(MD5Hash.digest(attributeString).quarterDigest()));
        return builder.toString();
    }

    public SimpleInterval getStartPositionInterval() {
        return new SimpleInterval(getContig(), getStart(), getStart());
    }

    public Locatable getStopPositionInterval() {
        return new SimpleInterval(getContig(), getEnd(), getEnd());
    }
}