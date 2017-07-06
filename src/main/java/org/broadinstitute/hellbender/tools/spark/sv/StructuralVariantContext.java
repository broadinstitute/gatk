package org.broadinstitute.hellbender.tools.spark.sv;

import htsjdk.samtools.*;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
import htsjdk.variant.variantcontext.VariantContext;
import org.bdgenomics.formats.avro.StructuralVariant;
import org.bdgenomics.formats.avro.StructuralVariantType;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.CigarUtils;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.seqdoop.hadoop_bam.VariantContextWritable;

import java.io.IOException;
import java.util.*;
import java.util.function.Supplier;

/**
 * Created by valentin on 4/26/17.
 */
public class StructuralVariantContext extends VariantContext {

    private static final long serialVersionUID = 1L;

    private transient List<Integer> assemblyIDs;

    protected StructuralVariantContext(final VariantContext other) {
        super(other);
    }

    public boolean isStructural() {
        return getNAlleles() > 1 && getAlternateAlleles().stream()
                .anyMatch(StructuralVariantAllele::isStructural);
    }

    /**
     * Returns the assembly ids for this context's structural variant.
     * <p>
     *     The list returned is an immutable list.
     * </p>
     * @throws IllegalStateException if the {@link GATKSVVCFHeaderLines#ASSEMBLY_IDS} annotation contains invalid
     *
     * @return never {@code null}, an empty list if no structural variant is specified.
     */
    public List<Integer> assemblyIDs() {
        if (assemblyIDs == null) {
            if (!hasAttribute(GATKSVVCFHeaderLines.ASSEMBLY_IDS)) {
                assemblyIDs = Collections.emptyList();
            } else {
                try {
                    final List<Integer> result = getAttributeAsIntList(GATKSVVCFHeaderLines.ASSEMBLY_IDS, -1);
                    if (result.stream().anyMatch(i -> i < 0)) {
                        throw new IllegalStateException("this variant context contains unspecified or negative assembly ids");
                    }
                    assemblyIDs = Collections.unmodifiableList(result);
                } catch (final NumberFormatException ex) {
                    throw new IllegalArgumentException(String.format("the variant context '%s' annotation contains non-integer entries: %s", GATKSVVCFHeaderLines.ASSEMBLY_IDS,
                            String.valueOf(getAttribute(GATKSVVCFHeaderLines.ASSEMBLY_IDS))));
                }

            }
        }
        return assemblyIDs;
    }

    public Haplotype composeHaplotype(final int index, final int paddingSize, final ReferenceMultiSource reference)  {
        if (index < 0 || index > 1)
            throw new IllegalArgumentException("wrong index must be 0 (ref) or 1 (alt)");
        Utils.nonNull(reference);

        final SimpleInterval referenceInterval = new SimpleInterval(getContig(), getStart() - paddingSize - 1, getStart() + paddingSize + Math.abs(Math.min(0,getStructuralVariantLength())));
        final ReferenceBases bases;
        try {
            bases = reference.getReferenceBases(null, referenceInterval);
        } catch (final IOException ex) {
            throw new UserException.CouldNotReadInputFile("could not read reference file");
        }
        if (index == 0) {
            final Haplotype result = new Haplotype(bases.getBases(), true);
            result.setCigar(new Cigar(Collections.singletonList(new CigarElement(referenceInterval.size(), CigarOperator.M))));
            result.setGenomeLocation(referenceInterval);
            return result;
        } else { //index == 1
            switch (getStructuralVariantAllele()) {
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
        final int deletionSize = - getStructuralVariantLength();
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
        System.arraycopy(referenceBaseBytes, rightPaddingSize, resultBases, leftPaddingSize + insertedSequence.length, rightPaddingSize);
        final Cigar cigar = new Cigar(Arrays.asList(new CigarElement(leftPaddingSize, CigarOperator.M),
                                                    new CigarElement(insertedSequence.length, CigarOperator.I),
                                                    new CigarElement(rightPaddingSize, CigarOperator.M)));
        final Haplotype result = new Haplotype(resultBases, false);
        result.setCigar(cigar);
        result.setGenomeLocation(referenceBases.getInterval());
        return result;
    }

    public StructuralVariantAllele getStructuralVariantAllele() {
        if (!isStructural()) {
            throw new IllegalStateException("this is not in fact a structural variant context");
        } else {
            return StructuralVariantAllele.valueOf(alleles.get(1));
        }
    }

    public byte[] getInsertedSequence() {
        if (hasAttribute(GATKSVVCFHeaderLines.INSERTED_SEQUENCE)) {
            final String asString = getAttributeAsString(GATKSVVCFHeaderLines.INSERTED_SEQUENCE, null);
            if (asString == null) {
                return null;
            } else {
                return asString.getBytes();
            }
        } else {
            return null;
        }
    }

    public int getStructuralVariantLength() {

        final Supplier<List<String>> s = ArrayList<String>::new;
        if (hasAttribute(GATKSVVCFHeaderLines.SVLEN)) {
            return getAttributeAsInt(GATKSVVCFHeaderLines.SVLEN, 0);
        } else {
            throw new IllegalStateException("missing insertion length");
        }
    }
}
