package org.broadinstitute.hellbender.tools.sv;

import com.google.common.collect.Sets;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.cluster.CanonicalSVCollapser;
import org.broadinstitute.hellbender.tools.sv.cluster.PloidyTable;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static org.broadinstitute.hellbender.tools.sv.SVCallRecord.UNDEFINED_LENGTH;

public final class SVCallRecordUtils {

    /**
     * Create a builder for a variant from an {@link SVCallRecord} for VCF interoperability
     * @param record variant to convert
     * @return variant context builder
     */
    public static VariantContextBuilder getVariantBuilder(final SVCallRecord record) {
        Utils.nonNull(record);
        final int end;
        if (record.getType().equals(StructuralVariantType.INS) || record.getType().equals(StructuralVariantType.BND)) {
            end = record.getPositionA();
        } else {
            end = record.getPositionB();
        }
        final int end2;
        if (record.getType().equals(StructuralVariantType.INS)) {
            end2 = record.getPositionA();
        } else {
            end2 = record.getPositionB();
        }

        final List<Allele> altAlleles = record.getAltAlleles();
        final Allele refAllele = record.getRefAllele();
        final int numAlleles = 1 + altAlleles.size();
        final List<Allele> alleles = new ArrayList<>(numAlleles);
        if (refAllele == null) {
            // If no reference allele, use N
            alleles.add(Allele.REF_N);
        } else {
            alleles.add(refAllele);
        }
        alleles.addAll(altAlleles);

        final VariantContextBuilder builder = new VariantContextBuilder(record.getId(), record.getContigA(), record.getPositionA(),
                end, alleles);
        final StructuralVariantType svtype = record.getType();
        builder.id(record.getId());
        builder.attributes(record.getAttributes());
        builder.attribute(VCFConstants.END_KEY, end);
        builder.attribute(GATKSVVCFConstants.SVTYPE, svtype);
        builder.attribute(GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE, record.getAlgorithms());
        if (svtype.equals(StructuralVariantType.BND)) {
            builder.attribute(GATKSVVCFConstants.CONTIG2_ATTRIBUTE, record.getContigB());
            builder.attribute(GATKSVVCFConstants.END2_ATTRIBUTE, end2);
        }
        if (svtype.equals(StructuralVariantType.INS)) {
            builder.attribute(GATKSVVCFConstants.SVLEN, record.getLength());
        }
        if (svtype.equals(StructuralVariantType.BND) || svtype.equals(StructuralVariantType.INV)) {
            builder.attribute(GATKSVVCFConstants.STRANDS_ATTRIBUTE, getStrandString(record));
        }
        final GenotypesContext genotypes = GenotypesContext.create(record.getGenotypes().size());
        for (final Genotype g : record.getGenotypes()) {
            genotypes.add(sanitizeEmptyGenotype(g));
        }
        // htsjdk vcf encoder does not allow genotypes to have empty alleles
        builder.genotypes(record.getGenotypes().stream().map(SVCallRecordUtils::sanitizeEmptyGenotype).collect(Collectors.toList()));
        return builder;
    }

    /**
     * Adds NO_CALL allele if empty
     */
    private static Genotype sanitizeEmptyGenotype(final Genotype g) {
        if (g.getAlleles().isEmpty()) {
            return new GenotypeBuilder(g).alleles(Collections.singletonList(Allele.NO_CALL)).make();
        } else {
            return g;
        }
    }

    /**
     * Populates genotypes for samples not present in the given record. Samples with existing genotypes are not
     * touched. Samples without genotypes are assigned one according to the provided default reference allele
     * and ploidy, specified by the {@link GATKSVVCFConstants#EXPECTED_COPY_NUMBER_FORMAT} value. If the record
     * represents a CNV, the {@link GATKSVVCFConstants#COPY_NUMBER_FORMAT} is also set.
     * @param record record containing the genotypes
     * @param samples samples which the resulting genotypes must contain (existing samples are ignored)
     * @param refAlleleDefault default allele to use for samples without genotypes
     * @param ploidyTable ploidy table, which must contain at least all samples with missing genotypes
     * @return genotypes augmented with missing samples
     */
    public static GenotypesContext populateGenotypesForMissingSamplesWithAlleles(final SVCallRecord record,
                                                                                 final Set<String> samples,
                                                                                 final boolean refAlleleDefault,
                                                                                 final PloidyTable ploidyTable) {
        Utils.nonNull(record);
        Utils.nonNull(samples);
        final GenotypesContext genotypes = record.getGenotypes();
        final Set<String> missingSamples = Sets.difference(samples, genotypes.getSampleNames());
        if (missingSamples.isEmpty()) {
            return genotypes;
        }
        final ArrayList<Genotype> newGenotypes = new ArrayList<>(genotypes.size() + missingSamples.size());
        newGenotypes.addAll(genotypes);
        final String contig = record.getContigA();
        final List<Allele> altAlleles = record.getAltAlleles();
        final Allele refAllele = record.getRefAllele();
        final boolean isCNV = record.isSimpleCNV();
        for (final String sample : missingSamples) {
            final GenotypeBuilder genotypeBuilder = new GenotypeBuilder(sample);
            final int ploidy = ploidyTable.get(sample, contig);
            genotypeBuilder.attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, ploidy);
            if (isCNV) {
                genotypeBuilder.attribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT, ploidy);
                genotypeBuilder.alleles(CanonicalSVCollapser.getCNVGenotypeAllelesFromCopyNumber(altAlleles, refAllele, ploidy, ploidy));
            } else {
                genotypeBuilder.alleles(Collections.nCopies(ploidy, refAlleleDefault ? refAllele : Allele.NO_CALL));
            }
            newGenotypes.add(genotypeBuilder.make());
        }
        return GenotypesContext.create(newGenotypes);
    }

    /**
     * Creates shallow copy of the given record with genotypes replaced.
     * @param record base record
     * @param genotypes replacement genotypes
     * @return new record
     */
    public static <T extends SVCallRecord> SVCallRecord copyCallWithNewGenotypes(final T record, final GenotypesContext genotypes) {
        return new SVCallRecord(record.getId(), record.getContigA(), record.getPositionA(), record.getStrandA(), record.getContigB(),
                record.getPositionB(), record.getStrandB(), record.getType(), record.getLength(), record.getAlgorithms(), record.getAlleles(),
                genotypes, record.getAttributes());
    }

    /**
     * Creates shallow copy of the given record with attributes replaced.
     * @param record base record
     * @param attr replacement attributes
     * @return new record
     */
    public static <T extends SVCallRecord> SVCallRecord copyCallWithNewAttributes(final T record, final Map<String, Object> attr) {
        return new SVCallRecord(record.getId(), record.getContigA(), record.getPositionA(), record.getStrandA(), record.getContigB(),
                record.getPositionB(), record.getStrandB(), record.getType(), record.getLength(), record.getAlgorithms(), record.getAlleles(),
                record.getGenotypes(), attr);
    }

    /**
     * Get string representation for the given record's strands.
     */
    private static String getStrandString(final SVCallRecord record) {
        return getStrandString(record.getStrandA()) + getStrandString(record.getStrandB());
    }

    /**
     * Get string representation for a single strand.
     */
    private static String getStrandString(final boolean forwardStrand) {
        return forwardStrand ? SVCallRecord.STRAND_PLUS : SVCallRecord.STRAND_MINUS;
    }

    /**
     * Gets a comparator for {@link SVCallRecord} objects. Note we do not implement SVCallRecord as {@link Comparable}
     * due to the need for a common sequence dictionary. This would not only incur a cost for storing the dictionary
     * reference, but also for checking that the dictionaries of two records match. Note this is the same pattern used for
     * {@link IntervalUtils#compareLocatables(Locatable, Locatable, SAMSequenceDictionary) compareLocatables}.
     * @param dictionary sequence dictionary pertaining to both records (not validated)
     * @param <T> record class
     * @return comparator
     */
    public static <T extends SVCallRecord> Comparator<T> getCallComparator(final SAMSequenceDictionary dictionary) {
        return (o1, o2) -> compareCalls(o1, o2, dictionary);
    }

    /**
     * Gets a comparator for objects implementing {@link SVLocatable}.
     * @param dictionary sequence dictionary pertaining to both records (not validated)
     * @param <T> record class
     * @return comparator
     */
    public static <T extends SVLocatable> Comparator<T> getSVLocatableComparator(final SAMSequenceDictionary dictionary) {
        return (o1, o2) -> compareSVLocatables(o1, o2, dictionary);
    }

    /**
     * Compares two objects based on start and end positions.
     */
    public static int compareSVLocatables(final SVLocatable first, final SVLocatable second, final SAMSequenceDictionary dictionary) {
        Utils.nonNull(first);
        Utils.nonNull(second);
        // First locus
        final Comparator<Locatable> locatableComparator = IntervalUtils.getDictionaryOrderComparator(dictionary);
        final int compareA = locatableComparator.compare(new SimpleInterval(first.getContigA(), first.getPositionA(), first.getPositionA()),
                new SimpleInterval(second.getContigA(), second.getPositionA(), second.getPositionA()));
        if (compareA != 0) return compareA;
        // Second locus
        final int compareB = locatableComparator.compare(new SimpleInterval(first.getContigB(), first.getPositionB(), first.getPositionB()),
                new SimpleInterval(second.getContigB(), second.getPositionB(), second.getPositionB()));
        return compareB;
    }

    /**
     * Compare two records based by performing comparisons in the following order: start position, end position, start
     * strand, end strand, length, SV type. The first non-equal comparison determines the result.
     */
    public static int compareCalls(final SVCallRecord first, final SVCallRecord second, final SAMSequenceDictionary dictionary) {
        final int compareLocatables = compareSVLocatables(first, second, dictionary);
        if (compareLocatables != 0) return compareLocatables;

        // Type
        final int compareType = first.getType().compareTo(second.getType());
        if (compareType != 0) return compareType;

        //Strands
        if (first.getStrandA() == null && second.getStrandA() != null) {
            return -1;
        } else if (first.getStrandA() != null && second.getStrandA() == null) {
            return 1;
        } else if (first.getStrandA() != null && second.getStrandA() != null) {
            final int compareStartStrand = Boolean.compare(first.getStrandA(), second.getStrandA());
            if (compareStartStrand != 0) return compareStartStrand;
        }

        if (first.getStrandB() == null && second.getStrandB() != null) {
            return -1;
        } else if (first.getStrandB() != null && second.getStrandB() == null) {
            return 1;
        } else if (first.getStrandB() != null && second.getStrandB() != null) {
            final int compareEndStrand = Boolean.compare(first.getStrandB(), second.getStrandB());
            if (compareEndStrand != 0) return compareEndStrand;
        }

        // Length
        if (first.getLength() == null && second.getLength() != null) {
            return -1;
        } else if (first.getLength() != null && second.getLength() == null) {
            return 1;
        } else if (first.getLength() != null && second.getLength() != null) {
            return Integer.compare(first.getLength(), second.getLength());
        } else {
            return 0;
        }
    }

    /**
     * Converts the given record into a pair of BND records with "++" and "--" strandedness.
     * @param record inversion record
     * @return stream of BND records pair, or the original record if not an INV
     */
    public static Stream<SVCallRecord> convertInversionsToBreakends(final SVCallRecord record, final SAMSequenceDictionary dictionary) {
        if (!record.getType().equals(StructuralVariantType.INV)) {
            return Stream.of(record);
        }
        Utils.validateArg(record.isIntrachromosomal(), "Inversion " + record.getId() + " is not intrachromosomal");
        final SVCallRecord positiveBreakend = new SVCallRecord(record.getId(), record.getContigA(),
                record.getPositionA(), true, record.getContigB(), record.getPositionB(), true, StructuralVariantType.BND, null,
                record.getAlgorithms(), record.getAlleles(), record.getGenotypes(), record.getAttributes(), dictionary);
        final SVCallRecord negativeBreakend = new SVCallRecord(record.getId(), record.getContigA(),
                record.getPositionA(), false, record.getContigB(), record.getPositionB(), false, StructuralVariantType.BND, null,
                record.getAlgorithms(), record.getAlleles(), record.getGenotypes(), record.getAttributes(), dictionary);
        return Stream.of(positiveBreakend, negativeBreakend);
    }

    /**
     * Creates a new {@link SVCallRecord} from the given {@link VariantContext}, keeping any variant fields.
     * @see SVCallRecordUtils#create(VariantContext, boolean)
     */
    public static SVCallRecord create(final VariantContext variant) {
        return create(variant, true);
    }

    /**
     * Creates a new {@link SVCallRecord} from the given {@link VariantContext}.
     * @param variant SV record
     * @param keepVariantAttributes retain variant attribute fields
     * @return converted record
     */
    public static SVCallRecord create(final VariantContext variant, boolean keepVariantAttributes) {
        Utils.nonNull(variant);
        final String id = variant.getID();
        final String contigA = variant.getContig();
        final int positionA = variant.getStart();

        final StructuralVariantType type = inferStructuralVariantType(variant);
        final List<String> algorithms = getAlgorithms(variant);

        final String strands;
        if (type.equals(StructuralVariantType.DEL) || type.equals(StructuralVariantType.INS) ||
                type.equals(StructuralVariantType.CNV) || type.equals(StructuralVariantType.DUP)) {
            // SVCallRecord class can resolve these
            strands = null;
        } else {
            strands = getStrands(variant, type);
        }
        final Boolean strand1 = strands == null ? null : strands.startsWith(SVCallRecord.STRAND_PLUS);
        final Boolean strand2 = strands == null ? null : strands.endsWith(SVCallRecord.STRAND_PLUS);

        final Integer length;
        if (type.equals(StructuralVariantType.BND) || type.equals(StructuralVariantType.DEL) ||
                type.equals(StructuralVariantType.DUP) || type.equals(StructuralVariantType.CNV) ||
                type.equals(StructuralVariantType.INV)) {
            // SVCallRecord class can resolve these
            length = null;
        } else {
            length = getLength(variant, type);
        }

        final Map<String, Object> attributes = keepVariantAttributes ? variant.getAttributes() : Collections.emptyMap();

        final String contigB;
        final int positionB;
        if (type.equals(StructuralVariantType.BND)) {
            // If END2 and CHROM2 are both defined, use those.
            // If neither is defined, use start contig and position
            final boolean hasContig2 = variant.hasAttribute(GATKSVVCFConstants.CONTIG2_ATTRIBUTE);
            final boolean hasEnd2 = variant.hasAttribute(GATKSVVCFConstants.END2_ATTRIBUTE);
            if (!(hasContig2 && hasEnd2)) {
                throw new UserException.BadInput("Attributes " + GATKSVVCFConstants.END2_ATTRIBUTE +
                        " and " + GATKSVVCFConstants.CONTIG2_ATTRIBUTE + " are required for BND records (variant " +
                        variant.getID() + ").");
            }
            contigB = variant.getAttributeAsString(GATKSVVCFConstants.CONTIG2_ATTRIBUTE, null);
            positionB = variant.getAttributeAsInt(GATKSVVCFConstants.END2_ATTRIBUTE, 0);
        } else {
            contigB = contigA;
            positionB = variant.getEnd();
        }
        final Map<String, Object> sanitizedAttributes = sanitizeAttributes(attributes);
        return new SVCallRecord(id, contigA, positionA, strand1, contigB, positionB, strand2, type, length, algorithms,
                variant.getAlleles(), variant.getGenotypes(), sanitizedAttributes);
    }

    private static Map<String, Object> sanitizeAttributes(final Map<String, Object> attributes) {
        final Map<String, Object> newAttributes = new HashMap<>(attributes);
        for (final String key : SVCallRecord.INVALID_ATTRIBUTES) {
            newAttributes.remove(key);
        }
        return newAttributes;
    }

    private static Integer getLength(final VariantContext variant, final StructuralVariantType type) {
        Utils.nonNull(variant);
        Utils.validateArg(variant.hasAttribute(GATKSVVCFConstants.SVLEN), "Expected " + GATKSVVCFConstants.SVLEN + " field" + " for variant " + variant.getID());
        final int length = variant.getAttributeAsInt(GATKSVVCFConstants.SVLEN, UNDEFINED_LENGTH);
        if (length == UNDEFINED_LENGTH) {
            return null;
        }
        Utils.validate(length >= 0, "Length must be non-negative or " + UNDEFINED_LENGTH + " for variant " + variant.getID());
        return length;
    }

    public static List<String> getAlgorithms(final VariantContext variant) {
        Utils.nonNull(variant);
        Utils.validateArg(variant.hasAttribute(GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE), "Expected " + GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE + " field for variant " + variant.getID());
        return variant.getAttributeAsStringList(GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE, null);
    }

    private static String getStrands(final VariantContext variant, final StructuralVariantType type) {
        Utils.nonNull(variant);
        Utils.nonNull(type);
        final String strandsAttr = variant.getAttributeAsString(GATKSVVCFConstants.STRANDS_ATTRIBUTE, null);
        if (variant.getStructuralVariantType() == StructuralVariantType.INV && strandsAttr == null) {
            return null;
        }
        Utils.validateArg(strandsAttr != null, "Strands field not found for variant " + variant.getID() + " of type " + type);
        if (strandsAttr.length() != 2) {
            throw new IllegalArgumentException("Strands field is not 2 characters long for variant " + variant.getID());
        }
        final String strand1Char = strandsAttr.substring(0, 1);
        if (!strand1Char.equals(SVCallRecord.STRAND_PLUS) && !strand1Char.equals(SVCallRecord.STRAND_MINUS)) {
            throw new IllegalArgumentException("Valid start strand not found for variant " + variant.getID());
        }
        final String strand2Char = strandsAttr.substring(1, 2);
        if (!strand2Char.equals(SVCallRecord.STRAND_PLUS) && !strand2Char.equals(SVCallRecord.STRAND_MINUS)) {
            throw new IllegalArgumentException("Valid end strand not found for variant " + variant.getID());
        }
        return strandsAttr;
    }

    /**
     * Attempts to determine SV type from of a variant. If it is not explicitly available (i.e. through
     * {@link VariantContext#getStructuralVariantType()}) then the type is inferred from the alt alleles. The only
     * supported multi-allelic type is CNV when the alleles are DEL/DUP. Otherwise, a single symbolic alt allele is expected.
     */
    public static StructuralVariantType inferStructuralVariantType(final VariantContext variant) {
        final StructuralVariantType type = variant.getStructuralVariantType();
        if (type != null) {
            return type;
        }
        final List<Allele> alleles = variant.getAlternateAlleles();
        Utils.validate(!alleles.isEmpty(), "Missing alt allele for variant " + variant.getID());
        if (alleles.size() == 2 && alleles.contains(GATKSVVCFConstants.DEL_ALLELE) && alleles.contains(GATKSVVCFConstants.DUP_ALLELE)) {
            return StructuralVariantType.CNV;
        }
        Utils.validate(alleles.size() == 1, "Non-CNV multiallelic variants not supported (variant " + variant.getID() + ")");
        final Allele allele = alleles.get(0);
        Utils.validate(allele.isSymbolic(), "Expected symbolic alt allele");
        // TODO use htsjdk (see https://github.com/samtools/htsjdk/issues/18)
        return StructuralVariantType.valueOf(allele.getDisplayString().replace("<", "").replace(">", ""));
    }

    public static boolean containsAltAllele(final Genotype g) {
        return g.getAlleles().stream().anyMatch(SVCallRecordUtils::isAltAllele);
    }

    public static boolean isAltGenotype(final Genotype g) {
        return g.getAlleles().stream().anyMatch(SVCallRecordUtils::isAltAllele);
    }

    public static boolean isAltAllele(final Allele allele) {
        return allele != null && !allele.isNoCall() && !allele.isReference();
    }

    public static boolean isNonRefAllele(final Allele allele) {
        return allele != null && !allele.isReference();
    }

    // TODO this is sort of hacky but the Allele compareTo() method doesn't give stable ordering
    public static List<Allele> sortAlleles(final Collection<Allele> alleles) {
        return alleles.stream().sorted(Comparator.nullsFirst(Comparator.comparing(Allele::getDisplayString))).collect(Collectors.toList());
    }
}
