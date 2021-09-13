package org.broadinstitute.hellbender.tools.sv.cluster;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.variant.variantcontext.*;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVCallRecordUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.reference.ReferenceUtils;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Class for collapsing a collection of similar {@link SVCallRecord} objects, such as clusters produced by
 * {@link CanonicalSVLinkage}, into a single representative call.
 */
public class CanonicalSVCollapser implements SVCollapser<SVCallRecord> {

    /**
     * Determines whether a given SV type represents a CNV.
     */
    public static boolean isCnvType(final StructuralVariantType type) {
        return type == StructuralVariantType.DEL || type == StructuralVariantType.DUP || type == StructuralVariantType.CNV;
    }

    /**
     * Define strategies for collapsing variant intervals.
     */
    public enum BreakpointSummaryStrategy {
        /**
         * Use the (first) middle value to summarize cluster starts and ends, such that the start and end were seen in the data
         */
        MEDIAN_START_MEDIAN_END,

        /**
         * A conservative strategy to summarize a cluster by its smallest extent
         */
        MIN_START_MAX_END,

        /**
         * A permissive strategy to summarize a cluster by it largest extent
         */
        MAX_START_MIN_END,

        /**
         * Summarize a cluster using the mean value for each end, even if that value was not represented in any sample
         */
        MEAN_START_MEAN_END

    }

    /**
     * Define strategies for collapsing alt alleles with different subtypes.
     */
    public enum AltAlleleSummaryStrategy {
        /**
         * Use the most specific subtype that doesn't conflict with any of the other alleles.
         * For example, (&lt;INS&gt, &lt;INS:MEI:SVA&gt, &lt;INS:MEI:LINE&gt) results in &lt;INS:MEI&gt.
         */
        MOST_SPECIFIC_SUBTYPE,

        /**
         * Use subtypes in common among all alleles.
         * For example, (&lt;INS&gt, &lt;INS:MEI:SVA&gt, &lt;INS:MEI&gt) results in &lt;INS&gt.
         */
        COMMON_SUBTYPE

    }

    /**
     * Strategies for collapsing the insertion LENGTH attribute. Only applied when there are multiple defined values.
     */
    public enum InsertionLengthSummaryStrategy {
        MEDIAN,
        MEAN,
        MIN,
        MAX,
        UNDEFINED
    }

    private final AltAlleleSummaryStrategy altAlleleSummaryStrategy;
    private final BreakpointSummaryStrategy breakpointSummaryStrategy;
    private final InsertionLengthSummaryStrategy insertionLengthSummaryStrategy;
    private final ReferenceSequenceFile reference;

    private static final AlleleCollectionCollapserComparator ALLELE_COMPARATOR = new AlleleCollectionCollapserComparator();

    public CanonicalSVCollapser(final ReferenceSequenceFile reference,
                                final AltAlleleSummaryStrategy altAlleleSummaryStrategy,
                                final BreakpointSummaryStrategy breakpointSummaryStrategy,
                                final InsertionLengthSummaryStrategy insertionLengthSummaryStrategy) {
        this.reference = Utils.nonNull(reference);
        this.altAlleleSummaryStrategy = altAlleleSummaryStrategy;
        this.breakpointSummaryStrategy = breakpointSummaryStrategy;
        this.insertionLengthSummaryStrategy = insertionLengthSummaryStrategy;
    }

    @Override
    public SVCallRecord collapse(final Collection<SVCallRecord> items) {
        final String id = collapseIds(items);
        final List<String> algorithms = collapseAlgorithms(items);
        final StructuralVariantType type = collapseTypes(items);
        final Map<String, Object> attributes = collapseVariantAttributes(items);

        // Prefer using variants generated with PESR callers, which tend to generate more precise breakpoints
        final Collection<SVCallRecord> mostPreciseCalls = getRecordsWithMostPreciseBreakpoints(items);
        final SVCallRecord exampleCall = mostPreciseCalls.iterator().next();

        final Pair<Integer, Integer> coordinates = collapseInterval(mostPreciseCalls);
        final int start = coordinates.getKey();
        final int end = coordinates.getValue();
        final int length = collapseLength(mostPreciseCalls, start, end, type);

        final Allele refAllele = collapseRefAlleles(exampleCall.getContigA(), start);
        final List<Allele> altAlleles = collapseAltAlleles(items, type);
        final int numAlleles = (refAllele == null ? 0 : 1) + altAlleles.size();
        final List<Allele> alleles = new ArrayList<>(numAlleles);
        if (refAllele != null) {
            alleles.add(refAllele);
        }
        alleles.addAll(altAlleles);
        final List<Genotype> genotypes = collapseAllGenotypes(items, refAllele);

        return new SVCallRecord(id, exampleCall.getContigA(), start, exampleCall.getStrandA(),
                exampleCall.getContigB(), end, exampleCall.getStrandB(), type, length, algorithms, alleles, genotypes, attributes);
    }

    /**
     * Finds representative expected copy number for a sample. This method simply checks whether this is defined in all
     * genotypes and verifies that there is only one unique value and returns it.
     * @param sampleGenotypes genotypes belonging to a single sample
     * @return representative copy number for the reference state
     */
    protected int collapseExpectedCopyNumber(final Collection<Genotype> sampleGenotypes) {
        Utils.nonNull(sampleGenotypes);
        Utils.nonEmpty(sampleGenotypes);
        final List<Integer> expectedCopyNumberValues = sampleGenotypes.stream()
                .map(SVCallRecord::getExpectedCopyNumber)
                .distinct()
                .collect(Collectors.toList());
        Utils.validate(expectedCopyNumberValues.size() == 1,
                "Expected 1 unique expected copy number but found " + expectedCopyNumberValues.size());
        return expectedCopyNumberValues.get(0);
    }

    /**
     * Returns the ref allele at the given locus.
     */
    protected Allele collapseRefAlleles(final String contig, final int pos) {
        final byte[] bases = ReferenceUtils.getRefBaseAtPosition(reference, contig, pos);
        Utils.validate(bases != null && bases.length == 1, "Invalid reference locus " + contig + ":" + pos);
        return Allele.create(bases[0], true);
    }

    /**
     * Parses allele for the base type string, e.g. "INS" for "&lt;INS:MEI&gt;"
     */
    @VisibleForTesting
    protected static final String[] getSymbolicAlleleSymbols(final Allele allele) {
        return allele.getDisplayString()
                .replace("<", "")
                .replace(">", "")
                .split(":");
    }

    /**
     * Collapses alternate alleles into a list of representative alleles. Note this supports sub-typed alleles such as
     * &lt;INS:MEI&gt;. If multiple alt alleles are found, the variant must either be a CNV or sub-typed alleles with the
     * same base symbol (e.g. &lt;INS:MEI&gt; and &lt;INS:MEI:SVA&gt; would result in &lt;INS&gt;).
     * @param items records whose alt alleles should be collapsed
     * @param type type of the collapsed records
     * @return collapsed alt alleles
     */
    protected List<Allele> collapseAltAlleles(final Collection<SVCallRecord> items, final StructuralVariantType type) {
        Utils.nonNull(items);
        final List<Allele> altAlleles = items.stream().map(SVCallRecord::getAltAlleles)
                .flatMap(List::stream)
                .distinct()
                .collect(Collectors.toList());
        if (altAlleles.isEmpty()) {
            return Collections.emptyList();
        } else if (altAlleles.size() == 1) {
            return Collections.singletonList(altAlleles.get(0));
        } else {
            // Multiple non-ref alleles need collapsing
            // TODO does not search for subtypes e.g. <DUP:TANDEM>
            if (altAlleles.size() == 2 && altAlleles.contains(Allele.SV_SIMPLE_DEL) && altAlleles.contains(Allele.SV_SIMPLE_DUP)) {
                // CNVs
                if (type == StructuralVariantType.CNV) {
                    return altAlleles;
                } else {
                    throw new IllegalArgumentException("Encountered multi-allelic with DEL/DUP for non-CNV");
                }
            }

            final String[] collapsedAlleleTokens;
            if (altAlleleSummaryStrategy == AltAlleleSummaryStrategy.COMMON_SUBTYPE) {
                collapsedAlleleTokens = collapseAltAllelesCommon(altAlleles);
            } else if (altAlleleSummaryStrategy == AltAlleleSummaryStrategy.MOST_SPECIFIC_SUBTYPE) {
                collapsedAlleleTokens = collapseAltAllelesMostSpecific(altAlleles);
            } else {
                throw new UnsupportedOperationException("Unimplemented alt allele summary strategy: " + altAlleleSummaryStrategy.name());
            }
            Utils.validate(collapsedAlleleTokens.length > 0, "Encountered multiple symbolic allele base symbols for non-CNV");
            return Collections.singletonList(Allele.create("<" + String.join(":", collapsedAlleleTokens) + ">", false));
        }
    }

    private String[] collapseAltAllelesCommon(final List<Allele> alleles) {
        final List<String[]> alleleTokens = alleles.stream()
                .map(CanonicalSVCollapser::getSymbolicAlleleSymbols)
                .collect(Collectors.toList());
        final String[] firstTokens = alleleTokens.get(0);
        int alleleSize = 0;
        for (int i = 0; i < firstTokens.length; i++) {
            final String subtype = firstTokens[i];
            for (int j = 1; j < alleleTokens.size(); j++) {
                final String[] tokens = alleleTokens.get(j);
                if (i < tokens.length && subtype.equals(tokens[i])) {
                    alleleSize = i + 1;
                } else {
                    return Arrays.copyOf(firstTokens, alleleSize);
                }
            }
        }
        return Arrays.copyOf(firstTokens, alleleSize);
    }

    private String[] collapseAltAllelesMostSpecific(final List<Allele> alleles) {
        final List<String[]> alleleTokens = alleles.stream()
                .map(CanonicalSVCollapser::getSymbolicAlleleSymbols)
                .collect(Collectors.toList());
        final int maxSize = alleleTokens.stream().mapToInt(arr -> arr.length).max().getAsInt();
        int alleleIndex = 0;
        int alleleSize = 0;
        outerloop:
        for (int i = 0; i < maxSize; i++) {
            String subtype = null;
            for (int j = 0; j < alleleTokens.size(); j++) {
                final String[] tokens = alleleTokens.get(j);
                if (i < tokens.length) {
                    if (subtype == null) {
                        subtype = tokens[i];
                        alleleIndex = j;
                        alleleSize = i + 1;
                    } else if (!subtype.equals(tokens[i])) {
                        alleleIndex = j;
                        alleleSize = i;
                        break outerloop;
                    }
                }
            }
        }
        return Arrays.copyOf(alleleTokens.get(alleleIndex), alleleSize);
    }

    private List<Genotype> collapseAllGenotypes(final Collection<SVCallRecord> items,
                                                final Allele refAllele) {
        return items.stream()
                .map(SVCallRecord::getGenotypes)
                .flatMap(GenotypesContext::stream)
                .collect(Collectors.groupingBy(Genotype::getSampleName))
                .values()
                .stream()
                .map(g -> collapseSampleGenotypes(g, refAllele))
                .collect(Collectors.toList());
    }

    /**
     * Collapses a set of genotypes from a single sample into a list of alleles for a representative genotype. Ploidy
     * is determined with the {@link #collapseExpectedCopyNumber(Collection)} method. The result is the most common non-ref genotype
     * alleles list, potentially augmented with additional ref alleles to match the ploidy.
     * @param genotypes genotypes from a single variant and single sample
     * @param refAllele ref allele for the genotype
     * @param sampleAltAlleles unique alt alleles
     * @return list of alleles for a representative genotype (may be empty)
     */
    @VisibleForTesting
    protected List<Allele> collapseSampleGenotypeAlleles(final Collection<Genotype> genotypes,
                                                         final int expectedCopyNumber,
                                                         final Allele refAllele,
                                                         final List<Allele> sampleAltAlleles) {
        Utils.nonNull(genotypes);
        Utils.nonEmpty(genotypes);

        // Ploidy 0
        if (expectedCopyNumber == 0) {
            return Collections.emptyList();
        }

        // No alt alleles, so return all ref
        //if (sampleAltAlleles.isEmpty()) {
        //    return Collections.nCopies(expectedCopyNumber, refAllele);
        //}

        // Determine frequency of each genotype
        final Map<List<Allele>, Integer> genotypeCounts = genotypes.stream().map(Genotype::getAlleles)
                .collect(Collectors.groupingBy(l -> l, Collectors.collectingAndThen(Collectors.toList(), List::size)));
        // Sort keys somehow, for stability
        final List<List<Allele>> sortedAlleleLists = genotypeCounts.keySet().stream().sorted(ALLELE_COMPARATOR).collect(Collectors.toList());
        // Find most common non-ref genotype
        // Ties go to the first result in the list, as determined by the comparator
        // TODO : implement an allele comparator that can take into account GQ if available
        List<Allele> bestGenotypeAltAlleles = null;
        int bestGenotypeFrequency = 0;
        for (final List<Allele> alleles : sortedAlleleLists) {
            final List<Allele> altAlleles = alleles.stream().filter(SVCallRecordUtils::isAltAllele).collect(Collectors.toList());
            final int genotypeFrequency = genotypeCounts.get(alleles);
            if (!altAlleles.isEmpty() && (bestGenotypeAltAlleles == null || genotypeFrequency > bestGenotypeFrequency)) {
                bestGenotypeAltAlleles = altAlleles;
                bestGenotypeFrequency = genotypeFrequency;
            }
        }

        // Create alleles list with the selected alt alleles, and fill in remaining with ref
        final List<Allele> alleles = new ArrayList<>(expectedCopyNumber);
        final int numCollapsedRefAlleles = expectedCopyNumber - bestGenotypeAltAlleles.size();
        Utils.validate(numCollapsedRefAlleles >= 0, "Expected copy number is less than number of alt alleles");
        for (int i = 0; i < numCollapsedRefAlleles; i++) {
            alleles.add(refAllele);
        }
        alleles.addAll(bestGenotypeAltAlleles);
        return alleles;
    }

    protected Genotype collapseSampleGenotypes(final Collection<Genotype> genotypes,
                                               final Allele refAllele) {
        final GenotypeBuilder builder = new GenotypeBuilder(genotypes.iterator().next());
        final List<Allele> sampleAltAlleles = genotypes.stream().map(Genotype::getAlleles)
                .flatMap(List::stream)
                .filter(SVCallRecordUtils::isAltAllele)
                .distinct()
                .collect(Collectors.toList());
        final int expectedCopyNumber = collapseExpectedCopyNumber(genotypes);
        final List<Allele> collapsedAlleles = collapseSampleGenotypeAlleles(genotypes, expectedCopyNumber, refAllele, sampleAltAlleles);
        builder.alleles(collapsedAlleles);
        builder.noAttributes();
        builder.attributes(collapseGenotypeAttributes(genotypes, collapsedAlleles, expectedCopyNumber));
        return builder.make();
    }

    protected Map<String, Object> collapseGenotypeAttributes(final Collection<Genotype> genotypes,
                                                             final List<Allele> collapsedAlleles,
                                                             final int expectedCopyNumber) {
        Utils.nonNull(genotypes);
        Utils.nonEmpty(genotypes);
        final List<Allele> collapsedAltAlleles = collapsedAlleles.stream().filter(SVCallRecordUtils::isAltAllele).collect(Collectors.toList());
        final Map<String, Object> collapsedAttributes = new HashMap<>();
        final Map<String, Set<Object>> genotypeFields = genotypes.stream().map(Genotype::getExtendedAttributes)
                .map(Map::entrySet)
                .flatMap(Set::stream)
                .collect(Collectors.groupingBy(Map.Entry::getKey, Collectors.mapping(Map.Entry::getValue, Collectors.toSet())));
        for (final Map.Entry<String, Set<Object>> entry : genotypeFields.entrySet()) {
            if (entry.getKey().equals(GATKSVVCFConstants.COPY_NUMBER_FORMAT)) {
                collapsedAttributes.put(entry.getKey(), collapseSampleCopyNumber(expectedCopyNumber, collapsedAltAlleles));
            } else if (entry.getKey().equals(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT)) {
                continue;
            } else {
                collapsedAttributes.put(entry.getKey(), collapseSampleGenotypeAttribute(entry.getKey(), entry.getValue()));
            }
        }
        collapsedAttributes.put(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, expectedCopyNumber);
        return collapsedAttributes;
    }

    protected Object collapseSampleGenotypeAttribute(final String key, final Set<Object> values) {
        if (values.size() == 1) {
            return values.iterator().next();
        } else {
            return null;
        }
    }

    /**
     * Special handling for CNV copy number attribute
     */
    protected Object collapseSampleCopyNumber(final int expectedCopyNumber,
                                              final List<Allele> collapsedAltAlleles) {
        int copyNumber = expectedCopyNumber;
        for (final Allele allele : collapsedAltAlleles) {
            if (allele.equals(Allele.SV_SIMPLE_DEL)) {
                copyNumber--;
            } else if (allele.equals(Allele.SV_SIMPLE_DUP)) {
                copyNumber++;
            }
        }
        return new Integer(copyNumber);
    }

    protected Map<String, Object> collapseVariantAttributes(final Collection<SVCallRecord> items) {
        Utils.nonNull(items);
        Utils.nonEmpty(items);
        final Map<String, Object> collapsedAttributes = new HashMap<>();
        final Map<String, Set<Object>> attributes = items.stream().map(SVCallRecord::getAttributes)
                .map(Map::entrySet)
                .flatMap(Set::stream)
                .collect(Collectors.groupingBy(Map.Entry::getKey, Collectors.mapping(Map.Entry::getValue, Collectors.toSet())));
        for (final Map.Entry<String, Set<Object>> entry : attributes.entrySet()) {
            collapsedAttributes.put(entry.getKey(), collapseSingleVariantAttribute(entry.getValue()));
        }
        collapsedAttributes.put(GATKSVVCFConstants.CLUSTER_MEMBER_IDS_KEY, items.stream().map(SVCallRecord::getId).sorted().collect(Collectors.toList()));
        return collapsedAttributes;
    }

    protected Object collapseSingleVariantAttribute(final Set<Object> values) {
        if (values.size() == 1) {
            return values.iterator().next();
        } else {
            return null;
        }
    }

    protected final int collapseLength(final Collection<SVCallRecord> items, final int newStart, final int newEnd,
                                       final StructuralVariantType newType) {
        Utils.nonNull(items);
        Utils.nonEmpty(items);
        final SVCallRecord exampleCall = items.iterator().next();
        if (!exampleCall.isIntrachromosomal()) {
            return -1;
        } else if (newType.equals(StructuralVariantType.INS)) {
            return collapseInsertionLength(items);
        } else {
            return newEnd - newStart + 1;
        }
    }

    protected Integer collapseInsertionLength(final Collection<SVCallRecord> items) {
        // No need for strategy when the answer is obvious
        if (items.size() == 1) {
            return items.iterator().next().getLength();
        }

        // We only want to use strategies when there are multiple defined lengths
        final int[] definedLengths = items.stream()
                .map(SVCallRecord::getLength)
                .filter(len -> len != null)
                .mapToInt(Integer::intValue)
                .toArray();

        // All lengths were undefined
        if (definedLengths.length == 0) {
            return null;
        }

        // Only one length is defined - no ambiguity here
        if (definedLengths.length == 1) {
            return definedLengths[0];
        }

        // Employ strategies when there are two or more defined lengths
        if (insertionLengthSummaryStrategy == InsertionLengthSummaryStrategy.MEDIAN) {
            return MathUtils.median(definedLengths, Percentile.EstimationType.R_3);
        } else if (insertionLengthSummaryStrategy == InsertionLengthSummaryStrategy.MEAN) {
            return (int) Math.ceil(MathUtils.sum(definedLengths) / (double) definedLengths.length);
        } else if (insertionLengthSummaryStrategy == InsertionLengthSummaryStrategy.MIN) {
            return MathUtils.arrayMin(definedLengths);
        } else if (insertionLengthSummaryStrategy == InsertionLengthSummaryStrategy.MAX) {
            return MathUtils.arrayMax(definedLengths);
        } else if (insertionLengthSummaryStrategy == InsertionLengthSummaryStrategy.UNDEFINED) {
            return null;
        } else {
            throw new UnsupportedOperationException("Unimplemented insertion summary strategy: " + insertionLengthSummaryStrategy.name());
        }
    }

    protected String collapseIds(final Collection<SVCallRecord> records) {
        Utils.nonNull(records);
        Utils.nonEmpty(records);
        return records.stream().map(SVCallRecord::getId).sorted().collect(Collectors.toList()).get(0);
    }

    /**
     * Gets records likely to have the most accurate breakpoint position. These usually are supported by PE/SR/AS support,
     * whereas depth-only calls tend to be approximate.
     */
    protected Collection<SVCallRecord> getRecordsWithMostPreciseBreakpoints(final Collection<SVCallRecord> items) {
        if (items.stream().allMatch(call -> call.isDepthOnly())) {
            return items;
        } else {
            return items.stream().filter(call -> !call.isDepthOnly()).collect(Collectors.toList());
        }
    }

    /**
     * @param items
     * @return (key, value) entry of (start, end)
     */
    protected Pair<Integer,Integer> collapseInterval(final Collection<SVCallRecord> items) {
        Utils.nonNull(items);
        Utils.nonEmpty(items);
        final SVCallRecord exampleCall = items.iterator().next();
        if (items.size() > 1) {
            final List<String> contigsA = items.stream().map(SVCallRecord::getContigA).distinct().collect(Collectors.toList());
            Utils.validate(contigsA.size() == 1, "Cannot collapse intervals with multiple position A contigs");
            final List<String> contigsB = items.stream().map(SVCallRecord::getContigB).distinct().collect(Collectors.toList());
            Utils.validate(contigsB.size() == 1, "Cannot collapse intervals with multiple position B contigs");
        }
        final int[] startPositions = items.stream().mapToInt(SVCallRecord::getPositionA).sorted().toArray();
        final int[] endPositions = items.stream().mapToInt(SVCallRecord::getPositionB).sorted().toArray();
        //use the mid value of the sorted list so the start and end represent real breakpoint observations
        final int medianStart = MathUtils.median(startPositions, Percentile.EstimationType.R_3);
        final int medianEnd = MathUtils.median(endPositions, Percentile.EstimationType.R_3);
        final int newStart;
        final int newEnd;
        if (exampleCall.getType().equals(StructuralVariantType.INS)) {
            // Insertions should be a single locus; also fixes case where end-supporting split reads are to the
            // left of start-supporting split reads
            final int mean = (medianStart + medianEnd) / 2;
            newStart = mean;
            newEnd = mean;
        } else {
            switch (breakpointSummaryStrategy) {
                case MEDIAN_START_MEDIAN_END:
                    newStart = medianStart;
                    newEnd = medianEnd;
                    break;
                case MIN_START_MAX_END:
                    newStart = startPositions[0];
                    newEnd = endPositions[endPositions.length - 1];
                    break;
                case MAX_START_MIN_END:
                    newStart = startPositions[startPositions.length - 1];
                    newEnd = endPositions[0];
                    break;
                case MEAN_START_MEAN_END:
                    newStart = (int)Math.round(MathUtils.sum(startPositions) / (double) startPositions.length);
                    newEnd = (int)Math.round(MathUtils.sum(endPositions) / (double) startPositions.length);
                    break;
                default:
                    throw new UnsupportedOperationException("Unknown breakpoint summary strategy: " + breakpointSummaryStrategy.name());
            }
        }
        return Pair.of(newStart, newEnd);
    }

    protected StructuralVariantType collapseTypes(final Collection<SVCallRecord> records) {
        final Set<StructuralVariantType> types = records.stream().map(SVCallRecord::getType).collect(Collectors.toSet());
        if (types.size() == 1) {
            return types.iterator().next();
        }
        if (types.stream().allMatch(CanonicalSVCollapser::isCnvType)) {
            return StructuralVariantType.CNV;
        }
        final List<String> typeStrings = types.stream().map(StructuralVariantType::name).collect(Collectors.toList());
        throw new IllegalArgumentException("Incompatible SV types found in cluster: " + String.join(", ", typeStrings));
    }

    protected List<String> collapseAlgorithms(final Collection<SVCallRecord> records) {
        return records.stream()
                .map(SVCallRecord::getAlgorithms)
                .flatMap(Collection::stream)
                .distinct()
                .sorted()
                .collect(Collectors.toList());
    }

    public static final class AlleleCollectionCollapserComparator implements Comparator<Collection<Allele>> {
        @Override
        public int compare(final Collection<Allele> o1, final Collection<Allele> o2) {
            if (o1 == o2) {
                return 0;
            }
            Utils.nonNull(o1);
            Utils.nonNull(o2);
            final List<Allele> sorted1 = SVCallRecordUtils.sortAlleles(o1);
            final List<Allele> sorted2 = SVCallRecordUtils.sortAlleles(o2);
            // Resort to display strings
            for (int i = 0; i < sorted1.size() && i < sorted2.size(); i++) {
                final int comp = sorted1.get(i).getDisplayString().compareTo(sorted2.get(i).getDisplayString());
                if (comp != 0) {
                    return comp;
                }
            }
            if (o1.size() < o2.size()) {
                return -1;
            } else if (o1.size() > o2.size()) {
                return 1;
            } else {
                return 0;
            }
        }
    }
}
