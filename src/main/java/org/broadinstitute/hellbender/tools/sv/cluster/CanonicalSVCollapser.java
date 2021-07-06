package org.broadinstitute.hellbender.tools.sv.cluster;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.primitives.Doubles;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.variant.variantcontext.*;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVCallRecordUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
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

    private final BreakpointSummaryStrategy breakpointSummaryStrategy;
    private final ReferenceSequenceFile reference;

    private static final AlleleCollectionCollapserComparator ALLELE_COMPARATOR = new AlleleCollectionCollapserComparator();

    public CanonicalSVCollapser(final ReferenceSequenceFile reference, final BreakpointSummaryStrategy breakpointSummaryStrategy) {
        this.reference = reference;
        this.breakpointSummaryStrategy = Utils.nonNull(breakpointSummaryStrategy);
    }

    @Override
    public SVCallRecord collapse(final Collection<SVCallRecord> items) {
        final String id = collapseIds(items);
        final List<String> algorithms = collapseAlgorithms(items);
        final StructuralVariantType type = collapseTypes(items);
        final Map<String, Object> attributes = collapseVariantAttributes(items);

        // Prefer using variants generated with PESR callers, which tend to generate more precise breakpoints
        final Collection<SVCallRecord> mostPreciseCalls = getMostPreciseCalls(items);
        final SVCallRecord exampleCall = mostPreciseCalls.iterator().next();

        final Map.Entry<Integer, Integer> coordinates = collapseInterval(mostPreciseCalls);
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
     * Finds representative ploidy for a sample. Note the ploidy of each genotype is determined by the number of alleles.
     * The maximum ploidy is returned, as in most cases this is likely to be the best result.
     * @param sampleGenotypes genotypes belonging to a single sample
     * @return ploidy
     */
    protected int collapsePloidy(final Collection<Genotype> sampleGenotypes) {
        Utils.nonNull(sampleGenotypes);
        Utils.nonEmpty(sampleGenotypes);
        return sampleGenotypes.stream().mapToInt(Genotype::getPloidy).max().getAsInt();
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
    protected static final String getSymbolicAlleleBaseSymbol(final Allele allele) {
        return allele.getDisplayString()
                .replace("<", "")
                .replace(">", "")
                .split(":")[0];
    }

    /**
     * Collapses alternate alleles a list of representative alleles. Note this supports sub-typed alleles such as
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

            // Check that base symbols match
            final List<String> uniqueAlleleBaseSymbols = altAlleles.stream()
                    .map(CanonicalSVCollapser::getSymbolicAlleleBaseSymbol)
                    .distinct()
                    .collect(Collectors.toList());
            Utils.validate(uniqueAlleleBaseSymbols.size() == 1, "Encountered multiple symbolic allele base symbols for non-CNV");

            // Look for subtyped alts
            final List<Allele> subtypedAlleles = altAlleles.stream().filter(a -> a.getDisplayString().contains(":")).collect(Collectors.toList());
            if ( subtypedAlleles.size() == 1) {
                return Collections.singletonList(subtypedAlleles.get(0));
            } else {
                final String baseSymbol = uniqueAlleleBaseSymbols.get(0);
                return Collections.singletonList(Allele.create("<" + baseSymbol + ">", false));
            }
        }
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

    @VisibleForTesting
    protected List<Allele> collapseSampleAlleles(final Collection<Genotype> genotypes,
                                                 final Allele refAllele,
                                                 final List<Allele> sampleAltAlleles) {
        Utils.nonNull(genotypes);
        Utils.nonEmpty(genotypes);
        final int inferredPloidy = collapsePloidy(genotypes);
        if (inferredPloidy == 0) {
            return Collections.emptyList();
        }
        if (sampleAltAlleles.isEmpty()) {
            return Collections.nCopies(inferredPloidy, refAllele);
        }

        final Map<List<Allele>, Integer> genotypeCounts = genotypes.stream().map(Genotype::getAlleles)
                .collect(Collectors.groupingBy(l -> l, Collectors.collectingAndThen(Collectors.toList(), List::size)));
        List<Allele> bestGenotypeAltAlleles = null;
        int bestGenotypeFrequency = 0;
        // Sort keys somehow, for stability
        final List<List<Allele>> sortedAlleleLists = genotypeCounts.keySet().stream().sorted(ALLELE_COMPARATOR).collect(Collectors.toList());
        // TODO : break ties
        for (final List<Allele> alleles : sortedAlleleLists) {
            final List<Allele> altAlleles = alleles.stream().filter(SVCallRecordUtils::isAltAllele).collect(Collectors.toList());
            final int genotypeFrequency = genotypeCounts.get(alleles);
            if (altAlleles.size() > 0 && (bestGenotypeAltAlleles == null || genotypeFrequency > bestGenotypeFrequency)) {
                bestGenotypeAltAlleles = altAlleles;
                bestGenotypeFrequency = genotypeFrequency;
            }
        }

        final List<Allele> alleles = new ArrayList<>(inferredPloidy);
        final int numCollapsedRefAlleles = inferredPloidy - bestGenotypeAltAlleles.size();
        Utils.validate(numCollapsedRefAlleles >= 0, "Ploidy is less than number of alt alleles");
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
                .filter(a -> a != null && !a.isNoCall() && !a.isReference())
                .distinct()
                .collect(Collectors.toList());
        builder.alleles(collapseSampleAlleles(genotypes, refAllele, sampleAltAlleles));
        builder.noAttributes();
        builder.attributes(collapseGenotypeAttributes(genotypes));
        return builder.make();
    }

    protected Map<String, Object> collapseGenotypeAttributes(final Collection<Genotype> genotypes) {
        Utils.nonNull(genotypes);
        Utils.nonEmpty(genotypes);
        final Map<String, Object> collapsedAttributes = new HashMap<>();
        final Map<String, Set<Object>> genotypeFields = genotypes.stream().map(Genotype::getExtendedAttributes)
                .map(Map::entrySet)
                .flatMap(Set::stream)
                .collect(Collectors.groupingBy(Map.Entry::getKey, Collectors.mapping(Map.Entry::getValue, Collectors.toSet())));
        for (final Map.Entry<String, Set<Object>> entry : genotypeFields.entrySet()) {
            collapsedAttributes.put(entry.getKey(), collapseSampleGenotypeAttribute(entry.getValue()));
        }
        return collapsedAttributes;
    }

    protected Object collapseSampleGenotypeAttribute(final Set<Object> values) {
        if (values.size() == 1) {
            return values.iterator().next();
        } else {
            return null;
        }
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
            // Median
            final int[] lengths = items.stream().mapToInt(SVCallRecord::getLength).toArray();
            final int midIndex = lengths.length / 2;
            if (lengths.length % 2 == 0) {
                return (int) Math.ceil((lengths[midIndex - 1] + lengths[midIndex]) / 2.0);
            } else {
                return lengths[midIndex];
            }
        } else {
            return newEnd - newStart + 1;
        }
    }

    protected String collapseIds(final Collection<SVCallRecord> records) {
        Utils.nonNull(records);
        Utils.nonEmpty(records);
        return records.stream().map(SVCallRecord::getId).sorted().collect(Collectors.toList()).get(0);
    }

    protected Collection<SVCallRecord> getMostPreciseCalls(final Collection<SVCallRecord> items) {
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
    protected Map.Entry<Integer,Integer> collapseInterval(final Collection<SVCallRecord> items) {
        Utils.nonNull(items);
        Utils.nonEmpty(items);
        final SVCallRecord exampleCall = items.iterator().next();
        if (items.size() > 1) {
            final List<String> contigsA = items.stream().map(SVCallRecord::getContigA).distinct().collect(Collectors.toList());
            Utils.validate(contigsA.size() == 1, "Cannot collapse intervals with multiple position A contigs");
            final List<String> contigsB = items.stream().map(SVCallRecord::getContigB).distinct().collect(Collectors.toList());
            Utils.validate(contigsB.size() == 1, "Cannot collapse intervals with multiple position B contigs");
        }
        final List<Integer> startPositions = items.stream().map(SVCallRecord::getPositionA).sorted().collect(Collectors.toList());
        final List<Integer> endPositions = items.stream().map(SVCallRecord::getPositionB).sorted().collect(Collectors.toList());
        //use the mid value of the sorted list so the start and end represent real breakpoint observations
        final int medianStart = startPositions.get(startPositions.size() / 2);
        final int medianEnd = endPositions.get(endPositions.size() / 2);
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
                    newStart = startPositions.get(0);
                    newEnd = endPositions.get(endPositions.size() - 1);
                    break;
                case MAX_START_MIN_END:
                    newStart = startPositions.get(startPositions.size() - 1);
                    newEnd = endPositions.get(0);
                    break;
                case MEAN_START_MEAN_END:
                    newStart = (int)Math.round(new Mean().evaluate(Doubles.toArray(startPositions)));
                    newEnd = (int)Math.round(new Mean().evaluate(Doubles.toArray(endPositions)));
                    break;
                default:
                    throw new UnsupportedOperationException("Unknown breakpoint summary strategy: " + breakpointSummaryStrategy.name());
            }
        }
        return new AbstractMap.SimpleImmutableEntry<>(newStart, newEnd);
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
