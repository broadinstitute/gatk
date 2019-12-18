package org.broadinstitute.hellbender.tools.sv;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.gcnv.GermlineCNVSegmentVariantComposer;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public final class SVCallRecordUtils {

    public final static List<String> nonDepthCallerAttributes = Arrays.asList(
            VCFConstants.END_KEY,
            GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE,
            GATKSVVCFConstants.STRANDS_ATTRIBUTE,
            GATKSVVCFConstants.SVLEN,
            GATKSVVCFConstants.SVTYPE
    );

    /**
     * Create a variant from a call for VCF interoperability
     *
     * @param call variant to convert
     * @return
     */
    public static VariantContextBuilder getVariantBuilder(final SVCallRecord call) {
        Utils.nonNull(call);
        final Allele altAllele = Allele.create("<" + call.getType().name() + ">", false);
        final Allele refAllele = Allele.REF_N;
        final int end;
        if (call.getType().equals(StructuralVariantType.INS) || call.getType().equals(StructuralVariantType.BND)) {
            end = call.getPositionA() + 1;
        } else {
            end = call.getPositionB();
        }
        final int end2;
        if (call.getType().equals(StructuralVariantType.INS)) {
            end2 = call.getPositionA() + 1;
        } else {
            end2 = call.getPositionB();
        }
        final VariantContextBuilder builder = new VariantContextBuilder(call.getId(), call.getContigA(), call.getPositionA(),
                end, Lists.newArrayList(refAllele, altAllele));
        builder.id(call.getId());
        builder.attribute(VCFConstants.END_KEY, end);
        builder.attribute(GATKSVVCFConstants.CONTIG2_ATTRIBUTE, call.getContigB());
        builder.attribute(GATKSVVCFConstants.END2_ATTRIBUTE, end2);
        builder.attribute(GATKSVVCFConstants.SVLEN, call.getLength());
        builder.attribute(GATKSVVCFConstants.SVTYPE, call.getType());
        builder.attribute(GATKSVVCFConstants.STRANDS_ATTRIBUTE, getStrandString(call));
        builder.attribute(GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE, call.getAlgorithms());
        builder.genotypes(call.getGenotypes());
        return builder;
    }

    public static GenotypesContext fillMissingSamplesWithEmptyGenotypes(final GenotypesContext genotypes, final Set<String> samples) {
        Utils.nonNull(genotypes);
        Utils.nonNull(samples);
        final Set<String> missingSamples = Sets.difference(samples, genotypes.getSampleNames());
        if (missingSamples.isEmpty()) {
            return genotypes;
        }
        final List<Genotype> newGenotypes = new ArrayList<>(genotypes.size() + missingSamples.size());
        newGenotypes.addAll(genotypes);
        for (final String sample : missingSamples) {
            final GenotypeBuilder genotypeBuilder = new GenotypeBuilder(sample);
            newGenotypes.add(genotypeBuilder.make());
        }
        return GenotypesContext.copy(newGenotypes);
    }

    public static SVCallRecord copyCallWithNewGenotypes(final SVCallRecord record, final GenotypesContext genotypes) {
        return new SVCallRecord(record.getId(), record.getContigA(), record.getPositionA(), record.getStrandA(), record.getContigB(),
                record.getPositionB(), record.getStrandB(), record.getType(), record.getLength(), record.getAlgorithms(),
                genotypes);
    }

    public static SVCallRecordWithEvidence copyCallWithNewGenotypes(final SVCallRecordWithEvidence record, final GenotypesContext genotypes) {
        return new SVCallRecordWithEvidence(record.getId(), record.getContigA(), record.getPositionA(), record.getStrandA(), record.getContigB(),
                record.getPositionB(), record.getStrandB(), record.getType(), record.getLength(), record.getAlgorithms(),
                genotypes, record.getStartSplitReadSites(), record.getEndSplitReadSites(), record.getDiscordantPairs(),
                record.getCopyNumberDistribution());
    }

    public static GenotypesContext filterAndAddGenotypeAttributes(final GenotypesContext genotypes,
                                                                  final Predicate<Genotype> genotypeFilter,
                                                                  final Function<Genotype, Map<String,Object>> attributeGenerator,
                                                                  final boolean clearGenotypes) {
        Utils.nonNull(genotypes);
        Utils.nonNull(genotypeFilter);
        Utils.nonNull(attributeGenerator);
        final ArrayList<Genotype> newGenotypes = new ArrayList<>();
        for (final Genotype genotype : genotypes) {
            if (genotypeFilter.test(genotype)) {
                final GenotypeBuilder builder = clearGenotypes ? new GenotypeBuilder(genotype.getSampleName()) : new GenotypeBuilder(genotype);
                builder.attributes(attributeGenerator.apply(genotype));
                newGenotypes.add(builder.make());
            }
        }
        return GenotypesContext.create(newGenotypes);
    }

    public static GenotypesContext predicateGenotypeAlleles(final GenotypesContext genotypes,
                                                            final Predicate<Genotype> predicate,
                                                            final List<Allele> trueAlleles,
                                                            final List<Allele> falseAlleles) {
        Utils.nonNull(genotypes);
        Utils.nonNull(predicate);
        Utils.nonNull(trueAlleles);
        Utils.nonNull(falseAlleles);
        final ArrayList<Genotype> newGenotypes = new ArrayList<>();
        for (final Genotype genotype : genotypes) {
            final GenotypeBuilder builder = new GenotypeBuilder(genotype);
            if (predicate.test(genotype)) {
                builder.alleles(trueAlleles);
            } else {
                builder.alleles(falseAlleles);
            }
            newGenotypes.add(builder.make());
        }
        return GenotypesContext.create(newGenotypes);
    }

    public static VariantContextBuilder createBuilderWithEvidence(final SVCallRecordWithEvidence call) {
        final VariantContextBuilder builder = getVariantBuilder(call);
        final boolean includeEvidence = !SVClusterEngine.isDepthOnlyCall(call);
        final SplitReadSite startSplitReadCounts = includeEvidence ? getSplitReadCountsAtPosition(call.getStartSplitReadSites(), call.getPositionA()) : null;
        final SplitReadSite endSplitReadCounts = includeEvidence ? getSplitReadCountsAtPosition(call.getEndSplitReadSites(), call.getPositionB()) : null;
        final Map<String,Integer> discordantPairCounts = includeEvidence ? getDiscordantPairCountsMap(call.getDiscordantPairs()) : null;
        final List<Genotype> genotypes = builder.getGenotypes();
        final List<Genotype> newGenotypes = new ArrayList<>(genotypes.size());
        for (final Genotype genotype : genotypes) {
            final String sample = genotype.getSampleName();
            final GenotypeBuilder genotypeBuilder = new GenotypeBuilder(genotype);
            if (includeEvidence) {
                final Integer startCount = startSplitReadCounts.getCount(sample);
                final Integer endCount = endSplitReadCounts.getCount(sample);
                final Integer pairedEndCount = discordantPairCounts.getOrDefault(sample, 0);
                genotypeBuilder.attribute(GATKSVVCFConstants.START_SPLIT_READ_COUNT_ATTRIBUTE, startCount);
                genotypeBuilder.attribute(GATKSVVCFConstants.END_SPLIT_READ_COUNT_ATTRIBUTE, endCount);
                genotypeBuilder.attribute(GATKSVVCFConstants.DISCORDANT_PAIR_COUNT_ATTRIBUTE, pairedEndCount);
            }
            newGenotypes.add(genotypeBuilder.make());
        }
        builder.genotypes(newGenotypes);
        return builder;
    }

    private static String getStrandString(final SVCallRecord call) {
        return getStrandString(call.getStrandA()) + getStrandString(call.getStrandB());
    }

    private static String getStrandString(final boolean strand) {
        return strand ? SVCallRecord.STRAND_PLUS : SVCallRecord.STRAND_MINUS;
    }

    private static SplitReadSite getSplitReadCountsAtPosition(final List<SplitReadSite> sites, final int pos) {
        Utils.nonNull(sites);
        Utils.validateArg(pos > 0, "Non-positive position");
        if (sites.stream().map(SplitReadSite::getPosition).distinct().count() != sites.size()) {
            throw new IllegalArgumentException("Sites did not have unique positions");
        }
        return sites.stream()
                .filter(s -> s.getPosition() == pos)
                .findAny()
                .orElse(new SplitReadSite(pos, Collections.emptyMap()));
    }

    private static Map<String,Integer> getDiscordantPairCountsMap(final Collection<DiscordantPairEvidence> discordantPairs) {
        Utils.nonNull(discordantPairs);
        return discordantPairs.stream()
                .collect(Collectors.groupingBy(DiscordantPairEvidence::getSample,
                        Collectors.reducing(0, e -> 1, Integer::sum)));
    }

    public static <T extends SVLocatable> Comparator<T> getSVLocatableComparator(final SAMSequenceDictionary dictionary) {
        return (o1, o2) -> compareSVLocatables(o1, o2, dictionary);
    }

    public static <T extends SVCallRecord> Comparator<T> getCallComparator(final SAMSequenceDictionary dictionary) {
        return (o1, o2) -> compareCalls(o1, o2, dictionary);
    }

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

    public static int compareCalls(final SVCallRecord first, final SVCallRecord second, final SAMSequenceDictionary dictionary) {
        final int compareLocatables = compareSVLocatables(first, second, dictionary);
        if (compareLocatables != 0) return compareLocatables;

        //Strands
        final int compareStartStrand = Boolean.compare(first.getStrandA(), second.getStrandA());
        if (compareStartStrand != 0) return compareStartStrand;
        final int compareEndStrand = Boolean.compare(first.getStrandB(), second.getStrandB());
        if (compareEndStrand != 0) return compareEndStrand;

        // Length
        final int compareLength = Integer.compare(first.getLength(), second.getLength());
        if (compareLength != 0) return compareLength;

        // Type
        final int compareType = first.getType().compareTo(second.getType());
        return compareType;
    }

    public static boolean isValidSize(final SVCallRecord call, final int minEventSize) {
        return call.getType().equals(StructuralVariantType.BND) || call.getLength() >= minEventSize;
    }

    public static <T> boolean intervalIsIncluded(final SVCallRecord call, final Map<String, IntervalTree<T>> includedIntervalTreeMap,
                                                 final double minDepthOnlyIncludeOverlap) {
        if (SVDepthOnlyCallDefragmenter.isDepthOnlyCall(call)) {
            return intervalIsIncludedDepthOnly(call, includedIntervalTreeMap, minDepthOnlyIncludeOverlap);
        }
        return intervalIsIncludedNonDepthOnly(call, includedIntervalTreeMap);
    }

    private static <T> boolean intervalIsIncludedNonDepthOnly(final SVCallRecord call, final Map<String,IntervalTree<T>> includedIntervalTreeMap) {
        final IntervalTree<T> startTree = includedIntervalTreeMap.get(call.getContigA());
        if (startTree == null) {
            return false;
        }
        final IntervalTree<T> endTree = includedIntervalTreeMap.get(call.getContigB());
        if (endTree == null) {
            return false;
        }
        return startTree.overlappers(call.getPositionA(), call.getPositionA() + 1).hasNext()
                && endTree.overlappers(call.getPositionB(), call.getPositionB() + 1).hasNext();
    }

    private static <T> boolean intervalIsIncludedDepthOnly(final SVCallRecord call, final Map<String,IntervalTree<T>> includedIntervalTreeMap,
                                                           final double minDepthOnlyIncludeOverlap) {
        final IntervalTree<T> tree = includedIntervalTreeMap.get(call.getContigA());
        if (tree == null) {
            return false;
        }
        final double overlapFraction = totalOverlap(call.getPositionA(), call.getPositionB(), tree) / (double) call.getLength();
        return overlapFraction >= minDepthOnlyIncludeOverlap;
    }

    private static <T> long totalOverlap(final int start, final int end, final IntervalTree<T> tree) {
        final Iterator<IntervalTree.Node<T>> iter = tree.overlappers(start, end);
        long overlap = 0;
        while (iter.hasNext()) {
            final IntervalTree.Node<T> node = iter.next();
            overlap += intersectionLength(start, end, node.getStart(), node.getEnd());
        }
        return overlap;
    }

    private static long intersectionLength(final int start1, final int end1, final int start2, final int end2) {
        return Math.max(0, Math.min(end1, end2) - Math.max(start1, start2) + 1);
    }

    public static Stream<SVCallRecord> convertInversionsToBreakends(final SVCallRecord call) {
        if (!call.getType().equals(StructuralVariantType.INV)) {
            return Stream.of(call);
        }
        Utils.validateArg(isIntrachromosomal(call), "Inversion is not intrachromosomal");
        final SVCallRecord positiveBreakend = new SVCallRecord(call.getId(), call.getContigA(),
                call.getPositionA(), true, call.getContigB(), call.getPositionB(), true, StructuralVariantType.BND, -1,
                call.getAlgorithms(), call.getGenotypes());
        final SVCallRecord negativeBreakend = new SVCallRecord(call.getId(), call.getContigA(),
                call.getPositionA(), false, call.getContigB(), call.getPositionB(), false, StructuralVariantType.BND, -1,
                call.getAlgorithms(), call.getGenotypes());
        return Stream.of(positiveBreakend, negativeBreakend);
    }

    public static void validateCoordinates(final SVCallRecord call) {
        Utils.nonNull(call);
        Utils.validateArg(call.getPositionA() >= 1, "Call start non-positive");
        if (isIntrachromosomal(call)) {
            Utils.validateArg(call.getPositionA() <= call.getPositionB(), "Second position before end on same contig");
        } else {
            Utils.validateArg(call.getPositionB() >= 1, "Call second position non-positive");
        }
    }

    public static void validateCoordinatesWithDictionary(final SVCallRecord call, final SAMSequenceDictionary dictionary) {
        Utils.nonNull(dictionary);
        validateCoordinates(call);
        final SAMSequenceRecord contigARecord = dictionary.getSequence(call.getContigA());
        Utils.validateArg(contigARecord != null, "Call first contig " + call.getContigA() + " not in dictionary");
        final SAMSequenceRecord contigBRecord = dictionary.getSequence(call.getContigB());
        Utils.validateArg(contigBRecord != null, "Call second contig " + call.getContigB() + " not in dictionary");
        Utils.validateArg(call.getPositionA() <= contigARecord.getSequenceLength(), "Call first position greater than contig length");
        Utils.validateArg(call.getPositionB() <= contigBRecord.getSequenceLength(), "Call second position greater than contig length");
    }

    public static boolean isIntrachromosomal(final SVCallRecord call) {
        return call.getContigA().equals(call.getContigB());
    }

    public static SVCallRecord create(final VariantContext variant) {
        Utils.nonNull(variant);
        Utils.validate(variant.getAttributes().keySet().containsAll(nonDepthCallerAttributes), "Call is missing attributes");
        final String id = variant.getID();
        final String contigA = variant.getContig();
        final int positionA = variant.getStart();
        final int end = variant.getEnd();
        final String contigB;
        final int positionB;

        // If END2 and CONTIG2 are both defined, use those.
        // If neither is defined, use start contig and position.
        // If only CONTIG2 is defined, END2 is taken as END
        // Having only END2 is unacceptable
        final boolean hasContig2 = variant.hasAttribute(GATKSVVCFConstants.CONTIG2_ATTRIBUTE);
        final boolean hasEnd2 = variant.hasAttribute(GATKSVVCFConstants.END2_ATTRIBUTE);
        if (hasContig2 && hasEnd2) {
            contigB = variant.getAttributeAsString(GATKSVVCFConstants.CONTIG2_ATTRIBUTE, null);
            positionB = variant.getAttributeAsInt(GATKSVVCFConstants.END2_ATTRIBUTE, 0);
        } else if (!hasContig2 && !hasEnd2) {
            contigB = contigA;
            positionB = positionA;
        } else if (hasContig2) {
            contigB = variant.getAttributeAsString(GATKSVVCFConstants.CONTIG2_ATTRIBUTE, null);
            positionB = end;
        } else {
            throw new UserException.BadInput("Attribute " + GATKSVVCFConstants.END2_ATTRIBUTE +
                    " cannot be defined without " + GATKSVVCFConstants.CONTIG2_ATTRIBUTE);
        }

        final StructuralVariantType type = variant.getStructuralVariantType();
        Utils.validateArg(variant.hasAttribute(GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE), "Attribute " + GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE + " is required");
        final List<String> algorithms = variant.getAttributeAsStringList(GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE, null);
        Utils.validateArg(variant.hasAttribute(GATKSVVCFConstants.STRANDS_ATTRIBUTE), "Attribute " + GATKSVVCFConstants.STRANDS_ATTRIBUTE + " is required");
        final String strands = variant.getAttributeAsString(GATKSVVCFConstants.STRANDS_ATTRIBUTE, null);
        if (strands.length() != 2) {
            throw new IllegalArgumentException("Strands field is not 2 characters long");
        }
        final String strand1Char = strands.substring(0, 1);
        if (!strand1Char.equals(SVCallRecord.STRAND_PLUS) && !strand1Char.equals(SVCallRecord.STRAND_MINUS)) {
            throw new IllegalArgumentException("Valid start strand not found");
        }
        final String strand2Char = strands.substring(1, 2);
        if (!strand2Char.equals(SVCallRecord.STRAND_PLUS) && !strand2Char.equals(SVCallRecord.STRAND_MINUS)) {
            throw new IllegalArgumentException("Valid end strand not found");
        }
        final boolean strand1 = strand1Char.equals(SVCallRecord.STRAND_PLUS);
        final boolean strand2 = strand2Char.equals(SVCallRecord.STRAND_PLUS);
        Utils.validateArg(variant.hasAttribute(GATKSVVCFConstants.SVLEN), "Attribute " + GATKSVVCFConstants.SVLEN + " is required");
        final int length = variant.getAttributeAsInt(GATKSVVCFConstants.SVLEN, 0);
        return new SVCallRecord(id, contigA, positionA, strand1, contigB, positionB, strand2, type, length, algorithms, variant.getGenotypes());
    }

    /**
     *
     * @param variant single-sample variant from a gCNV segments VCF
     * @param minQuality drop events with quality lower than this
     * @return
     */
    public static SVCallRecord createDepthOnlyFromGCNVWithOriginalGenotypes(final VariantContext variant, final double minQuality) {
        Utils.nonNull(variant);

        if (variant.getGenotypes().size() == 1) {
            //only cluster good variants
            final Genotype g = variant.getGenotypes().get(0);
            if (g.isHomRef() || (g.isNoCall() && !g.hasExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT))
                    || Integer.valueOf((String) g.getExtendedAttribute(GermlineCNVSegmentVariantComposer.QS)) < minQuality) {
                return null;
            }
        }


        final List<String> algorithms = Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM);

        boolean isDel = false;
        for (final Genotype g : variant.getGenotypes()) {
            if (g.isHomRef() || (g.isNoCall() && !g.hasExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT))) {
                continue;
            }
            if (variant.getReference().equals(Allele.REF_N)) {  //old segments VCFs had ref Ns and genotypes that didn't reflect ploidy accurately
                if (g.getAlleles().stream().anyMatch(a -> a.equals(GATKSVVCFConstants.DEL_ALLELE))) {
                    isDel = true;
                } else if (g.getAlleles().stream().anyMatch(a -> a.equals(GATKSVVCFConstants.DUP_ALLELE))) {
                    isDel = false;
                } else if (g.getAlleles().stream().allMatch(a -> a.isNoCall())) {
                    if (g.hasExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT)) {
                        isDel = (Integer.parseInt(g.getExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT).toString()) < g.getPloidy());
                    } else {
                        throw new IllegalStateException("Genotype for sample " + g.getSampleName() + " at " + variant.getContig() + ":" + variant.getStart() + " had no CN attribute and will be dropped.");
                    }
                } else {
                    throw new IllegalArgumentException("Segment VCF schema expects <DEL>, <DUP>, and no-call allele, but found "
                            + g.getAllele(0) + " at " + variant.getContig() + ":" + variant.getStart());
                }
            } else {  //for spec-compliant VCFs (i.e. with non-N ref allele) we can just trust the ALT
                isDel = (variant.getAlternateAlleles().contains(GATKSVVCFConstants.DEL_ALLELE)
                        && !variant.getAlternateAlleles().contains(GATKSVVCFConstants.DUP_ALLELE));
            }
        }

        final boolean startStrand = isDel ? true : false;
        final boolean endStrand = isDel ? false : true;
        final StructuralVariantType type;
        if (!variant.getReference().equals(Allele.REF_N) && variant.getAlternateAlleles().contains(GATKSVVCFConstants.DUP_ALLELE)
                && variant.getAlternateAlleles().contains(GATKSVVCFConstants.DEL_ALLELE)) {
            type = StructuralVariantType.CNV;
        } else {
            type = isDel ? StructuralVariantType.DEL : StructuralVariantType.DUP;
        }

        final String id = variant.getID();
        final String startContig = variant.getContig();
        final int start = variant.getStart();
        final int end = variant.getEnd();
        final int length = end - start;
        return new SVCallRecord(id, startContig, start, startStrand, startContig, end, endStrand, type, length, algorithms, variant.getGenotypes());
    }

    public static SVCallRecord deduplicateWithRawCallAttribute(final Collection<SVCallRecord> items) {
        if (items.isEmpty()) {
            return null;
        }
        final List<Genotype> genotypes = collapseRecordGenotypesWithRawCallAttribute(items);
        final List<String> algorithms = collapseAlgorithms(items);
        final SVCallRecord example = items.iterator().next();
        return new SVCallRecord(
                example.getId(),
                example.getContigA(),
                example.getPositionA(),
                example.getStrandA(),
                example.getContigB(),
                example.getPositionB(),
                example.getStrandB(),
                example.getType(),
                example.getLength(),
                algorithms,
                genotypes);
    }

    public static SVCallRecordWithEvidence deduplicateWithRawCallAttributeWithEvidence(final Collection<SVCallRecordWithEvidence> items) {
        if (items.isEmpty()) {
            return null;
        }
        final List<Genotype> genotypes = collapseRecordGenotypesWithRawCallAttribute(items);
        final List<String> algorithms = collapseAlgorithms(items);
        final SVCallRecordWithEvidence example = items.iterator().next();
        return new SVCallRecordWithEvidence(
                example.getId(),
                example.getContigA(),
                example.getPositionA(),
                example.getStrandA(),
                example.getContigB(),
                example.getPositionB(),
                example.getStrandB(),
                example.getType(),
                example.getLength(),
                algorithms,
                genotypes,
                example.getStartSplitReadSites(),
                example.getEndSplitReadSites(),
                example.getDiscordantPairs(),
                example.getCopyNumberDistribution());
    }

    private static List<Genotype> collapseRecordGenotypesWithRawCallAttribute(final Collection<? extends SVCallRecord> records) {
        return records.stream()
                .map(SVCallRecord::getGenotypes)
                .flatMap(g -> g.stream())
                .collect(Collectors.groupingBy(Genotype::getSampleName))
                .values()
                .stream()
                .map(SVCallRecordUtils::collapseSampleGenotypesWithRawCallAttribute)
                .collect(Collectors.toList());
    }

    private static Genotype collapseSampleGenotypesWithRawCallAttribute(final Collection<Genotype> genotypes) {
        final GenotypeBuilder builder = new GenotypeBuilder(genotypes.iterator().next().getSampleName());
        if (genotypes.stream().anyMatch(SVCallRecord::isRawCall)) {
            builder.attribute(GATKSVVCFConstants.RAW_CALL_ATTRIBUTE, GATKSVVCFConstants.RAW_CALL_ATTRIBUTE_TRUE);
        } else {
            builder.attribute(GATKSVVCFConstants.RAW_CALL_ATTRIBUTE, GATKSVVCFConstants.RAW_CALL_ATTRIBUTE_FALSE);
        }
        return builder.make();
    }

    private static List<String> collapseAlgorithms(final Collection<? extends SVCallRecord> records) {
        return records.stream()
                .map(SVCallRecord::getAlgorithms)
                .flatMap(Collection::stream)
                .distinct()
                .collect(Collectors.toList());
    }
}
