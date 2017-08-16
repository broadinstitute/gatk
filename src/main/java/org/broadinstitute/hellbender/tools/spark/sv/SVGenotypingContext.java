package org.broadinstitute.hellbender.tools.spark.sv;


import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVContext;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVHaplotype;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.AlleleList;
import org.broadinstitute.hellbender.utils.genotyper.IndexedAlleleList;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Collects relevant and convenient information for genotyping a structural variant.
 */
@DefaultSerializer(SVGenotypingContext.Serializer.class)
final class SVGenotypingContext {

    final SVContext variant;
    final List<Template> templates;
    final List<SVHaplotype> haplotypes;
    private final Map<String, List<GATKRead>> sampleTemplatesAsUnmappedReads;
    private final SampleList sampleList;
    final SVHaplotype refHaplotype;
    final SVHaplotype altHaplotype;
    final int[] refBreakPoints;
    final int[] altBreakPoints;
    final Allele refAllele;
    final Allele altAllele;
    final int refHaplotypeIndex;
    final int altHaplotypeIndex;
    private final AlleleList<Allele> alleles;
    final int refAlleleIndex;
    final int altAlleleIndex;
    final int numberOfHaplotypes;
    final int numberOfTemplates;

    /**
     * Used by kryo's deserilization.
     */
    private SVGenotypingContext(final SVContext variant, final List<SVHaplotype> haplotypes,
                                final List<Template> templates, final SampleList sampleList,
                                final AlleleList<Allele> alleles,
                                final Map<String,List<GATKRead>> sampleTemplateAsReads,
                                final SVHaplotype refHaplotype, final SVHaplotype altHaplotype,
                                final Allele refAllele, final Allele altAllele,
                                final int refHaplotypeIndex, final int altHaplotypeIndex,
                                final int refAlleleIndex, final int altAlleleIndex,
                                final int numberOfHaplotypes, final int numberOfTemplates,
                                final int[] refBreakPoints, final int[] altBreakPoints) {
        this.haplotypes = haplotypes;
        this.templates = templates;
        this.numberOfTemplates = numberOfTemplates;
        this.numberOfHaplotypes = numberOfHaplotypes;
        this.sampleTemplatesAsUnmappedReads = sampleTemplateAsReads;
        this.sampleList = sampleList;
        this.alleles = alleles;
        this.refHaplotype = refHaplotype;
        this.altHaplotype = altHaplotype;
        this.refHaplotypeIndex = refHaplotypeIndex;
        this.altHaplotypeIndex = altHaplotypeIndex;
        this.refAlleleIndex = refAlleleIndex;
        this.altAlleleIndex = altAlleleIndex;
        this.altAllele = altAllele;
        this.refAllele = refAllele;
        this.variant = variant;
        this.refBreakPoints = refBreakPoints;
        this.altBreakPoints = altBreakPoints;
    }

    static class Allele extends htsjdk.variant.variantcontext.Allele {

        private static final long serialVersionUID = 1L;

        private static Allele of(final SVHaplotype haplotype, final SVContext context) {
            if (haplotype.isNeitherReferenceNorAlternative()) {
                return new Allele(haplotype, "<" + haplotype.getName() + ">", false);
            } else if (haplotype.isReference()) {
                return new Allele(haplotype, context.getReference().getBaseString(), true);
            } else { // assume is "altHaplotype".
                return new Allele(haplotype, context.getAlternateAlleles().get(0).getDisplayString(), false);
            }
        }

        protected Allele(final SVHaplotype haplotype, final String basesString, final boolean isRef) {
            super(basesString, isRef);
        }
    }

    SVGenotypingContext(final SVContext variant,
                               final Iterable<SVHaplotype> haplotypes,
                               final Iterable<Template> templates,
                               final String sampleName,
                               final SAMSequenceDictionary dictionary) {
        Utils.nonNull(sampleName);
        this.variant = Utils.nonNull(variant);
        this.haplotypes = removeRedundantHaplotypes(Utils.stream(haplotypes).collect(Collectors.toList()));
        this.templates = Utils.stream(templates)
                .distinct()
                .filter(t -> t.fragments().stream().anyMatch(f -> f.getMappingQuality() > 0))
                .collect(Collectors.toList());
        final List<GATKRead> templatesAsUnmappedReads = this.templates.stream()
                .map(t -> t.asUnmappedRead(null)).collect(Collectors.toList());
        sampleTemplatesAsUnmappedReads = Collections.singletonMap(sampleName, templatesAsUnmappedReads);

        final Map<SVHaplotype, Integer> haplotypeToIndex = IntStream.range(0, this.haplotypes.size())
                .boxed()
                .collect(Collectors.toMap(this.haplotypes::get, i -> i));
        this.refHaplotype = this.haplotypes.stream()
                .filter(SVHaplotype::isReference)
                .findFirst()
                .orElseThrow(() -> new IllegalArgumentException("no reference haplotype provided for " + variant.getUniqueID() + "; haplotypes: " + this.haplotypes.stream().map(SVHaplotype::getName).collect(Collectors.joining(","))));
        this.altHaplotype = this.haplotypes.stream()
                .filter(SVHaplotype::isAlternative)
                .findFirst()
                .orElseThrow(() -> new IllegalArgumentException("no alternative haplotype provided"));
        this.refHaplotypeIndex = haplotypeToIndex.get(refHaplotype);
        this.altHaplotypeIndex = haplotypeToIndex.get(altHaplotype);
        this.refBreakPoints = calculateBreakPoints(refHaplotype, variant, dictionary);
        this.altBreakPoints = calculateBreakPoints(altHaplotype, variant, dictionary);
        refAllele = Allele.of(refHaplotype, variant);
        altAllele = Allele.of(altHaplotype, variant);
        alleles = new IndexedAlleleList<>(refAllele, altAllele);
        refAlleleIndex = alleles.indexOfAllele(refAllele); // would be 0
        altAlleleIndex = alleles.indexOfAllele(altAllele); // would be 1
        numberOfHaplotypes = this.haplotypes.size();
        numberOfTemplates = this.templates.size();
        sampleList = SampleList.singletonSampleList(sampleName);
    }

    private static List<SVHaplotype> removeRedundantHaplotypes(final List<SVHaplotype> in) {
        if (in.size() <= 2) {
            return in;
        } else {
            final Set<String> sequencesSoFar = new HashSet<>(in.size());
            final List<SVHaplotype> out = new ArrayList<>(in.size());
            for (final SVHaplotype haplotype : in) {
                if (haplotype.isReference() || haplotype.isAlternative()) {
                    out.add(haplotype);
                } else if (sequencesSoFar.add(new String(haplotype.getBases()))) {
                    out.add(haplotype);
                }
            }
            return in.size() == out.size() ? in : out;
        }
    }

    void reduceNumberOfTemplatesTo(final int maximum) {
        if (templates.size() > maximum) {
          final Random rdn = new Random(variant.getUniqueID().hashCode());
          Collections.shuffle(templates, rdn);
          templates.subList(maximum, templates.size()).clear();
        }
    }

    ReadLikelihoods<Allele> newLikelihoods() {
        return new ReadLikelihoods<>(sampleList, alleles, sampleTemplatesAsUnmappedReads);
    }

    private static int[] calculateBreakPoints(final SVHaplotype haplotype, final SVContext context, final SAMSequenceDictionary dictionary) {
        final List<AlignmentInterval> intervals = haplotype.getReferenceAlignment();
        final List<SimpleInterval> breakPoints = context.getBreakPointIntervals(0, dictionary, false);
        final List<Integer> result = new ArrayList<>(breakPoints.size());
        for (final SimpleInterval breakPoint : breakPoints) {
            for (final AlignmentInterval interval : intervals) {
                final Cigar cigar = interval.cigarAlongReference();
                if (interval.referenceSpan.overlaps(breakPoint)) {
                    int refPos = interval.referenceSpan.getStart();
                    int hapPos = interval.startInAssembledContig;
                    for (final CigarElement element : cigar) {
                        final CigarOperator operator = element.getOperator();
                        final int length = element.getLength();
                        if (operator.consumesReferenceBases() && breakPoint.getStart() >= refPos && breakPoint.getStart() <= refPos + length) {
                            if (operator.isAlignment()) {
                                result.add(hapPos + breakPoint.getStart() - refPos);
                            } else { // deletion.
                                result.add(hapPos);
                            }
                        } else if (!operator.consumesReferenceBases() && breakPoint.getStart() == refPos - 1) {
                            if (operator.consumesReadBases()) {
                                result.add(hapPos + length);
                            }
                        }
                        if (operator.consumesReferenceBases()) {
                            refPos += length;
                        }
                        if (operator.consumesReadBases() || operator.isClipping()) {
                            hapPos += length;
                        }
                    }
                }
            }
        }
        Collections.sort(result);
        return result.stream().mapToInt(i -> i).toArray();
    }

    public static class Serializer extends com.esotericsoftware.kryo.Serializer<SVGenotypingContext> {
        @Override
        public void write(final Kryo kryo, final Output output, final SVGenotypingContext context) {
            output.writeInt(context.numberOfHaplotypes);
            output.writeInt(context.numberOfTemplates);
            output.writeInt(context.refBreakPoints.length);
            output.writeInt(context.altBreakPoints.length);
            output.writeInt(context.refHaplotypeIndex);
            output.writeInt(context.altHaplotypeIndex);
            output.writeString(context.sampleList.getSample(0));
            kryo.writeObject(output, context.variant);
            for (final SVHaplotype haplotype : context.haplotypes) {
                kryo.writeClassAndObject(output, haplotype);
            }
            for (final Template template : context.templates) {
                kryo.writeObject(output, template);
            }
            output.writeInts(context.refBreakPoints);
            output.writeInts(context.altBreakPoints);
        }

        @Override
        public SVGenotypingContext read(Kryo kryo, Input input, Class<SVGenotypingContext> type) {
            final int numberOfHaplotypes = input.readInt();
            final int numberOfTemplates = input.readInt();
            final int numberOfReferenceBreakPoints = input.readInt();
            final int numberOfAlternativeBreakPoints = input.readInt();
            final int refHaplotypeIndex = input.readInt();
            final int altHaplotypeIndex = input.readInt();
            final String sampleName = input.readString();
            final SVContext variant = kryo.readObject(input, SVContext.class);

            final List<SVHaplotype> haplotypes = new ArrayList<>(numberOfHaplotypes);
            for (int i = 0; i < numberOfHaplotypes; i++) {
                haplotypes.add((SVHaplotype) kryo.readClassAndObject(input));
            }
            final List<Template> templates = new ArrayList<>(numberOfTemplates);
            for (int i = 0; i < numberOfTemplates; i++) {
                templates.add(kryo.readObject(input, Template.class));
            }
            final List<GATKRead> templatesAsReads = templates.stream()
                    .map(t -> t.asUnmappedRead(null))
                    .collect(Collectors.toList());
            final Map<String, List<GATKRead>> sampleTemplateAsReads = Collections.singletonMap(sampleName, templatesAsReads);
            final SampleList sampleList = SampleList.singletonSampleList(sampleName);
            final SVHaplotype refHaplotypes = haplotypes.get(refHaplotypeIndex);
            final SVHaplotype altHaplotypes = haplotypes.get(altHaplotypeIndex);
            final int[] refBreakPoints = input.readInts(numberOfReferenceBreakPoints);
            final int[] altBreakPoints = input.readInts(numberOfAlternativeBreakPoints);

            final Allele refAllele = Allele.of(refHaplotypes, variant);
            final Allele altAllele = Allele.of(altHaplotypes, variant);
            final AlleleList<Allele> alleles = new IndexedAlleleList<>(refAllele, altAllele);
            final int refAlleleIndex = alleles.indexOfAllele(refAllele);
            final int altAlleleIndex = alleles.indexOfAllele(altAllele);

            return new SVGenotypingContext(variant, haplotypes, templates, sampleList, alleles, sampleTemplateAsReads,
                    refHaplotypes, altHaplotypes, refAllele, altAllele, refHaplotypeIndex, altHaplotypeIndex,
                    refAlleleIndex, altAlleleIndex, numberOfHaplotypes, numberOfTemplates, refBreakPoints, altBreakPoints);
        }
    }
}
