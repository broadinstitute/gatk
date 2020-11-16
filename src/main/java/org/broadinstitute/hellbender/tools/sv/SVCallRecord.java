package org.broadinstitute.hellbender.tools.sv;

import htsjdk.tribble.Feature;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.StructuralVariantType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.tools.copynumber.gcnv.GermlineCNVSegmentVariantComposer;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.codecs.SVCallRecordCodec;

import java.util.*;
import java.util.stream.Collectors;

public class SVCallRecord implements Feature {

    private final String startContig;
    private final int start;
    private final boolean startStrand;
    private final String endContig;
    private final int end;
    private final boolean endStrand;
    private final StructuralVariantType type;
    private int length;
    private final List<String> algorithms;
    private final List<Genotype> genotypes;
    private LinkedHashSet<String> samples;

    private final static List<String> nonDepthCallerAttributes = Arrays.asList(
            VCFConstants.END_KEY,
            GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE,
            GATKSVVCFConstants.STRANDS_ATTRIBUTE,
            GATKSVVCFConstants.SVLEN,
            GATKSVVCFConstants.SVTYPE
    );

    public static SVCallRecord create(final VariantContext variant) {
        Utils.nonNull(variant);
        Utils.validate(variant.getAttributes().keySet().containsAll(nonDepthCallerAttributes), "Call is missing attributes");
        final String startContig = variant.getContig();
        final int start = variant.getStart();
        final String endContig = variant.getAttributeAsString(GATKSVVCFConstants.END_CONTIG_ATTRIBUTE, "NA");
        final int end = variant.getAttributeAsInt(VCFConstants.END_KEY, variant.getStart());
        final StructuralVariantType type = variant.getStructuralVariantType();
        final List<String> algorithms = variant.getAttributeAsStringList(GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE, GATKSVVCFConstants.DEPTH_ALGORITHM);
        final String strands = variant.getAttributeAsString(GATKSVVCFConstants.STRANDS_ATTRIBUTE, "0");
        if (strands.length() != 2) {
            throw new IllegalArgumentException("Strands field is not 2 characters long");
        }
        final String startStrandChar = strands.substring(0, 1);
        if (!startStrandChar.equals(SVCallRecordCodec.STRAND_PLUS) && !startStrandChar.equals(SVCallRecordCodec.STRAND_MINUS)) {
            throw new IllegalArgumentException("Valid start strand not found");
        }
        final String endStrandChar = strands.substring(1, 2);
        if (!endStrandChar.equals(SVCallRecordCodec.STRAND_PLUS) && !endStrandChar.equals(SVCallRecordCodec.STRAND_MINUS)) {
            throw new IllegalArgumentException("Valid end strand not found");
        }
        final boolean startStrand = startStrandChar.equals(SVCallRecordCodec.STRAND_PLUS);
        final boolean endStrand = endStrandChar.equals(SVCallRecordCodec.STRAND_PLUS);
        final int length = variant.getAttributeAsInt(GATKSVVCFConstants.SVLEN, 0);
        return new SVCallRecord(startContig, start, startStrand, endContig, end, endStrand, type, length, algorithms, variant.getGenotypes());
    }

    /**
     *
     * @param variant single-sample variant from a gCNV segments VCF
     * @param minQuality drop events with quality lower than this
     * @return a new SVCallRecord, may be null
     */
    public static SVCallRecord createDepthOnlyFromGCNV(final VariantContext variant, final double minQuality) {
        Utils.nonNull(variant);

        if (variant.getGenotypes().size() == 1) {
            //only cluster good variants
            final Genotype g = variant.getGenotypes().get(0);
            if (g.isHomRef() || (g.isNoCall() && !g.hasExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT))
                    || Integer.valueOf((String) g.getExtendedAttribute(GermlineCNVSegmentVariantComposer.QS)) < minQuality
                    || isNullCall(g)) {
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

        final String startContig = variant.getContig();
        final int start = variant.getStart();
        final int end = variant.getEnd();
        final int length = end - start;
        return new SVCallRecord(startContig, start, startStrand, startContig, end, endStrand, type, length, algorithms,
                new ArrayList<>(variant.getGenotypes()));
    }

    /**
     *
     * @param g
     * @return true if this is a call on a missing contig
     */
    private static boolean isNullCall(final Genotype g) {
        return g.hasExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT)
                && Integer.parseInt(g.getExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT).toString()) == 0
                && g.isNoCall();

    }

    public SVCallRecord(final String startContig,
                        final int start,
                        final boolean startStrand,
                        final String endContig,
                        final int end,
                        final boolean endStrand,
                        final StructuralVariantType type,
                        final int length,
                        final List<String> algorithms,
                        final List<Genotype> genotypes) {
        Utils.nonNull(startContig);
        Utils.nonNull(endContig);
        Utils.nonNull(type);
        Utils.nonNull(algorithms);
        Utils.nonNull(genotypes);
        Utils.nonEmpty(algorithms);
        Utils.nonEmpty(genotypes);
        Utils.containsNoNull(algorithms, "Encountered null algorithm");
        Utils.containsNoNull(genotypes, "Encountered null genotype");
        this.startContig = startContig;
        this.start = start;
        this.startStrand = startStrand;
        this.endContig = endContig;
        this.end = end;
        this.endStrand = endStrand;
        this.type = type;
        this.length = length;
        this.algorithms = algorithms;
        this.genotypes = genotypes;
        this.samples = genotypes.stream()
                .filter(g -> g.hasExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT)) //skip no-calls that don't have copy number (i.e. not dupe no-calls)
                .map(Genotype::getSampleName)
                .collect(Collectors.toCollection(LinkedHashSet::new));
    }

    @Override
    public String getContig() {
        return startContig;
    }

    @Override
    public int getStart() {
        return start;
    }

    public boolean getStartStrand() {
        return startStrand;
    }

    public String getEndContig() {
        return endContig;
    }

    @Override
    public int getEnd() {
        return end;
    }

    public boolean getEndStrand() {
        return endStrand;
    }

    public StructuralVariantType getType() {
        return type;
    }

    public boolean isCNV() {
        return type == StructuralVariantType.CNV || type == StructuralVariantType.DEL || type == StructuralVariantType.DUP;
    }

    public int getLength() {
        return length;
    }

    public List<String> getAlgorithms() {
        return algorithms;
    }

    public Set<String> getSamples() {
        return samples;
    }

    public List<Genotype> getGenotypes() {
        return genotypes;
    }

    public SimpleInterval getStartAsInterval() {
        return new SimpleInterval(startContig, start, start);
    }

    public SimpleInterval getEndAsInterval() {
        return new SimpleInterval(endContig, end, end);
    }

    @Override
    public boolean equals(final Object obj) {
        if (obj == null) {
            return false;
        }
        if (this.getClass() != obj.getClass()) {
            return false;
        }
        final SVCallRecord b = (SVCallRecord) obj;

        //quick check
        if (!this.getContig().equals(b.getContig())) return false;
        if (this.getStart() != b.getStart()) return false;

        boolean areEqual = this.getStartStrand() == b.getStartStrand();

        areEqual &= this.getEndContig() == b.getEndContig();
        areEqual &= this.getEnd() == b.getEnd();
        areEqual &= this.getEndStrand() == b.getEndStrand();

        areEqual &= this.getType() == b.getType();
        areEqual &= this.getLength() == b.getLength();

        areEqual &= this.getAlgorithms().containsAll(b.getAlgorithms());
        areEqual &= b.getAlgorithms().containsAll(this.getAlgorithms());

        areEqual &= this.getSamples().containsAll(b.getSamples());
        areEqual &= b.getSamples().containsAll(this.getSamples());

        areEqual &= this.getGenotypes().containsAll(b.getGenotypes());
        areEqual &= b.getGenotypes().containsAll(this.getGenotypes());

        return areEqual;
    }

    @Override
    public int hashCode() {
        return Objects.hash(algorithms, end, endContig, endStrand, genotypes,
                length, samples, start, startContig, startStrand, type);
    }
}
