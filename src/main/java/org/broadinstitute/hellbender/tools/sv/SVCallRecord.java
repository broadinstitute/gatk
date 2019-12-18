package org.broadinstitute.hellbender.tools.sv;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeType;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.StructuralVariantType;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.VariantContextGetters;

import java.util.*;
import java.util.stream.Collectors;

public class SVCallRecord implements SVLocatable {

    public static final String STRAND_PLUS = "+";
    public static final String STRAND_MINUS = "-";

    private final String id;
    private final String contigA;
    private final int positionA;
    private final boolean strandA;
    private final String contigB;
    private final int positionB;
    private final boolean strandB;
    private final StructuralVariantType type;
    private int length;
    private final List<String> algorithms;
    private final GenotypesContext genotypes;

    public SVCallRecord(final String id,
                        final String contigA,
                        final int positionA,
                        final boolean strandA,
                        final String contigB,
                        final int positionB,
                        final boolean strandB,
                        final StructuralVariantType type,
                        final int length,
                        final List<String> algorithms,
                        final List<Genotype> genotypes) {
        Utils.nonNull(id);
        Utils.nonNull(contigA);
        Utils.nonNull(contigB);
        Utils.nonNull(type);
        Utils.nonNull(algorithms);
        Utils.nonNull(genotypes);
        Utils.containsNoNull(algorithms, "Encountered null algorithm");
        Utils.containsNoNull(genotypes, "Encountered null genotype");
        this.id = id;
        this.contigA = contigA;
        this.positionA = positionA;
        this.strandA = strandA;
        this.contigB = contigB;
        this.positionB = positionB;
        this.strandB = strandB;
        this.type = type;
        this.length = length;
        this.algorithms = Collections.unmodifiableList(algorithms);
        this.genotypes = GenotypesContext.copy(genotypes).immutable();
    }

    public String getId() {
        return id;
    }

    @Override
    public String getContigA() {
        return contigA;
    }

    @Override
    public int getPositionA() {
        return positionA;
    }

    @Override
    public String getContigB() {
        return contigB;
    }

    @Override
    public int getPositionB() {
        return positionB;
    }

    public boolean getStrandA() {
        return strandA;
    }

    public boolean getStrandB() {
        return strandB;
    }

    public StructuralVariantType getType() {
        return type;
    }

    public int getLength() {
        return length;
    }

    public List<String> getAlgorithms() {
        return algorithms;
    }

    public Set<String> getAllSamples() {
        return genotypes.stream().map(Genotype::getSampleName)
                .collect(Collectors.toCollection(LinkedHashSet::new));
    }

    public Set<String> getCalledSamples() {
        return genotypes.stream()
                .filter(g -> !g.getType().equals(GenotypeType.NO_CALL) || g.hasExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT)) //skip no-calls that don't have copy number (i.e. not dupe no-calls)
                .map(Genotype::getSampleName)
                .collect(Collectors.toCollection(LinkedHashSet::new));
    }

    public Set<String> getCarrierSamples() {
        return genotypes.stream()
                .filter(SVCallRecord::isCarrier)
                .map(Genotype::getSampleName)
                .collect(Collectors.toCollection(LinkedHashSet::new));
    }

    public static boolean isCarrier(final Genotype g) {
        Utils.nonNull(g);
        return g.getType().equals(GenotypeType.HET) || g.getType().equals(GenotypeType.HOM_VAR);
    }

    public static boolean isRawCall(final Genotype g) {
        Utils.nonNull(g);
        return VariantContextGetters.getAttributeAsInt(g, GATKSVVCFConstants.RAW_CALL_ATTRIBUTE, GATKSVVCFConstants.RAW_CALL_ATTRIBUTE_FALSE) == GATKSVVCFConstants.RAW_CALL_ATTRIBUTE_TRUE;
    }

    public GenotypesContext getGenotypes() {
        return genotypes;
    }

    public SimpleInterval getPositionAInterval() {
        return new SimpleInterval(contigA, positionA, positionA);
    }

    public SimpleInterval getPositionBInterval() {
        return new SimpleInterval(contigB, positionB, positionB);
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof SVCallRecord)) return false;
        SVCallRecord that = (SVCallRecord) o;
        return positionA == that.positionA &&
                strandA == that.strandA &&
                positionB == that.positionB &&
                strandB == that.strandB &&
                length == that.length &&
                id.equals(that.id) &&
                contigA.equals(that.contigA) &&
                contigB.equals(that.contigB) &&
                type == that.type &&
                algorithms.equals(that.algorithms) &&
                genotypes.size() == that.genotypes.size() &&
                genotypes.containsAll(that.genotypes);  // TODO : List comparison expensive like this
    }

    @Override
    public int hashCode() {
        return Objects.hash(id, contigA, positionA, strandA, contigB, positionB, strandB, type, length, algorithms, genotypes);
    }
}
