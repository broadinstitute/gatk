package org.broadinstitute.hellbender.tools.walkers.sv.GqRecalibrator;

import org.apache.commons.math3.util.FastMath;

import java.util.List;

class FilterSummary {
    final MinGq minGq;
    final long numMendelianTriosPassed;
    final long numMendelianTrios;
    final long numTriosPassed;
    final long numVariants;
    final long numVariantsPassed;
    final long numTruePositive;
    final long numFalsePositive;
    final long numFalseNegative;
    final long numTrueNegative;
    final boolean isLargeAlleleFrequency;
    final double truthWeight;
    final double inheritanceWeight;
    final String label;

    FilterSummary(final MinGq minGq,
                  final long numMendelianTriosPassed, final long numMendelianTrios, final long numTriosPassed,
                  final long numVariants, final long numVariantsPassed, final long numTruePositive,
                  final long numFalsePositive, final long numFalseNegative, final long numTrueNegative,
                  final boolean isLargeAlleleFrequency, final double truthWeight, final double inheritanceWeight,
                  final String label) {
        if(minGq == null) {
            throw new IllegalArgumentException("Null minGq in FilterSummary constructor");
        }
        this.minGq = minGq;
        this.numMendelianTriosPassed = numMendelianTriosPassed;
        this.numMendelianTrios = numMendelianTrios;
        this.numTriosPassed = numTriosPassed;
        this.numVariants = numVariants;
        this.numVariantsPassed = numVariantsPassed;
        this.numTruePositive = numTruePositive;
        this.numFalsePositive = numFalsePositive;
        this.numFalseNegative = numFalseNegative;
        this.numTrueNegative = numTrueNegative;
        this.isLargeAlleleFrequency = isLargeAlleleFrequency;
        this.truthWeight = truthWeight;
        this.inheritanceWeight = inheritanceWeight;
        this.label = label;
        this.check();
    }

    FilterSummary(final FilterSummary other) {
        this.minGq = other.minGq;
        this.numMendelianTriosPassed = other.numMendelianTriosPassed;
        this.numMendelianTrios = other.numMendelianTrios;
        this.numTriosPassed = other.numTriosPassed;
        this.numVariants = other.numVariants;
        this.numVariantsPassed = other.numVariantsPassed;
        this.numTruePositive = other.numTruePositive;
        this.numFalsePositive = other.numFalsePositive;
        this.numFalseNegative = other.numFalseNegative;
        this.numTrueNegative = other.numTrueNegative;
        this.isLargeAlleleFrequency = other.isLargeAlleleFrequency;
        this.truthWeight = other.truthWeight;
        this.inheritanceWeight = other.inheritanceWeight;
        this.label = other.label;
    }

    static final FilterSummary EMPTY = new FilterSummary(
            MinGq.Empty,
            0L, 0L, 0L,
            0L, 0L,
            0L, 0L, 0L, 0L,
            false, 0.0, 0.0, null
    );

    void check() {
        if(numMendelianTrios < 0) {
            throw new IllegalArgumentException("numMendelianTrios (" + numMendelianTrios + ") < 0");
        }
        if(numMendelianTriosPassed < 0) {
            throw new IllegalArgumentException("numMendelianTriosPassed  (" + numMendelianTriosPassed + ") < 0");
        }
        if(numTriosPassed < 0) {
            throw new IllegalArgumentException("numTriosPassed  (" + numTriosPassed + ") < 0");
        }
        if(numMendelianTriosPassed > numMendelianTrios) {
            throw new IllegalArgumentException("numMendelianTriosPassed (" + numMendelianTriosPassed + ") > numMendelianTrios (" + numMendelianTrios + ")");
        }
        if(numMendelianTriosPassed > numTriosPassed) {
            throw new IllegalArgumentException("numMendelianTriosPassed (" + numMendelianTriosPassed + ") > numTriosPassed (" + numTriosPassed + ")");
        }
        if(numVariantsPassed > numVariants) {
            throw new IllegalArgumentException("numVariantsPassed (" + numVariantsPassed + ") > numVariants (" + numVariants + ")");
        }
        if(numVariants < 0) {
            throw new IllegalArgumentException("numVariants (" + numVariants + ") < 0");
        }
        if(numVariantsPassed < 0) {
            throw new IllegalArgumentException("numVariantsPassed (" + numVariantsPassed + ") < 0");
        }
        if(numTruePositive < 0) {
            throw new IllegalArgumentException("numTruePositive  (" + numTruePositive + ") < 0");
        }
        if(numTrueNegative < 0) {
            throw new IllegalArgumentException("numTrueNegative  (" + numTrueNegative + ") < 0");
        }
        if(numFalseNegative < 0) {
            throw new IllegalArgumentException("numFalseNegative  (" + numFalseNegative + ") < 0");
        }
        if(numFalsePositive < 0) {
            throw new IllegalArgumentException("numFalsePositive  (" + numFalsePositive + ") < 0");
        }
    }

    boolean hasInheritanceData() { return this.numMendelianTrios > 0; }

    boolean hasOverlapData() {
        return this.numFalsePositive + this.numTruePositive + this.numFalseNegative + this.numTrueNegative > 0;
    }

    boolean isNotEmpty() {
        return hasInheritanceData() || hasOverlapData();
    }

    FilterSummary setLabel(final String newLabel) {
        return new FilterSummary(minGq, numMendelianTriosPassed, numMendelianTrios, numTriosPassed, numVariants,
                numVariantsPassed, numTruePositive, numFalsePositive, numFalseNegative, numTrueNegative,
                isLargeAlleleFrequency, truthWeight, inheritanceWeight, newLabel);
    }

    FilterSummary add(final FilterSummary other) {
        return new FilterSummary(
            minGq == null || minGq.isEmpty() ? other.minGq : minGq,
            numMendelianTriosPassed + other.numMendelianTriosPassed,
            numMendelianTrios + other.numMendelianTrios,
            numTriosPassed + other.numTriosPassed,
            numVariants + other.numVariants,
            numVariantsPassed + other.numVariantsPassed,
            numTruePositive + other.numTruePositive,
            numFalsePositive + other.numFalsePositive,
            numFalseNegative + other.numFalseNegative,
            numTrueNegative + other.numTrueNegative,
            isLargeAlleleFrequency || other.isLargeAlleleFrequency,
            FastMath.max(truthWeight, other.truthWeight),
            FastMath.max(inheritanceWeight, other.inheritanceWeight),
            label == null ? other.label : label
        );
    }

    FilterSummary subtract(final FilterSummary other) {
        return new FilterSummary(
            minGq == null || minGq.isEmpty() ? other.minGq : minGq,
            numMendelianTriosPassed - other.numMendelianTriosPassed,
            numMendelianTrios - other.numMendelianTrios,
            numTriosPassed - other.numTriosPassed,
            numVariants - other.numVariants,
            numVariantsPassed - other.numVariantsPassed,
            numTruePositive - other.numTruePositive,
            numFalsePositive - other.numFalsePositive,
            numFalseNegative - other.numFalseNegative,
            numTrueNegative - other.numTrueNegative,
            isLargeAlleleFrequency,
            FastMath.max(truthWeight, other.truthWeight),
            FastMath.max(inheritanceWeight, other.inheritanceWeight),
            label == null ? other.label : label
        );
    }

    public static String tableHeader(final List<String> labelHeader) {
        return String.join("\t", labelHeader) +
            String.format("\t%9s\t%9s\t%9s", "nMendPass", "nMendTrio", "nTrioPass") +
            String.format("\t%9s\t%9s", "nVar", "nVarPass") +
            String.format("\t%6s\t%6s\t%6s\t%6s", "nTrue+", "nFalse+", "nTrue-", "nFalse-");
    }

    public String toTableLine() {
        return String.format("%s\t%9d\t%9d\t%9d", label, numMendelianTriosPassed, numMendelianTrios, numTriosPassed)
            + String.format("\t%9d\t%9d", numVariants, numVariantsPassed )
            + String.format("\t%6d\t%6d\t%6d\t%6d", numTruePositive, numFalsePositive, numTrueNegative, numFalseNegative);
    }

    @Override
    public String toString() {
        return label + ":{minGq:" + minGq
                + ",  numMendelian/MendPassed/Passed Trios:" + numMendelianTrios + "/" + numMendelianTriosPassed + "/" + numTriosPassed
                + ", numVariantsPassed:" + numVariantsPassed + "/" + numVariants
                + ", numTruePositive:" + numTruePositive
                + ", numFalsePositive:" + numFalsePositive + ", numFalseNegative: " + numFalseNegative
                + ", numTrueNegative:" + numTrueNegative + ", largeAlleleFrequency:" + isLargeAlleleFrequency + "}";
    }
}
