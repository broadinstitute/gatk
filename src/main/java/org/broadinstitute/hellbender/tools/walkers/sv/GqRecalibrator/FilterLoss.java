package org.broadinstitute.hellbender.tools.walkers.sv.GqRecalibrator;

import org.apache.commons.math3.util.FastMath;
import org.jetbrains.annotations.NotNull;

import java.util.Collection;
import java.util.List;
import java.util.stream.Collectors;

class FilterLoss implements Comparable<FilterLoss> {
    final double inheritanceLoss;
    final double truthLoss;
    final double inheritanceWeight;
    final double truthWeight;
    final String label;
    static final double targetPrecision = MinGqVariantFilterBase.targetPrecision;
    static final int minVariantsToEstimateRecall = MinGqVariantFilterBase.minVariantsToEstimateRecall;

    FilterLoss(final double inheritanceLoss, final double truthLoss,
               final double inheritanceWeight, final double truthWeight, final String label) {
        if(inheritanceWeight < 0) {
            throw new IllegalArgumentException("inheritanceWeight (" + inheritanceWeight + ") < 0");
        }
        if(truthWeight < 0) {
            throw new IllegalArgumentException("truthWeight (" + truthWeight + ") < 0");
        }
        this.inheritanceLoss = inheritanceLoss;
        this.truthLoss = truthLoss;
        this.inheritanceWeight = inheritanceWeight;
        this.truthWeight = truthWeight;
        this.label = label;
    }


    FilterLoss(final FilterSummary filterSummary) {
        this(1.0 - getInheritanceF1(filterSummary),
             1.0 - getTruthF1(filterSummary),
             filterSummary.inheritanceWeight * inheritanceAlleleFrequencyWeight(filterSummary),
             filterSummary.truthWeight,
             filterSummary.toTableLine());
    }

    FilterLoss(final FilterLoss copyLoss, final String label) {
        this(copyLoss.inheritanceLoss, copyLoss.truthLoss, copyLoss.inheritanceWeight, copyLoss.truthWeight, label);
    }

    FilterLoss(final Collection<FilterSummary> filterSummaries, final List<String> labelColumns) {
        this(
            filterSummaries.stream()
                .map(FilterLoss::new)
                .reduce(FilterLoss::add).orElse(EMPTY),
            getHeader(labelColumns) + "\n" +
                filterSummaries.stream()
                    .filter(FilterSummary::isNotEmpty)
                    .map(FilterLoss::new)
                    .map(FilterLoss::toString)
                    .collect(Collectors.joining("\n"))
                    + "\n" + filterSummaries.stream()
                        .filter(FilterSummary::isNotEmpty)
                        .reduce(FilterSummary::add).orElse(FilterSummary.EMPTY)
                        .setLabel(getOverallLabel(labelColumns))
                        .toTableLine()
        );
    }

    protected static String getHeader(final List<String> labelColumns){
        return FilterSummary.tableHeader(labelColumns) + "\tinheritLoss\ttruthLoss";
    }

    protected static String getOverallLabel(final List<String> labelColumns) {
        return labelColumns.stream()
            .map(col -> MinGqVariantFilterBase.padWidth("all", col.length()))
            .collect(Collectors.joining("\t"));
    }

    protected static double inheritanceAlleleFrequencyWeight(final FilterSummary filterSummary) {
        return filterSummary.isLargeAlleleFrequency ? 1e-4 : 1.0;
    }

    protected static double getInheritanceF1(final FilterSummary filterSummary) {
        // calculate f1 score:
        //     recall = numPassedAlleleCount / maxPassedAlleleCount
        //     precision = numMendelian / numDiscoverableTrios
        //     f1 = 2.0 / (1.0 / recall + 1.0 / precision)
        // No data to filter? -> no loss -> f1 = 1.0
        if(filterSummary.numMendelianTrios == 0) {
            // No way to score this if no Mendelian trios are discoverable
            return Double.NaN;
        } else if(filterSummary.numTriosPassed == 0) {
            return 0.0;  // Discovered nothing -> get bad score if you passed nothing
        }
        final double precision = filterSummary.numMendelianTriosPassed / (double)filterSummary.numTriosPassed;
        if(precision > targetPrecision) { // in high-precision regime, maximize f1
            final double recall = filterSummary.numMendelianTriosPassed / (double)filterSummary.numMendelianTrios;
            final double f1 = 2.0 / (1.0 / recall + 1.0 / precision);
            return f1 * (1.0 - targetPrecision) + targetPrecision;  // join continuously to low-precision regime
        } else {
            return precision; // in low-precision regime, just maximize precision
        }
    }

    protected static double getTruthF1(final FilterSummary filterSummary) {
        // NOTE: calls are compared against list of true variants and false variants
        // calculate f1 score:
        //     f1 = 2.0 / (1.0 / recall + 1.0 / precision)
        //     recall = numTruePositive / (numTruePositive + numFalseNegative)
        //     precision = numTruePositive / (numTruePositive + numFalsePositive)
        //     -> f1 = 2.0 * numTruePositive / (2 * numTruePositive + numFalseNegative + numFalsePositive)
        // No data to filter? -> no loss -> f1 = 1.0
        final long numPassedTruthKnown = filterSummary.numTruePositive + filterSummary.numFalsePositive;

        if(numPassedTruthKnown == 0) {
            // no truth information known, don't try to score
            return Double.NaN;
        }
        final double precision = filterSummary.numTruePositive / (double)numPassedTruthKnown;
        if(precision > targetPrecision) { // in high-precision regime, maximize f1
            final long numActuallyTrue = filterSummary.numTruePositive + filterSummary.numFalseNegative;
            // try to use only truth data to estimate recall, but if there isn't enough, fall back to using fraction
            // of variants passed
            final double recall = numActuallyTrue >= minVariantsToEstimateRecall ?
                filterSummary.numTruePositive / (double)numActuallyTrue :
                filterSummary.numVariantsPassed / (double)filterSummary.numVariants;

            final double f1 = 2.0 / (1.0 / recall + 1.0 / precision);
            return f1 * (1.0 - targetPrecision) + targetPrecision;  // join continuously to low-precision regime
        } else {
            return precision; // in low-precision regime, just maximize precision
        }
    }


    static FilterLoss add(final FilterLoss lossA, final FilterLoss lossB) {
        return new FilterLoss(addLosses(lossA.inheritanceLoss, lossA.inheritanceWeight,
                                        lossB.inheritanceLoss, lossB.inheritanceWeight),
                              addLosses(lossA.truthLoss, lossA.truthWeight,
                                        lossB.truthLoss, lossB.truthWeight),
                              lossA.inheritanceWeight + lossB.inheritanceWeight,
                              lossA.truthWeight + lossB.truthWeight,
                              lossA.label == null ? lossB.label : lossA.label);
    }

    static private double addLosses(final double lossA, final double weightA, final double lossB, final double weightB) {
        final double weight = (weightA + weightB);
        if(weight <= 0) {
            return Double.NaN;
        }
        if(Double.isFinite(lossA)) {
            if(Double.isFinite(lossB)) {
                return (lossA * weightA + lossB * weightB) / weight;
            } else {
                return lossA * weightA / weight;
            }
        } else {
            return lossB * weightB / weight;
        }
    }

    // NaN is a crap result, so treat it accordingly
    public boolean ge(final FilterLoss other) {
        return toDouble() >= other.toDouble();
    }

    public boolean gt(final FilterLoss other) {
        return toDouble() > other.toDouble();
    }

    public boolean le(final FilterLoss other) {
        return toDouble() <= other.toDouble();
    }

    public boolean lt(final FilterLoss other) {
        return toDouble() < other.toDouble();
    }

    @Override
    public int compareTo(final @NotNull FilterLoss other) {
        final double thisVal = this.toDouble();
        final double  otherVal = other.toDouble();
        if(Double.isNaN(thisVal)) {
            return Double.isNaN(otherVal) ? 0 : 1;
        } else if(Double.isNaN(otherVal)) {
            return -1;
        } else {
            return Double.compare(thisVal, otherVal);
        }
    }

    private static String getDescription(final double inheritanceLoss, final double truthLoss) {
        return String.format("%11.9f\t%11.9f", inheritanceLoss, truthLoss);
    }

    @Override
    public String toString() {
        return label == null ?
                getDescription(inheritanceLoss, truthLoss) :
                label + "\t" + getDescription(inheritanceLoss, truthLoss);
    }

    double toDouble() {
        final double weightedInheritanceLoss = Double.isFinite(inheritanceLoss) ?
            inheritanceWeight * inheritanceLoss :
            inheritanceWeight;
        final double weightedTruthLoss = Double.isFinite(truthLoss) ?
            truthWeight * truthLoss :
            truthWeight;
        return FastMath.sqrt(weightedInheritanceLoss * weightedInheritanceLoss +
                                 weightedTruthLoss * weightedTruthLoss);
    }

    float toFloat() {
        return (float) toDouble();
    }

    static final FilterLoss EMPTY = new FilterLoss(Double.NaN, Double.NaN, 0.0, 0.0, null);
}
