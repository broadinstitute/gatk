package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.STR;

import java.util.AbstractList;
import java.util.ArrayList;
import java.util.Collections;
import java.util.ListIterator;

/**
 * An STR Model covariate.
 */
public interface STRModelCovariate {

    int indexFor(final STRContext context);

    Object valueAt(final int index);
    int size();

    static STRModelCovariate parse(final String spec) {
        if (spec.matches("^UnitLength(\\(\\d+)\\)$")) {
            return new STRUnitLength(Integer.parseInt(spec.substring(11, spec.length() - 1)));
        } else if (spec.matches("^MeanRepeatCount\\((\\d+),(\\d+)\\)$")) {
            final String[] parameters = spec.substring(16, spec.length() - 1    ).split(",");
            return new STRMeanRepeatCount(Integer.parseInt(parameters[0]), Integer.parseInt(parameters[1]));
        } else {
            throw new IllegalArgumentException("unsupported covariate spec: " + spec);
        }
    }

    final class List extends AbstractList<STRModelCovariate> {

        private java.util.List<STRModelCovariate> covariates;
        private int combinations;

        List(final java.util.List<STRModelCovariate> covariates) {
            this.covariates = new ArrayList<>(covariates);
            this.combinations = this.covariates.stream().mapToInt(STRModelCovariate::size).reduce(1, (a, b) -> a * b);
        }

        @Override
        public STRModelCovariate get(final int index) {
            if (index < 0 || index >= covariates.size()) {
                throw new IllegalArgumentException("index");
            }
            return covariates.get(index);
        }

        @Override
        public int size() {
            return covariates.size();
        }

        public Object valueAt(final int index, final STRModelCovariate covariate) {
            int combinations = 1;
            final ListIterator<STRModelCovariate> it = covariates.listIterator(covariates.size());
            while (it.hasPrevious()) {
                final STRModelCovariate previous = it.previous();
                if (previous.equals(covariate)) {
                    return covariate.valueAt((index / combinations) % previous.size());
                }
                combinations *= previous.size();
            }
            throw new IllegalArgumentException("unknown covariate: " + covariate);
        }

        public java.util.List<Object> valuesAt(final int index) {
            int combinations = 1;
            final java.util.List<Object> result = new ArrayList<>(Collections.nCopies(covariates.size(), null));
            for (int i = covariates.size() - 1; i >= 0; i--) {
                result.set(i, covariates.get(i).valueAt((index / combinations) % covariates.get(i).size()));
                combinations *= covariates.get(i).size();
            }
            return result;
        }

        public int indexFor(final STRContext context) {
            if (context == null) {
                throw new IllegalArgumentException("the context cannot be null");
            }
            int result = 0;
            int combinations = 1;
            final ListIterator<STRModelCovariate> it = covariates.listIterator(covariates.size());
            while (it.hasPrevious()) {
                final STRModelCovariate previous = it.previous();
                final int previousIndex = previous.indexFor(context);
                if (previousIndex < 0) {
                    return -1;
                }
                result = previous.indexFor(context) * combinations + result;
                combinations *= previous.size();
            }
            return result;
        }

        public int combinationsCount() {
            return combinations;
        }
    }

    final class STRUnitLength implements STRModelCovariate {
        private final int maxValue;

        STRUnitLength(final int maxUnitLength) {
            this.maxValue = maxUnitLength;
            if (maxValue < 1) {
                throw new IllegalArgumentException("the max-value cannot be less than 1");
            }
        }

        public int indexFor(final STRContext context) {
            return Math.min(maxValue - 1, context.getAlleles().getRepeatUnitLength() - 1);
        }

        @Override
        public Object valueAt(final int index) {
            if (index < 0 || index >= size()) {
                throw new IllegalArgumentException("invalid index: " + index);
            }
            return index + 1;
        }

        public int size() {
            return maxValue;
        }

        @Override
        public String toString() {
            return String.format("%s(%d)",getClass().getSimpleName().replaceFirst("^STR",""), maxValue);
        }
    }

    final class STRMeanRepeatCount implements STRModelCovariate {
        private final int maxValue;
        private final int minValue;

        STRMeanRepeatCount(final int minValue, final int maxValue) {
            this.minValue = minValue;
            this.maxValue = maxValue;
            if (minValue > maxValue) {
                throw new IllegalArgumentException("the min-value cannot be larger thant he max-value");
            } if (minValue < 1) {
                throw new IllegalArgumentException("the min-value cannot be less than 1");
            }
        }

        @Override
        public int indexFor(final STRContext context) {
            final STRAlleleSet alleles = context.getAlleles();
            final int[] depths = context.getAlleleDepths();
            int sum = 0;
            int reads = 0;
            for (int i = 0; i < alleles.size(); i++) {
                reads += depths[i];
                sum += depths[i] * alleles.get(i).repeatCount;
            }
            return Math.min(maxValue - minValue, Math.max(0, (sum + 1) / (reads + 1) - minValue));
        }

        @Override
        public Object valueAt(final int index) {
            return minValue + index;
        }

        @Override
        public int size() {
            return maxValue - minValue + 1;
        }

        @Override
        public String toString() {
            return String.format("%s(%d,%d)",getClass().getSimpleName().replaceFirst("^STR",""), minValue, maxValue);
        }
    }
}
