package org.broadinstitute.hellbender.tools.copynumber.coverage.readcount.covariatebin;

import org.apache.commons.lang3.tuple.Pair;

import java.util.*;
import java.util.stream.Stream;

/**
 * Representation of a single multidimensional covariate bin
 *
 * @author Andrey Smirnov &lt;asmirnov@broadinstitute.org&gt;
 */
public class ReadCountCovariateBin {

    private static final String VALUE_SEPARATOR = ":";
    private static final String COVARIATE_SEPARATOR = ";";

    private final EnumMap<ReadCountCovariateBinningConfiguration, Pair<Double, Double>> covariateBinValues;

    protected ReadCountCovariateBin(final EnumMap<ReadCountCovariateBinningConfiguration, Pair<Double, Double>> covariateBinValues) {
        this.covariateBinValues = covariateBinValues;
    }

    @Override
    public boolean equals(Object o) {
        if (o == this) {
            return true;
        }
        if (o == null || this.getClass() != o.getClass()) {
            return false;
        }

        final ReadCountCovariateBin bin = (ReadCountCovariateBin) o;
        if (this.covariateBinValues.size() != bin.covariateBinValues.size()) {
            return false;
        }
        for (ReadCountCovariateBinningConfiguration config: this.covariateBinValues.keySet()) {
            if (!bin.covariateBinValues.containsKey(config)) {
                return false;
            } else {
                if (!this.covariateBinValues.get(config).equals(bin.covariateBinValues.get(config))) {
                    return false;
                }
            }

        }
        return true;
    }

    @Override
    public int hashCode() {
        int hash = 0;
        for (ReadCountCovariateBinningConfiguration config: this.covariateBinValues.keySet()) {
            hash += hash * 31 + covariateBinValues.get(config).getLeft().hashCode() + covariateBinValues.get(config).getRight().hashCode();
        }
        return hash;
    }

    //TODO document format
    @Override
    public String toString() {
        StringBuilder str = new StringBuilder();
        Stream.of(ReadCountCovariateBinningConfiguration.values()).forEach(
                config -> {
                    if(covariateBinValues.containsKey(config)) {
                        int index = config.getIndexFromBinStartAndEnd(covariateBinValues.get(config).getLeft(),
                                covariateBinValues.get(config).getRight());
                        str.append(config.getName());
                        str.append(VALUE_SEPARATOR);
                        str.append(index);
                        str.append(COVARIATE_SEPARATOR);
                    }
                }
        );
        return str.toString();
    }

}
