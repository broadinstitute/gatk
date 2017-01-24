package org.broadinstitute.hellbender.tools.copynumber.coverage.readcount.covariatebin;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;

import java.util.*;
import java.util.stream.Collector;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Collection of {@link ReadCountCovariateBin} objects that allows for queries given a {@link GATKRead}
 *
 * @author Andrey Smirnov &lt;asmirnov@broadinstitute.org&gt;
 */
public final class ReadCountCovariateBinCollection {

    /**
     * Ordered collection of binning configurations
     */
    private final List<ReadCountCovariateBinningConfiguration> binningConfigurations;

    /**
     * List of covariate bins in this collection
     */
    private final List<ReadCountCovariateBin> covariateBinList;

    /**
     * A map used for bookkeeping of the multidimensional representation of this collection
     */
    private final Map<ReadCountCovariateBinningConfiguration, Integer> productMap;

    /**
     * Map between covariate bin names and corresponding covariate bins
     */
    private final Map<String, ReadCountCovariateBin> nameToBinMap;

    /**
     * Public constructor
     *
     * @param binningConfigurations list of binning configurations
     */
    public ReadCountCovariateBinCollection(final List<ReadCountCovariateBinningConfiguration> binningConfigurations) {
        this.binningConfigurations = binningConfigurations;
        Utils.nonEmpty(binningConfigurations, "Binning configuration collection cannot be empty");

        productMap = new HashMap<>();
        int product = 1;
        for (ReadCountCovariateBinningConfiguration config: binningConfigurations) {
            productMap.put(config, product);
            product *= config.getNumberOfBins();
        }

        //initialize the bins
        OptionalInt numBins = binningConfigurations.stream().mapToInt(config -> config.getNumberOfBins()).reduce((a, b) -> a * b);
        covariateBinList = new ArrayList<>();
        IntStream.range(0, numBins.getAsInt()).forEach(index -> covariateBinList.add(buildNewCovariateBinInstance(this.binningConfigurations, index)));

        //initialize name to bin map
        nameToBinMap = new HashMap<>();
        covariateBinList.stream().forEach(bin -> nameToBinMap.put(bin.toString(), bin));
    }

    private ReadCountCovariateBin buildNewCovariateBinInstance(final List<ReadCountCovariateBinningConfiguration> binningConfigurations,
                                                               final int indexInArray) {
        final EnumMap<ReadCountCovariateBinningConfiguration, Pair<Double, Double>> binValues = new EnumMap<>(ReadCountCovariateBinningConfiguration.class);

        binningConfigurations.stream().forEach(config -> binValues.put(config,
                new ImmutablePair<>(config.getBinSize() * ((indexInArray / productMap.get(config)) % config.getNumberOfBins()),
                        config.getBinSize() * (((indexInArray / productMap.get(config)) % config.getNumberOfBins()) + 1))));
        return new ReadCountCovariateBin(binValues);
    }

    /**
     * Given a read return the covariate bin that it belongs to.
     *
     * @param read GATK read
     * @return bin characterized by given read
     */
    public ReadCountCovariateBin getReadCountCovariateBin(final GATKRead read) {
        int binIndex = 0;
        for (int i = binningConfigurations.size() - 1; i >= 0; i--) {
            ReadCountCovariateBinningConfiguration config = binningConfigurations.get(i);
            binIndex += config.getBinIndexFromRead(read) * productMap.get(config);
        }
        return covariateBinList.get(binIndex);
    }

    /**
     * Get list of covariate bins. Note that this list is ordered according to the multidimensional representation of the
     * bin collection
     *
     * @return list of covariate bins
     */
    public List<ReadCountCovariateBin> getCovariateBinList() {
        return covariateBinList;
    }

    /**
     * Get a {@link TableColumnCollection} with columns named according to {@link ReadCountCovariateBin#toString()} format
     *
     * @return table column collection
     */
    public TableColumnCollection getTableColumns() {
        return new TableColumnCollection(covariateBinList.stream().map(bin -> bin.toString()).collect(Collectors.toList()));
    }

    /**
     * Get a {@link ReadCountCovariateBin} object given its respective name according to the format described in
     * the {@link ReadCountCovariateBin#toString()} method
     *
     * @param binName name of the bin
     * @return {@code null} if no bin with a give name exist
     */
    public ReadCountCovariateBin getCovariateBinByName(final String binName) {
        return nameToBinMap.get(binName);
    }
}
