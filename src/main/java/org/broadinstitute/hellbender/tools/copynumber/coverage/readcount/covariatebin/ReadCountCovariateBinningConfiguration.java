package org.broadinstitute.hellbender.tools.copynumber.coverage.readcount.covariatebin;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Function;
import java.util.stream.Collectors;

/**
 * Configuration of a binning layout for a single covariate of read count coverage. Enum instances
 * are meant to represent different read count coverage covariates (such as GC or fragment length).
 *
 * Note: this implementation assumes that bins for a single covariate are spaced uniformly
 *
 * @author Andrey Smirnov &lt;asmirnov@broadinstitute.org&gt;
 */
public enum ReadCountCovariateBinningConfiguration implements ReadCovariateValueEvaluator {
    /**
     * Binning configuration for GC covariate.
     */
    GC_CONTENT("GC") {
        @Override
        public double getValueFromRead(GATKRead read) {
            final Nucleotide.Counter counter = new Nucleotide.Counter();
            counter.addAll(read.getBases());
            final long gcCount = counter.get(Nucleotide.C) + counter.get(Nucleotide.G);
            final long atCount = counter.get(Nucleotide.A) + counter.get(Nucleotide.T);
            final long totalCount = gcCount + atCount;
            Utils.validate(totalCount != 0, "Total number of bases in the read cannot be zero");
            return gcCount / (double) totalCount;
        }
    },
    /**
     * Binning configuration for fragment length covariate.
     */
    FRAGMENT_LENGTH("FL") {
        @Override
        public double getValueFromRead(GATKRead read) {
            //take absolute value since call to fragment length method of the read can return a negative value that indicates direction
            return Math.abs(read.getFragmentLength());
        }
    };

    private double startValue;
    private double endValue;
    private int numberOfBins;
    private double binSize;
    private final String name;

    private static final String VALUE_SEPARATOR = ":";
    private static final String COVARIATE_SEPARATOR = ";";
    private static final String VALUE_DOUBLE_FORMAT = "%.2f";

    private static final Map<String, ReadCountCovariateBinningConfiguration> nameToCovariateTypeMap = new HashMap<>();

    static {
        Arrays.asList(values()).stream().forEach(type -> nameToCovariateTypeMap.put(type.getName(), type));
    }

    /**
     * {@link ReadCountCovariateBinningConfiguration} constructor
     *
     * @param name configuration's string identifier
     */
    ReadCountCovariateBinningConfiguration(String name) {
        this.name = name;
    }

    /**
     * Get the starting value of the binning
     *
     * @return starting value
     */
    public double getStartingValue() {
        return startValue;
    }

    /**
     * Get the ending value of the binning
     *
     * @return ending value
     */
    public double getEndingValue() {
        return endValue;
    }

    /**
     * Get the number of bins in the binning configuration. Note that bins are assumed to be equally spaced
     *
     * @return number of bins
     */
    public int getNumberOfBins() {
        return numberOfBins;
    }

    /**
     * Helper method that computes the index of the bin in this configuration given a read
     *
     * @param read GATK read
     * @return index of the read
     */
    public int getBinIndexFromRead(final GATKRead read) {
        return getBinIndexFromValue(getValueFromRead(read));
    }

    /**
     * Helper method to get the size of each bin for this binning configuration
     *
     * @return size of the bin
     */
    public double getBinSize() {
        return binSize;
    }

    /**
     * Helper method that computes index of the bin given bin's start and end values
     *
     * @param binStartValue start value
     * @param binEndValue end value
     * @return index of the bin
     */
    public int getIndexFromBinStartAndEnd(final double binStartValue, final double binEndValue) {
        Utils.validateArg(binStartValue < binEndValue, "Bin start values has to be strictly less than end value");
        Utils.validateArg(binStartValue >= 0., "Bin start value should be positive");
        return getBinIndexFromValue((binStartValue + binEndValue) / 2.);
    }

    private int getBinIndexFromValue(final double value) {
        return Math.max(Math.min((int) ((value - getStartingValue()) / getBinSize()), getNumberOfBins() - 1), 0);
    }

    /**
     * Get string identifier of the covariate type
     *
     * @return identifier string
     */
    public String getName() {
        return name;
    }

    /**
     * Set and validate parameters of this binning configuration
     *
     * @param numberOfBins number of bins
     * @param startingValue start value
     * @param endingValue end value
     * @throws IllegalArgumentException if this combination of arguments is not legal
     */
    public ReadCountCovariateBinningConfiguration setParameters(final double startingValue, final double endingValue, final int numberOfBins) {
        Utils.validateArg(numberOfBins > 0, "Number of bins in the binning configuration has to be positive");
        Utils.validateArg(startingValue >= 0.0, "Start value of the binning configuration cannot be negative");
        Utils.validateArg(startingValue < endingValue, "Start value of the binning configuration has to be strictly less than end value");
        this.numberOfBins = numberOfBins;
        this.startValue = startingValue;
        this.endValue = endingValue;
        this.binSize = (endingValue - startingValue) / (double) numberOfBins;
        return this;
    }

    /**
     * Parse a string representation of a list of covariate configurations
     *
     * @param parameterString string representation
     * @return list of covariate configurations encoded in string representation in the same order
     * as they appear in the string
     */
    public static List<ReadCountCovariateBinningConfiguration> parseParameters(final String parameterString) {
        return Arrays.asList(parameterString.split(COVARIATE_SEPARATOR)).stream()
                .map(ReadCountCovariateBinningConfiguration::fromString).collect(Collectors.toList());
    }

    /**
     * Parse a string representation of a single covariate configuration
     *
     * @param parameterString
     * @return
     * @throws UserException.BadInput if parameter string is not formatted correctly
     */
    public static ReadCountCovariateBinningConfiguration fromString(final String parameterString) {
        String trimmedString = parameterString.replace(COVARIATE_SEPARATOR, "");
        List<String> values = Arrays.asList(trimmedString.split(VALUE_SEPARATOR));
        if (values.size() != 4) {
            throw new UserException.BadInput(String.format("The covariate bin descriptor %s is incorrectly formatted" +
                    "The correct format is 'NAME%sSTART%sSTOP%s'", parameterString, VALUE_SEPARATOR, VALUE_SEPARATOR, COVARIATE_SEPARATOR));
        }
        String covariateNameKey = values.get(0);
        double startPosition = Double.parseDouble(values.get(1));
        double endPosition = Double.parseDouble(values.get(2));
        int numberOfBins = Integer.parseInt(values.get(3));

        ReadCountCovariateBinningConfiguration config = ReadCountCovariateBinningConfiguration.getConfigurationByName(covariateNameKey);
        if (config == null) {
            throw new UserException.BadInput(String.format("Binning configuration with type %s doesn't exist", covariateNameKey));
        }
        config.setParameters(startPosition, endPosition, numberOfBins);
        return config;
    }

    /**
     * Get {@link ReadCountCovariateBinningConfiguration} enum instance given its name identifier
     *
     * @param name covariate configuration name
     * @return covariate binning configuration
     */
    public static ReadCountCovariateBinningConfiguration getConfigurationByName(String name) {
        return nameToCovariateTypeMap.get(name);
    }

    //TODO document toString format
    @Override
    public String toString() {
        StringBuilder str = new StringBuilder();
        str.append(this.name).append(VALUE_SEPARATOR);
        str.append(String.format(VALUE_DOUBLE_FORMAT, this.startValue)).append(VALUE_SEPARATOR);
        str.append(String.format(VALUE_DOUBLE_FORMAT, this.endValue)).append(VALUE_SEPARATOR);
        str.append(this.numberOfBins).append(COVARIATE_SEPARATOR);
        return str.toString();
    }

}
