package org.broadinstitute.hellbender.tools.copynumber.coverage.readcount;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.broadinstitute.barclay.utils.Utils;

import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.tools.copynumber.coverage.readcount.covariatebin.ReadCountCovariateBinningConfiguration;

import java.util.*;

/**
 * Enumeration of keys used in the read count table header to encode various parameters
 *
 * @author Andrey Smirnov &lt;asmirnov@broadinstitute.org&gt;
 */
public enum ReadCountFileHeaderKey {
    READ_COUNT_TYPE("readCountType"),
    BINNING_CONFIGURATION("covariateBinningConfiguration");

    private final String headerKeyName;

    public static final String KEY_VALUE_COMMENT_SEPARATOR = "=";
    private static final Map<String, ReadCountFileHeaderKey> stringToHeaderKeyMap = new HashMap<>();

    static {
        Arrays.asList(values()).stream().forEach(key -> stringToHeaderKeyMap.put(key.getHeaderKeyName(), key));
    }

    /**
     * @param name name of the key
     */
    ReadCountFileHeaderKey(final String name) {
        this.headerKeyName = Utils.nonNull(name);
    }

    /**
     * Get the name of the key
     *
     * never {@code null}.
     */
    public String getHeaderKeyName() {
        return headerKeyName;
    }

    /**
     * Parse the header comment line and return its key value pair if the key matches
     * one of the {@link ReadCountFileHeaderKey} instances
     *
     * @param headerComment a single header comment line
     * @return null if headerComment does not contain a known key
     */
    public static Pair<ReadCountFileHeaderKey, String> parseHeaderCommentKeyValuePair(final String headerComment) {
        final String[] splitCommentLine = headerComment.split(KEY_VALUE_COMMENT_SEPARATOR);
        if (splitCommentLine.length != 2) {
            throw new IllegalArgumentException("SAM header comment in read coverage file should have a single '" +
                    KEY_VALUE_COMMENT_SEPARATOR + "' separator.");
        }
        final ReadCountFileHeaderKey headerKey = stringToHeaderKeyMap.get(splitCommentLine[0]);
        final String value = splitCommentLine[1];
        return new ImmutablePair(headerKey, value);
    }

    /**
     * Parse the list of header comments that consist of key value pairs and return the value for a specific key
     *
     * @param headerComments list of comments from SAM header
     * @param key key to look for
     * @return a value for that key if it is found
     * @throws IllegalArgumentException if no such key is found in header comments
     */
    public static String getHeaderValueForKey(final List<String> headerComments, final ReadCountFileHeaderKey key) {
        final Optional<Pair<ReadCountFileHeaderKey, String>> optionalFileHeaderKeyValue =
                headerComments.stream().map(comment -> ReadCountFileHeaderKey.parseHeaderCommentKeyValuePair(comment))
                        .filter(pair -> pair.getKey() == key).findFirst();
        if (optionalFileHeaderKeyValue.isPresent()) {
            return optionalFileHeaderKeyValue.get().getValue();
        } else {
            throw new IllegalArgumentException("Header does not contain the key " + key.toString());
        }
    }

    /**
     *
     * @param readCountType
     * @return
     */
    public static String constructReadCountTypeCommentValue(final ReadCountType readCountType) {
        return readCountType.toString();
    }

    /**
     *
     * @param binningConfigurations
     * @return
     */
    public static String constructCovariateBinningConfigurationComment(
            final List<ReadCountCovariateBinningConfiguration> binningConfigurations) {
        final StringBuilder stringBuilder = new StringBuilder();
        binningConfigurations.stream().forEach(config -> stringBuilder.append(config.toString()));
        return stringBuilder.toString();
    }

    @Override
    public String toString() {
        return headerKeyName;
    }
}
