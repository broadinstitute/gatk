package org.broadinstitute.hellbender.engine;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.File;
import java.io.Serializable;
import java.net.URI;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.Map;

/**
 * Class to represent a Feature-containing input file. Tools should declare @Argument-annotated fields of
 * this type (or Collections of this type), and the Feature management system will automatically discover
 * them at runtime (provided that they are declared in the tool itself, a superclass of the tool, or an
 * ArgumentCollection of the tool).
 *
 * DO NOT ATTEMPT TO INSTANTIATE THIS CLASS DIRECTLY! FeatureInputs must be instantiated by the argument-parsing
 * system only in order to be recognized by the Feature management system. This is why the constructor is
 * marked as protected.
 *
 * FeatureInputs can be assigned logical names on the command line using the syntax:
 *
 *     --argument_name logical_name:feature_file
 *
 * These logical names can then be retrieved by the tool at runtime via {@link #getName}
 *
 * Furthermore, a list of comma-separated key=value pairs may be provided as follows:
 *
 *     --argument_name logical_name,key1=value1,key2=value2:feature_file
 *
 * the string value provided for a given key can be retrieved via {@link #getAttribute(String)}. Keys must be unique.
 *
 * @param <T> the type of Feature that this FeatureInput file contains (eg., VariantContext, BEDFeature, etc.)
 */
public final class FeatureInput<T extends Feature> implements Serializable {

    private static final long serialVersionUID = 1L;

    /**
     * Logical name for this source of Features optionally provided by the user on the command line
     * using the --argument_name logical_name:feature_file syntax. Defaults to the absolute path of
     * the underlying file if no logical name is specified
     */
    private final String name;

    private final Map<String, String> keyValueMap;

    /**
     * File containing Features as specified by the user on the command line
     */
    private final String featureFile;

    /**
     * Cache the codec for this feature input the first time we discover it, so we only do it once
     */
    private transient Class<FeatureCodec<T, ?>> featureCodecClass;

    /**
     * Delimiter between the logical name and the file name in the --argument_name logical_name:feature_file syntax
     */
    public static final String FEATURE_ARGUMENT_TAG_DELIMITER = ":";

    /**
     * Delimiter between key-value pairs in the --argument_name logical_name,key1=value1,key2=value2:feature_file syntax.
     */
    public static final String FEATURE_ARGUMENT_KEY_VALUE_PAIR_DELIMITER = ",";

    /**
     * Separator between keys and values in the --argument_name logical_name,key1=value1,key2=value2:feature_file syntax.
     */
    public static final String FEATURE_ARGUMENT_KEY_VALUE_SEPARATOR = "=";

    /**
     * Represents a parsed argument for the FeatureInput.
     * Always has a file and a name.
     * May have attributes.
     */
    private static final class ParsedArgument{
        private final Map<String, String> keyValueMap;
        private final String name;
        private final String file;
        private static final String URI_SCHEME_SEPARATOR = "//";
        public static final String USAGE = "Argument must either be a file, or of the form logical_name:file or logical_name(,key=value)*:feature_file";

        /**
         * Parses an argument value String of the forms:
         * "logical_name(,key=value)*:feature_file" or
         * "logical_name:feature_file" or
         * "feature_file"
         * into logical name and file name and key=value pairs.
         *
         * The absolute path of the file is used as the logical name if none is present.
         *
         * @param rawArgumentValue argument value from the command line to parse
         * @return The argument parsed from the provided string.
         */
        public static ParsedArgument of(final String rawArgumentValue) {
            //Use negative look ahead to avoid splitting URIs into multiple tokens
            //i.e. someName:file://somefile -> ["someName", "file://somefile"]
            final String MATCH_NAME_BUT_NOT_URI_OR_PORT = FEATURE_ARGUMENT_TAG_DELIMITER + "(?!" + URI_SCHEME_SEPARATOR + "|\\d+)";
            final String[] tokens = rawArgumentValue.split(MATCH_NAME_BUT_NOT_URI_OR_PORT, -1);

            if ( Arrays.stream(tokens).anyMatch(String::isEmpty)) {
                throw new CommandLineException.BadArgumentValue("", rawArgumentValue, "Empty name/file encountered. " + USAGE);
            }

            if (tokens.length == 0 || tokens.length > 2) {
                throw new CommandLineException.BadArgumentValue("", rawArgumentValue, USAGE);
            } else if (tokens.length == 1) {
                // No user-specified logical name for this FeatureInput, so use the absolute path to the File as its name
                final String featurePath = tokens[0];
                return new ParsedArgument(getDefaultName(featurePath), featurePath);
            } else {
                // User specified a logical name (and optional list of key-value pairs)
                // for this FeatureInput using name(,key=value)*:File syntax.
                // eg foo:file.vcf
                // eg foo,a=3,b=false,c=fred:file.vcf
                final String[] subtokens= tokens[0].split(FEATURE_ARGUMENT_KEY_VALUE_PAIR_DELIMITER, -1);
                if (subtokens[0].isEmpty()){
                    throw new CommandLineException.BadArgumentValue("", rawArgumentValue, USAGE);
                }
                final ParsedArgument pa= new ParsedArgument(subtokens[0], tokens[1]);
                //note: starting from 1 because 0 is the name
                for (int i = 1; i < subtokens.length; i++){
                    final String[] kv = subtokens[i].split(FEATURE_ARGUMENT_KEY_VALUE_SEPARATOR, -1);
                    if (kv.length != 2 || kv[0].isEmpty() || kv[1].isEmpty()){
                        throw new CommandLineException.BadArgumentValue("", rawArgumentValue, USAGE);
                    }
                    if (pa.containsKey(kv[0])){
                        throw new CommandLineException.BadArgumentValue("", rawArgumentValue, "Duplicate key " + kv[0] + "\n" + USAGE);
                    }
                    pa.addKeyValue(kv[0], kv[1]);
                }
                return pa;
            }
        }

        private static String getDefaultName(String featurePath) {
            return FeatureInput.makeIntoAbsolutePath(featurePath);
        }

        private ParsedArgument(final String name, final String file) {
            this.name=name;
            this.file=file;
            this.keyValueMap = new LinkedHashMap<>(2);
        }

        public String getFilePath(){
            return file;
        }

        public String getName() {
            return name;
        }

        /**
         * Returns an immutable view of the key-value map.
         */
        public Map<String, String> keyValueMap() {
            return Collections.unmodifiableMap(keyValueMap);
        }

        public void addKeyValue(final String k, final String v) {
            keyValueMap.put(k, v);
        }

        private boolean containsKey(final String k) {
            return keyValueMap.containsKey(k);
        }
    }

    /**
     * Construct a FeatureInput from a String argument value either of the form "logical_name:feature_file"
     * or simply "feature_file".
     *
     * Only meant to be called by the argument parsing system, and therefore marked as package-visible --
     * FeatureInputs constructed some other way will not be recognized by the engine.
     *
     * Note: cannot delegate to another constructor because Java only allows a call to "this" on the first line of a constructor.
     *
     * @param rawArgumentValue String of the form "logical_name:feature_file" or "feature_file"
     */
    FeatureInput(final String rawArgumentValue) {
        Utils.nonNull(rawArgumentValue, "rawArgumentValue");
        final ParsedArgument parsedArgument = ParsedArgument.of(rawArgumentValue);

        name = parsedArgument.getName();
        keyValueMap = parsedArgument.keyValueMap();
        this.featureFile = parsedArgument.getFilePath();
    }

    /**
     * Construct a FeatureInput from a path and a name
     *
     * This constructor is meant to be called only by the engine and test classes,
     * which is why it has package access.
     */
    FeatureInput(final String featurePath, final String name) {
        this(featurePath, name, Collections.emptyMap());
    }

    /**
     * Construct a FeatureInput from raw components: name, key value pairs and the file.
     *
     * This constructor is meant to be called by the engine and test classes --
     * FeatureInputs constructed some other way will not be recognized by the engine.
     */
    @VisibleForTesting
    public FeatureInput(final String featureFile, final String name, final Map<String, String> keyValueMap) {
        Utils.nonNull(name, "name");
        Utils.nonNull(keyValueMap, "kevValueMap");
        Utils.nonNull(featureFile, "feature-file");
        this.name = name;
        this.keyValueMap = Collections.unmodifiableMap(new LinkedHashMap<>(keyValueMap));   //make a unmodifiable copy
        this.featureFile = featureFile;
    }

    /**
     * Remember the FeatureCodec class for this input the first time it is discovered so we can bypass dynamic codec
     * discovery when multiple FeatureDataSources are created for the same input.
     */
    public void setFeatureCodecClass(final Class<FeatureCodec<T, ?>> featureCodecClass) {
        this.featureCodecClass = featureCodecClass;
    }

    /**
     * @return The previously established FeatureCodec class to use for this input, if any. May return {@code null}.
     */
    public Class<FeatureCodec<T, ?>> getFeatureCodecClass() {
        return this.featureCodecClass;
    }

    /**
     * creates a name from the given filePath by finding the absolute path of the given input
     */
    private static String makeIntoAbsolutePath(final String filePath){
        if(FeatureDataSource.isGenomicsDBPath(filePath)){
            return FeatureDataSource.GENOMIC_DB_URI_SCHEME + new File(filePath.replace(FeatureDataSource.GENOMIC_DB_URI_SCHEME,"")).getAbsolutePath();
        } else if (URI.create(filePath).getScheme() != null) {
            return IOUtils.getPath(filePath).toAbsolutePath().toUri().toString();
        } else {
            return new File(filePath).getAbsolutePath();
        }
    }

    /**
     * Gets the value for the given key associated with this Feature source or {@code null}
     * if no value is associated with a given key.
     * @throws IllegalArgumentException if the key is {@code null}.
     */
    public String getAttribute(final String key) {
        Utils.nonNull(key);
        return keyValueMap.get(key);
    }

    /**
     * Gets the logical name of this Feature source. This will be a user-provided value if the
     * --argument_name logical_name:feature_file was used on the command line, otherwise it will
     * default to the absolute path of the backing file
     *
     * @return logical name of this source of Features
     */
    public String getName() {
        return name;
    }

    /**
     *
     * @return true if the value for name does not match the default, indicating it was a user supplied name (i.e. foo:file.vcf)
     */
    public boolean hasUserSuppliedName() {
        return !ParsedArgument.getDefaultName(getFeaturePath()).equals(name);
    }

    /**
     * Gets the file backing this source of Features
     *
     * @return file backing this source of Features
     */
    public String getFeaturePath() {
        return featureFile;
    }

    /**
     * FeatureInputs will be hashed by the engine, so make an effort to produce a reasonable hash code
     *
     * @return hash code for this FeatureInput (combination of hash code of the name and file)
     */
    @Override
    public int hashCode() {
        return 31 * name.hashCode() + featureFile.hashCode();
    }

    /**
     * Test FeatureInputs for equality
     *
     * @param other object to compare this FeatureInput with
     * @return true if this FeatureInput equals other, otherwise false
     */
    @Override
    public boolean equals(final Object other) {
        if (! (other instanceof FeatureInput)) {
            return false;
        }

        final FeatureInput<?> otherFeature = (FeatureInput<?>)other;
        return name.equals(otherFeature.name) && featureFile.equals(otherFeature.featureFile);
    }

    /**
     * Returns a String representation of this FeatureInput. Will be the absolute path to
     * the featureFile if we have no logical name, or a String of the form
     * "logical_name:absolute_path_to_featureFile" if we do have a logical name.
     *
     * @return String representation of this FeatureInput
     */
    @Override
    public String toString() {
        final String featureFilePath = makeIntoAbsolutePath(featureFile);
        return name.equals(featureFilePath) ? featureFilePath :
                                              String.format("%s%s%s", name, FEATURE_ARGUMENT_TAG_DELIMITER, featureFilePath);
    }
}
