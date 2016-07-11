package org.broadinstitute.hellbender.engine;

import htsjdk.tribble.Feature;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.util.Map;
import java.util.Objects;

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
public final class FeatureInput<T extends Feature> {

    /**
     * File with a symbolic name and key-value mappings.
     */
    private final TaggedInputFileArgument taggedInputFileArgument;

    /**
     * Type of Feature in the featureFile. Set manually by the engine after construction.
     */
    private Class<? extends Feature> featureType;

    /**
     * Construct a FeatureInput from a String argument value either of the form "logical_name:feature_file"
     * or simply "feature_file", or "logical_name(,key=value)*:feature_file".
     *
     * Only meant to be called by the argument parsing system, and therefore marked as package-visible --
     * FeatureInputs constructed some other way will not be recognized by the engine.
     *
     * Note: cannot delegate to another constructor because Java only allows a call to "this" on the first line of a constructor.
     *
     * @param rawArgumentValue String of the form "logical_name:feature_file", "feature_file" or "logical_name(,key=value)*:feature_file".
     */
    FeatureInput(final String rawArgumentValue) {
        Utils.nonNull(rawArgumentValue, "rawArgumentValue");
        taggedInputFileArgument = new TaggedInputFileArgument(rawArgumentValue);
        featureType = null;  // Must be set after construction
    }

    /**
     * Construct a FeatureInput from raw components: name, key value pairs and the file.
     *
     * This constructor is meant to be called by the engine and test classes --
     * FeatureInputs constructed some other way will not be recognized by the engine.
     */
    public FeatureInput(final String name, final Map<String, String> keyValueMap, final File featureFile) {
        Utils.nonNull(name, "name");
        Utils.nonNull(keyValueMap, "keyValueMap");
        Utils.nonNull(featureFile, "featureFile");
        this.taggedInputFileArgument = new TaggedInputFileArgument(name, keyValueMap, featureFile);
        this.featureType = null;  // Must be set after construction
    }

    /**
     * Gets the value for the given key associated with this Feature source or {@code null}
     * if no value is associated with a given key.
     * @throws IllegalArgumentException if the key is {@code null}.
     */
    public String getAttribute(final String key) {
        Utils.nonNull(key);
        return taggedInputFileArgument.getAttribute(key);
    }

    /**
     * Gets the logical name of this Feature source. This will be a user-provided value if the
     * --argument_name logical_name:feature_file was used on the command line, otherwise it will
     * default to the absolute path of the backing file
     *
     * @return logical name of this source of Features
     */
    public String getName() {
        return taggedInputFileArgument.getName();
    }

    /**
     * Gets the file backing this source of Features
     *
     * @return file backing this source of Features
     */
    public File getFeatureFile() {
        return taggedInputFileArgument.getFile();
    }

    /**
     * Gets the type of Feature contained in our file
     *
     * @return the type of Feature contained in our file
     */
    public Class<? extends Feature> getFeatureType() {
        return featureType;
    }

    /**
     * Sets the type of Feature contained in our file. Called by the engine after construction time.
     *
     * @param featureType the type of Feature contained in our file
     */
    protected void setFeatureType(final Class<? extends Feature> featureType) {
        this.featureType = featureType;
    }

    /**
     * FeatureInputs will be hashed by the engine, so make an effort to produce a reasonable hash code
     *
     * @return hash code for this FeatureInput (combination of hash code of the name and file)
     */
    @Override
    public int hashCode() {
        //Note: the featureType does not count for hashCode
        return Objects.hashCode(taggedInputFileArgument);
    }

    /**
     * Test FeatureInputs for equality
     *
     * @param other object to compare this FeatureInput with
     * @return true if this FeatureInput equals other, otherwise false
     */
    @Override
    public boolean equals(final Object other) {
        if (other == this){
            return true;
        }
        if (!(other instanceof FeatureInput)){
            return false;
        }

        final FeatureInput<?> that = (FeatureInput<?>)other;
        //Note: the featureType does not count for equality
        return Objects.equals(this.taggedInputFileArgument, that.taggedInputFileArgument);
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
        return taggedInputFileArgument.toString();
    }
}
