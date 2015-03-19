package org.broadinstitute.hellbender.engine;

import htsjdk.tribble.Feature;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.File;

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
 * @param <T> the type of Feature that this FeatureInput file contains (eg., VariantContext, BEDFeature, etc.)
 */
public final class FeatureInput<T extends Feature> {

    /**
     * Logical name for this source of Features optionally provided by the user on the command line
     * using the --argument_name logical_name:feature_file syntax. Defaults to the absolute path of
     * the underlying file if no logical name is specified
     */
    private final String name;

    /**
     * File containing Features as specified by the user on the command line
     */
    private final File featureFile;

    /**
     * Type of Feature in the featureFile. Set manually by the engine after construction.
     */
    private Class<? extends Feature> featureType;

    /**
     * Delimiter between the logical name and the file name in the --argument_name logical_name:feature_file syntax
     */
    public static final String FEATURE_ARGUMENT_TAG_DELIMITER = ":";

    /**
     * Construct a FeatureInput from a String argument value either of the form "logical_name:feature_file"
     * or simply "feature_file".
     *
     * Only meant to be called by the argument parsing system, and therefore marked as protected --
     * FeatureInputs constructed some other way will not be recognized by the engine.
     *
     * @param rawArgumentValue String of the form "logical_name:feature_file" or "feature_file"
     */
    protected FeatureInput( final String rawArgumentValue ) {
        final Pair<String, File> parsedArgument = parseRawArgumentValue(rawArgumentValue);

        name = parsedArgument.getLeft();
        featureFile = parsedArgument.getRight();
        featureType = null;  // Must be set after construction
    }

    /**
     * Parses an argument value String of the form "logical_name:feature_file" or "feature_file"
     * into logical name and file name components. The absolute path of the file is used as the
     * logical name if none is present.
     *
     * @param rawArgumentValue argument value from the command line to parse
     * @return A Pair in which the first element is the parsed logical name, and the second element
     *         is the File containing Features
     */
    private Pair<String, File> parseRawArgumentValue( final String rawArgumentValue ) {
        final String[] tokens = rawArgumentValue.split(FEATURE_ARGUMENT_TAG_DELIMITER, -1);
        final String usage = "Argument must either be a file, or of the form logical_name:file";

        // Check for malformed argument values
        if ( tokens.length > 2 || tokens.length == 0 ) {
            throw new UserException.BadArgumentValue("", rawArgumentValue, usage);
        }
        for ( String token : tokens ) {
            if ( token.isEmpty() ) {
                throw new UserException.BadArgumentValue("", rawArgumentValue, "Empty name/file encountered. " + usage);
            }
        }

        if ( tokens.length == 1 ) {
            // No user-specified logical name for this FeatureInput, so use the absolute path to the File as its name
            final File featureFile = new File(tokens[0]);
            return Pair.of(featureFile.getAbsolutePath(), featureFile);
        }
        else {
            // User specified a logical name for this FeatureInput using Name:File syntax
            return Pair.of(tokens[0], new File(tokens[1]));
        }
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
     * Gets the file backing this source of Features
     *
     * @return file backing this source of Features
     */
    public File getFeatureFile() {
        return featureFile;
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
    public void setFeatureType( Class<? extends Feature> featureType ) {
        this.featureType = featureType;
    }

    /**
     * FeatureInputs will be hashed by the engine, so make an effort to produce a reasonable hash code
     *
     * @return hash code for this FeatureInput (combination of hash codec of the name and file)
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
    public boolean equals( Object other ) {
        if ( other == null || ! (other instanceof FeatureInput) ) {
            return false;
        }

        FeatureInput<?> otherFeature = (FeatureInput<?>)other;
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
        final String featureFilePath = featureFile.getAbsolutePath();
        return name.equals(featureFilePath) ? featureFilePath :
                                              String.format("%s%s%s", name, FEATURE_ARGUMENT_TAG_DELIMITER, featureFilePath);
    }
}
