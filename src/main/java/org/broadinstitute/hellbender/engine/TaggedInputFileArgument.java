package org.broadinstitute.hellbender.engine;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Objects;

/**
 * Represents a file argument that has a symbolic name and may have attributes key-value mappings.
 *
 * TaggedInputFileArgument can be assigned logical names on the command line using the syntax:
 *
 *     --argument_name logical_name:file
 *
 * These logical names can then be retrieved by the tool at runtime via {@link #getName}
 *
 * Furthermore, a list of comma-separated key=value pairs may be provided as follows:
 *
 *     --argument_name logical_name,key1=value1,key2=value2:file
 *
 * the string value provided for a given key can be retrieved via {@link #getAttribute(String)}. Keys must be unique.
 */
public final class TaggedInputFileArgument {

    /**
     * Delimiter between the logical name and the file name in the --argument_name logical_name:file syntax
     */
    public static final String ARGUMENT_TAG_DELIMITER = ":";

    /**
     * Delimiter between key-value pairs in the --argument_name logical_name,key1=value1,key2=value2:file syntax.
     */
    public static final String ARGUMENT_KEY_VALUE_PAIR_DELIMITER = ",";

    /**
     * Separator between keys and values in the --argument_name logical_name,key1=value1,key2=value2:file syntax.
     */
    public static final String ARGUMENT_KEY_VALUE_SEPARATOR = "=";

    private final Map<String, String> keyValueMap;
    private final String name;
    private final File file;

    /**
     * Parses an argument value String of the forms:
     * "logical_name(,key=value)*:file" or
     * "logical_name:file" or
     * "file"
     * into logical name and file name and key=value pairs.
     *
     * The absolute path of the file is used as the logical name if none is present.
     *
     * @param rawArgumentValue argument value from the command line to parse
     */
    public TaggedInputFileArgument(final String rawArgumentValue) {
        Utils.nonNull(rawArgumentValue);
        final String[] tokens = rawArgumentValue.split(ARGUMENT_TAG_DELIMITER, -1);
        final String usage = "Argument must either be a file, or of the form logical_name:file or logical_name(,key=value)*:file";

        // Check for malformed argument values
        if (tokens.length > 2 || tokens.length == 0) {
            throw new UserException.BadArgumentValue("", rawArgumentValue, usage);
        }
        for (final String token : tokens) {
            if (token.isEmpty()) {
                throw new UserException.BadArgumentValue("", rawArgumentValue, "Empty name/file encountered. " + usage);
            }
        }

        if (tokens.length == 1) {
            // No user-specified logical name for this TaggedArgument, so use the absolute path to the File as its name
            final File file = new File(tokens[0]);
            this.name = file.getAbsolutePath();
            this.file = file;
            this.keyValueMap = Collections.emptyMap();
        } else {
            // User specified a logical name (and optional list of key-value pairs)
            // for this TaggedArgument using name(,key=value)*:File syntax.
            // eg foo:file.vcf
            // eg foo,a=3,b=false,c=fred:file.vcf
            final String[] subtokens = tokens[0].split(ARGUMENT_KEY_VALUE_PAIR_DELIMITER, -1);
            if (subtokens[0].isEmpty()) {
                throw new UserException.BadArgumentValue("", rawArgumentValue, usage);
            }
            this.name = subtokens[0];
            this.file = new File(tokens[1]);
            final Map<String, String> keyValueMap = new LinkedHashMap<>(2);
            //note: starting from 1 because 0 is the name
            for (int i = 1; i < subtokens.length; i++) {
                final String[] kv = subtokens[i].split(ARGUMENT_KEY_VALUE_SEPARATOR, -1);
                if (kv.length != 2 || kv[0].isEmpty() || kv[1].isEmpty()) {
                    throw new UserException.BadArgumentValue("", rawArgumentValue, usage);
                }
                if (keyValueMap.containsKey(kv[0])) {
                    throw new UserException.BadArgumentValue("", rawArgumentValue, "Duplicate key " + kv[0] + "\n" + usage);
                }
                keyValueMap.put(kv[0], kv[1]);
            }
            this.keyValueMap = Collections.unmodifiableMap(keyValueMap);
        }
    }

    /**
     * Create a new argument from the given arguments.
     * No argument may be null.
     */
    public TaggedInputFileArgument(final String name, final Map<String, String> keyValueMap, final File file) {
        this.name = Utils.nonNull(name, "name");
        this.file = Utils.nonNull(file, "file");
        Utils.nonNull(keyValueMap, "keyValueMap");
        this.keyValueMap = Collections.unmodifiableMap(new LinkedHashMap<>(keyValueMap));//make an unmodifiable copy
    }

    /**
     * Returns the file represented by this argument.
     */
    public File getFile(){
        return file;
    }

    /**
     * Gets the path of the file represented by this argument.
     */
    public String getFilePath() {
        return file.getAbsolutePath();
    }

    /**
     * Returns the symbolic name of this argument.
     * If no explicit symbolic name was given when the argument was created,
     * then this is the name of the file that this argument represents.
     */
    public String getName() {
        return name;
    }

    /**
     * Gets the value for the given key associated with this argument or {@code null}
     * if no value is associated with a given key.
     * @throws IllegalArgumentException if the key is {@code null}.
     */
    public String getAttribute(final String key) {
        Utils.nonNull(key);
        return keyValueMap.get(key);
    }

    /**
     * Returns the key value pairings in this argument. The returned map may be empty and
     * it is always unmodifiable.
     */
    public Map<String, String> getAttributes() {
        return keyValueMap;
    }

    /**
     * Returns whether this input has a specific symbolic name that is different than the default.
     */
    public boolean hasSymbolicName(){
        return !Objects.equals(name, getFilePath());
    }

    /**
     * Returns a String representation of this ReadInput. Will be the absolute path to
     * the readFile if we have no logical name, or a String of the form
     * "logical_name:absolute_path_to_readFile" if we do have a logical name.
     *
     * @return String representation of this ReadInput
     */
    @Override
    public String toString() {
        final String readFilePath = file.getAbsolutePath();
        return !hasSymbolicName() ? readFilePath :
                String.format("%s%s%s", getName(), TaggedInputFileArgument.ARGUMENT_TAG_DELIMITER, readFilePath);
    }

    @Override
    public int hashCode() {
        //key-value pairs do not count for hashCode
        return Objects.hash(name, file);
    }

    @Override
    public boolean equals(final Object obj) {
        if (obj == this){
            return true;
        }
        if (!(obj instanceof TaggedInputFileArgument)){
            return false;
        }
        //key-value pairs do not count for equality
        final TaggedInputFileArgument that = (TaggedInputFileArgument)obj;
        return Objects.equals(this.name, that.name) &&  Objects.equals(this.file, that.file);
    }
}
