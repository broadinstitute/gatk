package org.broadinstitute.hellbender.cmdline.parser;

import htsjdk.samtools.util.StringUtil;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;


final class ArgumentDefinition {
    private static final int MINIMUM_VISUAL_SEPARATION = 2;
    private final Field field;
    private final String fieldName;
    private final String fullName;
    private final String shortName;
    private final String doc;
    private final String defaultValue;
    private final boolean isCommon;
    private final Set<String> mutuallyExclusive;
    private final Object parent;
    private final boolean isSensitive;
    private final boolean isCollection;

    private final boolean optional;


    private boolean hasBeenSet = false;

    public ArgumentDefinition(final Field field, final Argument annotation, final Object parent) {
        this.field = field;
        this.fieldName = field.getName();
        this.parent = parent;
        this.fullName = annotation.fullName();
        this.shortName = annotation.shortName();
        this.doc = annotation.doc();
        this.isCollection = CommandLineParser.isCollectionField(field);

        this.isCommon = annotation.common();
        this.isSensitive = annotation.sensitive();

        this.mutuallyExclusive = new HashSet<>(Arrays.asList(annotation.mutex()));

        this.defaultValue = getDefaultValueString(isCollection, getFieldValue());

        //null collections have been initialized by createCollection which is called in handleArgumentAnnotation
        //this is optional if it's specified as being optional or if there is a default value specified
        this.optional = annotation.optional() || !this.defaultValue.equals(CommandLineParser.NULL_STRING);
    }

    private static String getDefaultValueString(boolean isCollection, Object defaultValue) {
        if (defaultValue != null) {
            if (isCollection && ((Collection) defaultValue).isEmpty()) {
                //treat empty collections the same as uninitialized primitive types
                return CommandLineParser.NULL_STRING;
            } else {
                //this is an initialized primitive type or a non-empty collection
                return defaultValue.toString();
            }
        } else {
            return CommandLineParser.NULL_STRING;
        }
    }

    public String getArgumentParamUsage(final Map<String, ArgumentDefinition> argumentMap) {
        final StringBuilder output = new StringBuilder();

        final String argumentLabel = getPrettyNameAndType();
        output.append(argumentLabel);

        int numSpaces = CommandLineParser.ARGUMENT_COLUMN_WIDTH - argumentLabel.length();
        if (argumentLabel.length() > CommandLineParser.ARGUMENT_COLUMN_WIDTH - MINIMUM_VISUAL_SEPARATION) {
            output.append(System.lineSeparator());
            numSpaces = CommandLineParser.ARGUMENT_COLUMN_WIDTH;
        }
        output.append(StringUtils.repeat(' ', numSpaces));
        final String wrappedDescription = StringUtil.wordWrap(makeArgumentDescription(argumentMap), CommandLineParser.DESCRIPTION_COLUMN_WIDTH);
        final String[] descriptionLines = wrappedDescription.split("\n");

        final String offset = System.lineSeparator() + StringUtils.repeat(' ', CommandLineParser.ARGUMENT_COLUMN_WIDTH);
        output.append(String.join(offset, Arrays.asList(descriptionLines)));
        return output.toString();
    }

    public String getPrettyNameAndType() {
        final StringBuilder output = new StringBuilder();
        output.append("--").append(getLongName());

        if (!shortName.isEmpty() && !shortName.equals(getLongName())) {
            output.append(",-").append(shortName);
        }
        output.append(':').append(CommandLineParser.getUnderlyingType(field).getSimpleName());
        if (isCollection) {
            output.append("[]");
        }
        return output.toString();
    }

    private String makeArgumentDescription(Map<String, ArgumentDefinition> argumentMap) {
        final StringBuilder sb = new StringBuilder();
        if (!doc.isEmpty()) {
            sb.append(doc);
            sb.append("  ");
        }
        if (optional) {
            sb.append("Default: ");
            sb.append(defaultValue);
            sb.append(". ");
        }
        sb.append(CommandLineParser.getOptions(CommandLineParser.getUnderlyingType(field)));
        if (!mutuallyExclusive.isEmpty()) {
            sb.append(" Cannot be used with argument(s)");
            for (final String argument : mutuallyExclusive) {
                final ArgumentDefinition mutextArgumentDefinition = argumentMap.get(argument);

                if (mutextArgumentDefinition == null) {
                    throw new GATKException("Invalid argument definition in source code.  " + argument +
                            " doesn't match any known argument.");
                }

                sb.append(' ').append(mutextArgumentDefinition.fieldName);
                if (!mutextArgumentDefinition.shortName.isEmpty()) {
                    sb.append(" (").append(mutextArgumentDefinition.shortName).append(')');
                }
            }
        }
        return sb.toString();
    }

    @SuppressWarnings("unchecked")
    /**
     * Sets the value of the field associated to this argument definition
     */
    public void setArgument(final List<String> values) {
        //special treatment for flags
        if (isFlag() && values.isEmpty()) {
            hasBeenSet = true;
            setFieldValue(true);
            return;
        }

        if (!isCollection && (hasBeenSet || values.size() > 1)) {
            throw new UserException.CommandLineException("Argument '" + getNames() + "' cannot be specified more than once.");
        }

        for (String stringValue : values) {
            final Object value;
            if (stringValue.equals(CommandLineParser.NULL_STRING)) {
                //"null" is a special value that allows the user to override any default
                //value set for this arg
                if (optional) {
                    value = null;
                } else {
                    throw new UserException.CommandLineException("Non \"null\" value must be provided for '" + getNames() + "'.");
                }
            } else {
                value = CommandLineParser.constructFromString(CommandLineParser.getUnderlyingType(field), stringValue, getLongName());
            }

            if (isCollection) {
                @SuppressWarnings("rawtypes")
                final Collection c = (Collection) getFieldValue();
                if (value == null) {
                    //user specified this arg=null which is interpreted as empty list
                    c.clear();
                } else {
                    c.add(value);
                }
                hasBeenSet = true;
            } else {
                setFieldValue(value);
                hasBeenSet = true;
            }
        }
    }


    public Object getFieldValue() {
        try {
            field.setAccessible(true);
            return field.get(parent);
        } catch (IllegalAccessException e) {
            throw new GATKException.ShouldNeverReachHereException("This shouldn't happen since we setAccessible(true).", e);
        }
    }

    public void setFieldValue(final Object value) {
        try {
            field.setAccessible(true);
            field.set(parent, value);
        } catch (IllegalAccessException e) {
            throw new GATKException.ShouldNeverReachHereException("BUG: couldn't set field value. For "
                    + fieldName + " in " + parent + " with value " + value
                    + " This shouldn't happen since we setAccessible(true)", e);
        }
    }

    public boolean isFlag() {
        return field.getType().equals(boolean.class) || field.getType().equals(Boolean.class);
    }

    public List<String> getNames() {
        final List<String> names = new ArrayList<>();
        if (!shortName.isEmpty()) {
            names.add(shortName);
        }
        if (!fullName.isEmpty()) {
            names.add(fullName);
        } else {
            names.add(fieldName);
        }
        return names;
    }

    public String getLongName() {
        return !fullName.isEmpty() ? fullName : fieldName;
    }

    /**
     * Helper for pretty printing this option.
     *
     * @param value A value this argument was given
     * @return a string
     */
    private String prettyNameValue(Object value) {
        if (value != null) {
            if (isSensitive) {
                return String.format("--%s ***********", getLongName());
            } else {
                return String.format("--%s %s", getLongName(), value);
            }
        }
        return "";
    }

    /**
     * @return A string representation of this argument and it's value(s) which would be valid if copied and pasted
     * back as a command line argument
     */
    public String toCommandLineString() {
        final Object value = getFieldValue();
        if (this.isCollection) {
            final Collection<?> collect = (Collection<?>) value;
            return collect.stream()
                    .map(this::prettyNameValue)
                    .collect(Collectors.joining(" "));

        } else {
            return prettyNameValue(value);
        }
    }

    /**
     * @param argumentDefinitionMap a map from argument names to {@link ArgumentDefinition}
     * @return the set of arguments that this argument is mutually exclusive with (including itself in the set)
     */
    public Set<ArgumentDefinition> getMutuallyExclusiveArguments(Map<String, ArgumentDefinition> argumentDefinitionMap) {
        final Set<ArgumentDefinition> mutexArguments = mutuallyExclusive.stream().map(argumentDefinitionMap::get).collect(Collectors.toSet());
        mutexArguments.add(this);
        return mutexArguments;
    }

    public String composeMissingArgumentMessage(Map<String, ArgumentDefinition> argumentDefinitionMap) {
        if (mutuallyExclusive.isEmpty()) {
            return "the following required argument was missing:\n" + getArgumentParamUsage(argumentDefinitionMap) + '\n';
        } else {
            final Set<ArgumentDefinition> mutuallyExclusiveArguments = getMutuallyExclusiveArguments(argumentDefinitionMap);
            return mutuallyExclusiveArguments.stream()
                    .map(arg -> arg.getArgumentParamUsage(argumentDefinitionMap))
                    .sorted()
                    .collect(Collectors.joining("\n", "exactly one of the following arguments is required:\n", "\n"));
        }
    }

    public boolean isOptional() {
        return optional;
    }

    public boolean hasBeenSet() {
        return hasBeenSet;
    }

    public String getDoc() {
        return doc;
    }

    public String getDefaultValue() {
        return defaultValue;
    }

    public boolean isCommon() {
        return isCommon;
    }

    public Set<String> getMutuallyExclusive() {
        return mutuallyExclusive;
    }

}
