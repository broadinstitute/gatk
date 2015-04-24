package org.broadinstitute.hellbender.utils.report;

import java.util.EnumSet;
import java.util.HashMap;
import java.util.Map;

/**
 * The gatherable data types acceptable in a GATK report column.
 */
public enum GATKReportDataType {
    /**
     * The null type should not be used.
     */
    Null("Null"),

    /**
     * The default value when a format string is not present
     */
    Unknown("Unknown"),

    /**
     * Used for boolean values. Will display as true or false in the table.
     */
    Boolean("%[Bb]"),

    /**
     * Used for char values. Will display as a char so use printable values!
     */
    Character("%[Cc]"),

    /**
     * Used for float and double values. Will output a decimal with format %.8f unless otherwise specified.
     */
    Decimal("%.*[EeFf]"),

    /**
     * Used for int, byte, short, and long values. Will display the full number by default.
     */
    Integer("%[Dd]"),

    /**
     * Used for string values. Displays the string itself.
     */
    String("%[Ss]");

    private final java.lang.String dataTypeString;

    private GATKReportDataType(java.lang.String dataTypeString) {
        this.dataTypeString = dataTypeString;
    }

    private static final Map<java.lang.String, GATKReportDataType> lookup = new HashMap<>();

    static {
        for (GATKReportDataType s : EnumSet.allOf(GATKReportDataType.class))
            lookup.put(s.dataTypeString, s);
    }


    @Override
    public java.lang.String toString() {
        return this.dataTypeString;
    }

    /**
     * Returns a GATK report data type from the Object specified. It looks through the list of acceptable classes and
     * returns the appropriate data type.
     *
     * @param object the object ot derive the data type from
     * @return the appropriate data type
     */
    public static GATKReportDataType fromObject(Object object) {
        GATKReportDataType value;
        if (object instanceof java.lang.Boolean) {
            value = GATKReportDataType.Boolean;

        } else if (object instanceof java.lang.Character) {
            value = GATKReportDataType.Character;

        } else if (object instanceof Float ||
                object instanceof Double) {
            value = GATKReportDataType.Decimal;

        } else if (object instanceof java.lang.Integer ||
                object instanceof Long ||
                object instanceof Short ||
                object instanceof Byte) {
            value = GATKReportDataType.Integer;

        } else if (object instanceof java.lang.String) {
            value = GATKReportDataType.String;

        } else {
            value = GATKReportDataType.Unknown;
            //throw new UserException("GATKReport could not convert the data object into a GATKReportDataType. Acceptable data objects are found in the documentation.");
        }
        return value;
    }

    /**
     * Returns a GATK report data type from the format string specified. It uses regex matching from the enumerated
     * Strings.
     *
     * @param format the format string to derive the data type from
     * @return the appropriate data type
     */
    public static GATKReportDataType fromFormatString(java.lang.String format) {
        if (format.equals(""))
            return Unknown;
        for (GATKReportDataType type : lookup.values()) {
            if (format.matches(type.toString()) )
                return type;
        }
        return Unknown;
    }

    /**
     * Returns the default value of the data type. It returns an object that matches the class of the data type.
     *
     * @return an object that matches the data type
     */
    public Object getDefaultValue() {
        switch (this) {
            case Decimal:
                return 0.0D;
            case Boolean:
                return false;
            case Character:
                return '0';
            case Integer:
                return 0L;
            case String:
                return "";
            default:
                return null;
        }
    }

    /**
     * Checks if the two objects are equal using the appropriate test form the data types.
     *
     * @param a an object
     * @param b another object to check if equal
     * @return true - the objects are equal, false - the objects are nto equal
     */
    public boolean isEqual(Object a, Object b) {
        switch (this) {
            case Null:
                return true;
            case Decimal:
            case Boolean:
            case Integer:
                return a.toString().equals(b.toString());
            case Character:
            case String:
            default:
                return a.equals(b);
        }
    }

    /**
     * Converts an input String to the appropriate type using the data type. Used for parsing loading a GATK report from
     * file.
     *
     * @param obj The input string
     * @return an object that matches the data type.
     */
    Object Parse(Object obj) {
        if (obj instanceof java.lang.String) {
            java.lang.String str = obj.toString();
            switch (this) {
                case Decimal:
                    return Double.parseDouble(str);
                case Boolean:
                    return java.lang.Boolean.parseBoolean(str);
                case Integer:
                    return Long.parseLong(str);
                case String:
                    return str;
                case Character:
                    return str.toCharArray()[0];
                default:
                    return str;
            }
        } else
            return null;
    }

    /**
     * Returns a format string version of the value according to the data type.
     *
     * @return The printf string representation of the object according to data type.
     */
    public java.lang.String getDefaultFormatString() {
        switch (this) {
            case Decimal:
                return "%.8f";
            case Boolean:
                return "%b";
            case Integer:
                return "%d";
            case String:
                return "%s";
            case Character:
                return "%c";
            case Null:
            default:
                return "%s";
        }
    }
}
