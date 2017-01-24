package org.broadinstitute.hellbender.tools.copynumber.coverage.readcount;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

/**
 * Enum class representing different types of coverage collection
 *
 * @author Andrey Smirnov &lt;asmirnov@broadinstitute.org&gt;
 */
public enum ReadCountType {
    SIMPLE_COUNT("SIMPLE_COUNT"),
    BINNED("BINNED");

    private final String readCountTypeName;

    static final Map<String, ReadCountType> nameToReadCountTypeMap = new HashMap<>();

    // populate the map between read count type names and their corresponding enum instances
    static {
        Arrays.stream(ReadCountType.values()).forEach(key ->
                ReadCountType.nameToReadCountTypeMap.put(key.getReadCountTypeName(), key));
    }

    ReadCountType(String name) {
        readCountTypeName = name;
    }

    /**
     * Get the read count type by its name
     *
     * @param name name
     * @return read count type, {@code null} if there is no corresponding type
     */
    public static ReadCountType getReadCountTypeByName (String name) {
        return nameToReadCountTypeMap.get(name);
    }

    public String getReadCountTypeName() {
        return readCountTypeName;
    }
}
