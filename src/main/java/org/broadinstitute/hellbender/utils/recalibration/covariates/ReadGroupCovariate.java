package org.broadinstitute.hellbender.utils.recalibration.covariates;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.recalibration.RecalUtils;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationArgumentCollection;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.util.*;
import java.util.stream.Collectors;

/**
 * The Read Group covariate.
 */
public final class ReadGroupCovariate implements Covariate {
    private static final long serialVersionUID = 1L;

    //Note: these maps are initialized and made umodifiable at construction so the whole covariate is an immutable object once it's constructed.

    /*
     * Stores the mapping from read group id to a number.
     */
    private final Map<String, Integer> readGroupLookupTable;

    /*
     * Stores the reverse mapping, from number to read group id.
     */
    private final Map<Integer, String> readGroupReverseLookupTable;

    public ReadGroupCovariate(final RecalibrationArgumentCollection RAC, final List<String> readGroups){
        final Map<String, Integer> rgLookupTable = new LinkedHashMap<>();
        final Map<Integer, String> rgReverseLookupTable = new LinkedHashMap<>();

        readGroups.forEach(
                readGroupId -> {
                    if (!rgLookupTable.containsKey(readGroupId)) {
                        final int nextId = rgLookupTable.size();
                        rgLookupTable.put(readGroupId, nextId);
                        rgReverseLookupTable.put(nextId, readGroupId);
                    }
                }
        );
        readGroupLookupTable = Collections.unmodifiableMap(rgLookupTable);
        readGroupReverseLookupTable = Collections.unmodifiableMap(rgReverseLookupTable);
    }

    @Override
    public void recordValues(final GATKRead read, final SAMFileHeader header, final ReadCovariates values, final boolean recordIndelValues) {
        final SAMReadGroupRecord rg = ReadUtils.getSAMReadGroupRecord(read, header);
        final String readGroupId = getID(rg);
        final int key = keyForReadGroup(readGroupId);

        final int readLength = read.getLength();
        for (int i = 0; i < readLength; i++) {
            values.addCovariate(key, key, key, i);
        }
    }

    /**
     * Get the ID of the readgroup.
     */
    public static String getID(final SAMReadGroupRecord rg) {
        final String pu = rg.getPlatformUnit();
        return pu == null ? rg.getId() : pu;
    }

    @Override
    public String formatKey(final int key) {
        Utils.validate(readGroupReverseLookupTable.containsKey(key), () -> "missing key " + key);
        return readGroupReverseLookupTable.get(key);
    }

    @Override
    public int keyFromValue(final Object value) {
        return keyForReadGroup((String) value);
    }

    private int keyForReadGroup(final String readGroupId) {
        Utils.validate(readGroupLookupTable.containsKey(readGroupId), () -> "The covariates table is missing " + RecalUtils.READGROUP_COLUMN_NAME + " " + readGroupId + " in " + RecalUtils.READGROUP_REPORT_TABLE_TITLE);
        return readGroupLookupTable.get(readGroupId);
    }

    @Override
    public int maximumKeyValue() {
        return readGroupLookupTable.size() - 1;
    }

    public static List<String> getReadGroupIDs(final SAMFileHeader header) {
        return header.getReadGroups().stream().map(rg -> getID(rg)).collect(Collectors.toList());
    }
}
