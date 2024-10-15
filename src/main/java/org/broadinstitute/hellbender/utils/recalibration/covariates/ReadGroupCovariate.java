package org.broadinstitute.hellbender.utils.recalibration.covariates;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import org.broadinstitute.hellbender.utils.Utils;
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
    private boolean allowReadGroupsNotInRecalTable; // tsato: rename

    public static final int MISSING_READ_GROUP_KEY = -1;

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
        final Map<Integer, String> rgReverseLookupTable = new LinkedHashMap<>(); // tsato: this should be an array
        // allowReadGroupsNotInRecalTable = RAC.allowReadGroupsNotInRecalTable; // tsato: or should I just keep all of RAC...
        allowReadGroupsNotInRecalTable = false;

                readGroups.forEach(
                readGroupId -> {
                    if (!rgLookupTable.containsKey(readGroupId)) { // tsato: nextid can be 0
                        final int nextId = rgLookupTable.size(); // tsato: ok I guess this works too...
                        rgLookupTable.put(readGroupId, nextId); // tsato: setting the breakpoint here to see if the first id is 0
                        rgReverseLookupTable.put(nextId, readGroupId);
                    }
                }
        );
        readGroupLookupTable = Collections.unmodifiableMap(rgLookupTable);
        readGroupReverseLookupTable = Collections.unmodifiableMap(rgReverseLookupTable);
    }

    @Override
    public void recordValues(final GATKRead read, final SAMFileHeader header, final PerReadCovariateMatrix covariateTable, final boolean recordIndelValues) {
        final SAMReadGroupRecord rg = ReadUtils.getSAMReadGroupRecord(read, header);
        final String readGroupId = getReadGroupIdentifier(rg); // tsato: really this is a read group "name". ID is also standard.

        final int key = keyForReadGroup(readGroupId); // tsato: key = index, rename it.
        // tsato: I guess the ReadGroupCovariate instance can now whether it's called by ... recalibrate or apply
        // tsato: by the way, BaseRecalibrator should be called, CollectBQSRData or something, then ApplyBQSR.
        Utils.validate(key >= -1, "key must be -1 or a nonnegative integer but is " + key); // tsato: can it be 0...

        final int readLength = read.getLength();
        for (int i = 0; i < readLength; i++) {
            covariateTable.addCovariate(key, key, key, i); // tsato: (mismatch, insertion, deletion) for the first three args
        }
    }

    /**
     * If present, we use the Platform Unit (PU) as the identifier of a read group, rather than the read group ID.
     * PU has the format {FLOWCELL_BARCODE}.{LANE}.{SAMPLE_BARCODE}.
     *
     * See https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups
     *
     */
    public static String getReadGroupIdentifier(final SAMReadGroupRecord rg) {
        final String pu = rg.getPlatformUnit();
        return pu == null ? rg.getId() : pu;
    }


    @Override
    public String formatKey(final int key) {
        Utils.validate(readGroupReverseLookupTable.containsKey(key), () -> "missing key " + key);
        return readGroupReverseLookupTable.get(key);
    }

    @Override
    public int keyFromValue(final Object value) { // tsato: this use of key and value is problematic...
        return keyForReadGroup((String) value);
    }


    /**
     *
     * @param readGroupId
     * @return The index if the read group exists. -1 if the read group does not exist in the recal table,
     *         but the parameter is set. (REWORD)
     */
    private int keyForReadGroup(final String readGroupId) {
        if (readGroupLookupTable.containsKey(readGroupId)) { // tsato: should I also check "allowMissingReadGroup" here? Add it as an instance variable for this class...
            return readGroupLookupTable.get(readGroupId);
        } else { // tsato: perhaps we can return -1, and error
            return MISSING_READ_GROUP_KEY; // tsato: but we shouldn't return -1 when calling recalibration engine....ok with apply bqsr...caller should throw error as needed
        }
//        } else {
//            throw new GATKException("The covariates table is missing " + RecalUtils.READGROUP_COLUMN_NAME + " " + readGroupId + " in " + RecalUtils.READGROUP_REPORT_TABLE_TITLE);
//        }
    }

    @Override
    public int maximumKeyValue() {
        return readGroupLookupTable.size() - 1;
    }

    public static List<String> getReadGroupIDs(final SAMFileHeader header) {
        return header.getReadGroups().stream().map(rg -> getReadGroupIdentifier(rg)).collect(Collectors.toList());
    }
}
