package org.broadinstitute.hellbender.utils.read.markduplicates;

import htsjdk.samtools.ReservedTagConstants;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.Histogram;

import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Optional;

/**
 * A class to generate library Ids and keep duplication metrics by library IDs.
 *
 * @author nhomer
 */
public final class LibraryIdGenerator {

    public static final String UNKNOWN_LIBRARY = "Unknown Library";

    private final SAMFileHeader header;
    private final Map<String, Short> libraryIds = new LinkedHashMap<>(); // from library string to library id
    private short nextLibraryId = 1;
    private final Map<String, GATKDuplicationMetrics> metricsByLibrary = new LinkedHashMap<>();
    private final Histogram<Short> opticalDuplicatesByLibraryId = new Histogram<>();


    public LibraryIdGenerator(final SAMFileHeader header) {
        this.header = header;

        for (final SAMReadGroupRecord readGroup : header.getReadGroups()) {
            final String library = getReadGroupLibraryName(readGroup);
            GATKDuplicationMetrics metrics = metricsByLibrary.get(library);
            if (metrics == null) {
                metrics = new GATKDuplicationMetrics();
                metrics.LIBRARY = library;
                metricsByLibrary.put(library, metrics);
            }
        }
    }

    public Map<String, Short> getLibraryIdsMap() { return this.libraryIds; }

    public Map<String, GATKDuplicationMetrics> getMetricsByLibraryMap() { return this.metricsByLibrary; }

    public Histogram<Short> getOpticalDuplicatesByLibraryIdMap() { return this.opticalDuplicatesByLibraryId; }

	public static String getReadGroupLibraryName(final SAMReadGroupRecord readGroup) {
		return Optional.ofNullable(readGroup.getLibrary())
				.orElse(UNKNOWN_LIBRARY);
	}

    /**
     * Gets the library name from the header for the record. If the RG tag is not present on
     * the record, or the library isn't denoted on the read group, a constant string is
     * returned.
     */
    public static String getLibraryName(final SAMFileHeader header, final SAMRecord rec) {
        final String readGroupId = (String) rec.getAttribute(ReservedTagConstants.READ_GROUP_ID);
        return getLibraryName(header, readGroupId);
    }

    /**
     * Gets the library name from the header for the read group id. If the read group id is null
     * or the library isn't denoted on the read group, a constant string is
     * returned.
     */
    public static String getLibraryName(final SAMFileHeader header, String readGroupId) {
      if (readGroupId != null) {
            final SAMReadGroupRecord rg = header.getReadGroup(readGroupId);
            if (rg != null) {
                final String libraryName = rg.getLibrary();
                if (null != libraryName) return libraryName;
            }
        }

        return UNKNOWN_LIBRARY;
    }

    /** Get the library ID for the given SAM record. */
    public short getLibraryId(final SAMRecord rec) {
        final String library = getLibraryName(this.header, rec);
        Short libraryId = this.libraryIds.get(library);

        if (libraryId == null) {
            libraryId = this.nextLibraryId++;
            this.libraryIds.put(library, libraryId);
        }

        return libraryId;
    }

    public GATKDuplicationMetrics getMetricsByLibrary(final String library) {
        return this.metricsByLibrary.get(library);
    }

    public void addMetricsByLibrary(final String library, final GATKDuplicationMetrics metrics) {
        this.metricsByLibrary.put(library, metrics);
    }

    public long getNumberOfOpticalDuplicateClusters() {
        return (long) this.opticalDuplicatesByLibraryId.getSumOfValues();
    }
}
