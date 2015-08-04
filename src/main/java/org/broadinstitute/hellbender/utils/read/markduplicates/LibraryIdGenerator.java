package org.broadinstitute.hellbender.utils.read.markduplicates;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.Histogram;

import java.util.HashMap;
import java.util.Map;

/**
 * A class to generate library Ids and keep duplication metrics by library IDs.
 *
 * @author nhomer
 */
public final class LibraryIdGenerator {

    private final SAMFileHeader header;
    private final Map<String, Short> libraryIds = new HashMap<>(); // from library string to library id
    private short nextLibraryId = 1;
    private final Map<String, DuplicationMetrics> metricsByLibrary = new HashMap<>();
    private final Histogram<Short> opticalDuplicatesByLibraryId = new Histogram<>();


    public LibraryIdGenerator(final SAMFileHeader header) {
        this.header = header;

        for (final SAMReadGroupRecord readGroup : header.getReadGroups()) {
            final String library = readGroup.getLibrary();
            DuplicationMetrics metrics = metricsByLibrary.get(library);
            if (metrics == null) {
                metrics = new DuplicationMetrics();
                metrics.LIBRARY = library;
                metricsByLibrary.put(library, metrics);
            }
        }
    }

    public Map<String, Short> getLibraryIdsMap() { return this.libraryIds; }

    public Map<String, DuplicationMetrics> getMetricsByLibraryMap() { return this.metricsByLibrary; }

    public Histogram<Short> getOpticalDuplicatesByLibraryIdMap() { return this.opticalDuplicatesByLibraryId; }

    /**
     * Gets the library name from the header for the record. If the RG tag is not present on
     * the record, or the library isn't denoted on the read group, a constant string is
     * returned.
     */
    public static String getLibraryName(final SAMFileHeader header, final SAMRecord rec) {
        final String readGroupId = (String) rec.getAttribute("RG");
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

        return "Unknown Library";
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

    public DuplicationMetrics getMetricsByLibrary(final String library) {
        return this.metricsByLibrary.get(library);
    }

    public void addMetricsByLibrary(final String library, final DuplicationMetrics metrics) {
        this.metricsByLibrary.put(library, metrics);
    }

    public long getNumberOfOpticalDuplicateClusters() {
        return (long) this.opticalDuplicatesByLibraryId.getSumOfValues();
    }
}
