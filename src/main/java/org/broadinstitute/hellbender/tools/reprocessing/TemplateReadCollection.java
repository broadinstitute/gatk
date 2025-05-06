package org.broadinstitute.hellbender.tools.reprocessing;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.PeekableIterator;
import org.broadinstitute.hellbender.utils.runtime.ProgressLogger;

import java.util.ArrayList;
import java.util.List;

// This class is used to store and read in a set of records with the same query name from the BAM file, i.e.
// reads that come from the same template
public class TemplateReadCollection {
    private final String queryName;
    private final List<SAMRecord> records;
    private final SAMRecord firstInPairRecord;
    private final SAMRecord secondInPairRecord;

    /**
     * Note that this constructor changes the state of the iterator passed to it. Should only be used within this class.
     *
     * @param iterator PeekableIterator of SAMRecords
     */
    public TemplateReadCollection(final PeekableIterator<SAMRecord> iterator, final ProgressLogger progressLogger) {
        if (!iterator.hasNext()) {
            // there must be more records in this iterator, throw exception
            throw new IllegalArgumentException("No more records in the iterator");
        }

        final SAMRecord firstRecord = iterator.next();
        this.queryName = firstRecord.getReadName();
        this.records = new ArrayList<>();
        this.records.add(firstRecord);
        progressLogger.record(firstRecord);
        while (iterator.hasNext()) {
            final SAMRecord record = iterator.peek();
            if (record.getReadName().equals(queryName)) {
                records.add(iterator.next());
                progressLogger.record(record);
            } else {
                break;
            }
        }

        SAMRecord localFirstInPairRecord = null;
        SAMRecord localSecondInPairRecord = null;

        for (final SAMRecord record : records) {
            if (record.isSecondaryOrSupplementary()) {
                continue;
            }

            if (record.getFirstOfPairFlag()) {
                if (localFirstInPairRecord != null) {
                    throw new IllegalStateException("Multiple first in pair records found for query name: " + queryName);
                }
                localFirstInPairRecord = record;
            } else if (record.getSecondOfPairFlag()) {
                if (localSecondInPairRecord != null) {
                    throw new IllegalStateException("Multiple second in pair records found for query name: " + queryName);
                }
                localSecondInPairRecord = record;
            }
        }

        if (localFirstInPairRecord == null || localSecondInPairRecord == null) {
            throw new IllegalStateException("Missing first or second in pair record for query name: " + queryName);
        }
        this.firstInPairRecord = localFirstInPairRecord;
        this.secondInPairRecord = localSecondInPairRecord;
    }

    // Getters
    public String getQueryName() {
        return queryName;
    }

    public List<SAMRecord> getRecords() {
        return records;
    }

    public SAMRecord getFirstInPairRecord() {
        return firstInPairRecord;
    }

    public SAMRecord getSecondInPairRecord() {
        return secondInPairRecord;
    }
}
