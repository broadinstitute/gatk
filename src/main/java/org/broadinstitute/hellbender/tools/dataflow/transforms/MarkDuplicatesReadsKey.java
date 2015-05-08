package org.broadinstitute.hellbender.tools.dataflow.transforms;

import com.google.api.services.genomics.model.Read;
import com.google.cloud.genomics.gatk.common.GenomicsConverter;
import htsjdk.samtools.ReservedTagConstants;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;

/**
 * Encodes a unique key for read, read pairs and fragments. Used to identify duplicates for MarkDuplicates.
 */
public final class MarkDuplicatesReadsKey {

    /**
     * Makes a unique key for the fragment.
     */
    public static String keyForFragment(final SAMFileHeader header, final Read read) {
        //HACK: convert back to SAMRecord but fastest to code
        final SAMRecord samRecord = GenomicsConverter.makeSAMRecord(read, header);
        return String.format(
                "%s|%d|%d|%s",
                library(header, samRecord),
                index(header, samRecord.getReferenceName()),
                unclippedCoordinate(samRecord),
                orientation(samRecord));
    }

    /**
     * Makes a unique key for the paired reads.
     */
    public static String keyForPairedEnds(final SAMFileHeader header, final Read first, final Read second) {
        final String key = keyForFragment(header, first);
        if (second == null) {
            return key;
        }
        //HACK: convert back to SAMRecord but fastest to code
        final SAMRecord samRecordSecond = GenomicsConverter.makeSAMRecord(second, header);
        return String.format(
                "%s|%d|%d|%s",
                key,
                index(header, samRecordSecond.getReferenceName()),
                unclippedCoordinate(samRecordSecond),
                orientation(samRecordSecond));
    }

    /**
     * Makes a unique key for the pair.
     */
    public static String keyForPair(final SAMFileHeader header, final Read read) {
        return String.format(
                "%s|%s",
                read.getReadGroupId(),
                read.getFragmentName());
    }


    private static String library(final SAMFileHeader header, final SAMRecord record) {
        String lib = null;
        final String rg = record.getStringAttribute(ReservedTagConstants.READ_GROUP_ID);
        if (rg != null) {
            final SAMReadGroupRecord rgr = header.getReadGroup(rg);
            if (rgr != null) {
                lib = rgr.getLibrary();
            }
        }
        return (lib == null) ? "-" : lib;
    }

    public static int index(final SAMFileHeader header, final String ref) {
        return header.getSequenceDictionary().getSequenceIndex(ref);
    }

    static int unclippedCoordinate(final SAMRecord record) {
        return record.getReadNegativeStrandFlag()
                ? record.getUnclippedEnd()
                : record.getUnclippedStart();
    }

    private static String orientation(final SAMRecord record) {
        return record.getReadNegativeStrandFlag() ? "r" : "f";
    }
}

