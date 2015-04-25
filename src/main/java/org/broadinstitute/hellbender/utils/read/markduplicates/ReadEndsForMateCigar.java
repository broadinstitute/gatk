package org.broadinstitute.hellbender.utils.read.markduplicates;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.util.SamRecordWithOrdinal;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.util.List;
import java.util.Set;

/**
 * A class to store individual records for MarkDuplicatesWithMateCigar.  This aids in comparing records to determine which need to
 * be compared when we mark duplicates.  We also store the original SAMRecord and its ordinal in the input file (in SamRecordWithOrdinal) to
 * access optional tags (mate cigar) and other information.
 */
public class ReadEndsForMateCigar extends ReadEnds {
    // to see if either end is unmapped
    byte hasUnmapped = 0;

    // we need this reference so we can access the mate cigar among other things
    public SamRecordWithOrdinal samRecordWithOrdinal = null;

    /**
     * Physical locations used for optical duplicate tracking.  This is only stored for paired end reads where both ends are mapped,
     * and when we see the first mate.
     */
    private PhysicalLocationForMateCigarSet locationSet = null;

    /** Builds a read ends object that represents a single read. */
    public ReadEndsForMateCigar(final SAMFileHeader header, final SamRecordWithOrdinal samRecordWithOrdinal,
                                final OpticalDuplicateFinder opticalDuplicateFinder, final short libraryId) {

        this.readGroup = -1;
        this.tile = -1;
        this.x = this.y = -1;
        this.read2ReferenceIndex = this.read2Coordinate = -1;
        this.hasUnmapped = 0;

        this.samRecordWithOrdinal = samRecordWithOrdinal;

        final SAMRecord record = this.samRecordWithOrdinal.getRecord();

        this.read1ReferenceIndex = record.getReferenceIndex();
        this.read1Coordinate = record.getReadNegativeStrandFlag() ? record.getUnclippedEnd() : record.getUnclippedStart();
        if (record.getReadUnmappedFlag()) {
            throw new UserException("Found an unexpected unmapped read");
        }

        if (record.getReadPairedFlag() && !record.getReadUnmappedFlag() && !record.getMateUnmappedFlag()) {
            this.read2ReferenceIndex = record.getMateReferenceIndex();
            this.read2Coordinate = record.getMateNegativeStrandFlag() ? SAMUtils.getMateUnclippedEnd(record) : SAMUtils.getMateUnclippedStart(record);

            // set orientation
            this.orientation = ReadEnds.getOrientationByte(record.getReadNegativeStrandFlag(), record.getMateNegativeStrandFlag());

            // Set orientationForOpticalDuplicates, which always goes by the first then the second end for the strands.  NB: must do this
            // before updating the orientation later.
            if (record.getReadPairedFlag()) {
                if (record.getFirstOfPairFlag()) {
                    this.orientationForOpticalDuplicates = ReadEnds.getOrientationByte(record.getReadNegativeStrandFlag(), record.getMateNegativeStrandFlag());
                } else {
                    this.orientationForOpticalDuplicates = ReadEnds.getOrientationByte(record.getMateNegativeStrandFlag(), record.getReadNegativeStrandFlag());
                }
            }
        } else {
            this.orientation = record.getReadNegativeStrandFlag() ? ReadEndsForMateCigar.R : ReadEndsForMateCigar.F;
        }

        // Fill in the library ID
        this.libraryId = libraryId;

        // Is this unmapped or its mate?
        if (record.getReadUnmappedFlag() || (record.getReadPairedFlag() && record.getMateUnmappedFlag())) {
            this.hasUnmapped = 1;
        }

        // Fill in the location information for optical duplicates
        if (opticalDuplicateFinder.addLocationInformation(record.getReadName(), this)) {
            // calculate the RG number (nth in list)
            // NB: could this be faster if we used a hash?
            this.readGroup = 0;
            final String rg = (String) record.getAttribute("RG");
            final List<SAMReadGroupRecord> readGroups = header.getReadGroups();
            if (rg != null && readGroups != null) {
                for (final SAMReadGroupRecord readGroup : readGroups) {
                    if (readGroup.getReadGroupId().equals(rg)) break;
                    else this.readGroup++;
                }
            }
        }
    }

    /** A number of convenience functions */
    public SamRecordWithOrdinal getSamRecordIndex() { return this.samRecordWithOrdinal; }

    public SAMRecord getRecord() { return this.samRecordWithOrdinal.getRecord(); }

    public String getRecordReadName() { return this.samRecordWithOrdinal.getRecord().getReadName(); }

    @Override
    public boolean isPaired() { return this.getRecord().getReadPairedFlag(); }

    /** Gets the read ends for optical duplicate tracking */
    public Set<ReadEnds> getReadEndSetForOpticalDuplicates() {
        if (null == this.locationSet) throw new GATKException("Already called getReadEndSetForOpticalDuplicates");
        final Set<ReadEnds> locationSet = this.locationSet.getReadEnds();
        this.locationSet = null;
        return locationSet;
    }

    public PhysicalLocationForMateCigarSet getLocationSet() {
        return this.locationSet;
    }

    public PhysicalLocationForMateCigarSet removeLocationSet() {
        final PhysicalLocationForMateCigarSet locationSet = this.locationSet;
        this.locationSet = null;
        return locationSet;
    }

    public void setLocationSet(final PhysicalLocationForMateCigarSet locationSet) {
        this.locationSet = locationSet;
    }

}
