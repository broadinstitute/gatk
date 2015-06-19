package org.broadinstitute.hellbender.tools.dataflow;

import com.google.api.services.genomics.model.CigarUnit;
import com.google.api.services.genomics.model.Position;
import com.google.api.services.genomics.model.Read;
import htsjdk.samtools.*;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.StreamSupport;

import static java.lang.Boolean.TRUE;

import static java.lang.Boolean.TRUE;

/**
 * Utilities for reading with Reads (those from com.google.api.services.genomics.model).
 */
public final class GenomicsReadUtils {
    private GenomicsReadUtils(){}


    /**
     * Map form CIGAR operations as represented in the API to standard SAM ones.
     * HACK - copied from GenomicsConverter becuase it's private.
     */
    private static final Map<String, String> CIGAR_OPERATIONS = new HashMap<>();
    static {
        CIGAR_OPERATIONS.put("ALIGNMENT_MATCH", "M");
        CIGAR_OPERATIONS.put("CLIP_HARD", "H");
        CIGAR_OPERATIONS.put("CLIP_SOFT", "S");
        CIGAR_OPERATIONS.put("DELETE", "D");
        CIGAR_OPERATIONS.put("INSERT", "I");
        CIGAR_OPERATIONS.put("PAD", "P");
        CIGAR_OPERATIONS.put("SEQUENCE_MATCH", "=");
        CIGAR_OPERATIONS.put("SEQUENCE_MISMATCH", "X");
        CIGAR_OPERATIONS.put("SKIP", "N");
    }

    /**
     * Retruns true if the read is a primary alignment, false otherwise.
     */
    public static boolean isPrimaryAlignment(final Read read) {
        //Note: have to do TRUE dancing because of nulls.
        return read.getAlignment() != null && !TRUE.equals(read.getSecondaryAlignment()) && !TRUE.equals(read.getSupplementaryAlignment());
    }

    /**
     * @return the alignment start (1-based, inclusive) adjusted for clipped bases.  For example if the read
     * has an alignment start of 100 but the first 4 bases were clipped (hard or soft clipped)
     * then this method will return 96.
     *
     * Invalid to call on an unmapped read.
     */
    public static int getUnclippedStart(final Read read) {
        return SAMUtils.getUnclippedStart((int) getAlignmentStart(read), getCigar(read));
    }

    public static int unclippedCoordinate(final Read record) {
        return isNegativeStrand(record)
                ? getUnclippedEnd(record)
                : getUnclippedStart(record);
    }

    /**
     * @return the alignment end (1-based, inclusive) adjusted for clipped bases.  For example if the read
     * has an alignment end of 100 but the last 7 bases were clipped (hard or soft clipped)
     * then this method will return 107.
     *
     * Invalid to call on an unmapped read.
     */
    public static int getUnclippedEnd(final Read read) {
        return SAMUtils.getUnclippedEnd(getAlignmentEnd(read), getCigar(read));
    }

    public static Cigar getCigar(final Read read) {
        final List<CigarUnit> cigarUnits = read.getAlignment().getCigar();
        final StringBuilder buff = new StringBuilder(cigarUnits.size()*2);//guess
        for (final CigarUnit cu : cigarUnits){
            buff.append(cu.getOperationLength()).append(CIGAR_OPERATIONS.get(cu.getOperation()));
        }
        final String cigarString = buff.toString();
        return TextCigarCodec.decode(cigarString);
    }

    /**
     * @return 1-based inclusive rightmost position of the clipped sequence, or {@link SAMRecord.NO_ALIGNMENT_START} read if unmapped.
     */
    public static int getAlignmentEnd(final Read read) {
        if (isUnmapped(read)) {
            return SAMRecord.NO_ALIGNMENT_START;
        } else {
            return (int)getAlignmentStart(read) + getCigar(read).getReferenceLength() - 1;
        }
    }

    public static boolean isUnmapped(final Read read) {
        return read.getAlignment() == null || ! isMappedPosition(read.getAlignment().getPosition());
    }

    public static boolean isMappedPosition(final Position position) {
        return position != null && position.getPosition() != null;
    }

    public static boolean isPaired(final Read read) {
        return isMappedPosition(read.getNextMatePosition());
    }

    public static boolean isNegativeStrand(final Read record) {
        return record.getAlignment().getPosition().getReverseStrand();
    }

    public static String orientation(final Read record) {
        return isNegativeStrand(record) ? "r" : "f";
    }

    public static String library(final SAMFileHeader header, final Read record) {
        String lib = null;
        final String rg = record.getReadGroupId();
        if (rg != null) {
            final SAMReadGroupRecord rgr = header.getReadGroup(rg);
            if (rgr != null) {
                lib = rgr.getLibrary();
            }
        }
        return (lib == null) ? "-" : lib;
    }


    /**
     * Returns the mapping quality.
     */
    public static int getMappingQuality(final Read read){
        return read.getAlignment().getMappingQuality();
    }

    /**
     * Returns the name of the reference sequence for the given read (ie which contig the read is aligned to)
     * or null if that information is not available.
     */
    public static String getReferenceName(final Read read){
        if (read == null || read.getAlignment() == null || read.getAlignment().getPosition() == null ){
            return null;
        }
        return read.getAlignment().getPosition().getReferenceName();
    }

    public static long getAlignmentStart(final Read read) {
        return read.getAlignment().getPosition().getPosition();
    }

    public static String getMateReferenceName(final Read read) {
        if (read == null || read.getNextMatePosition() == null ){
            return null;
        }

        return read.getNextMatePosition().getReferenceName();
    }

    public static long getMateAlignmentStart(final Read read) {
        return read.getNextMatePosition().getPosition();
    }
}