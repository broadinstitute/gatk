package org.broadinstitute.hellbender.utils.read.markduplicates.sparksql;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.tools.spark.transforms.markduplicates.MarkDuplicatesSparkUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.read.markduplicates.LibraryIdGenerator;
import org.broadinstitute.hellbender.utils.read.markduplicates.MarkDuplicatesScoringStrategy;
import org.broadinstitute.hellbender.utils.read.markduplicates.ReadsKey;
import org.broadinstitute.hellbender.utils.read.markduplicates.sparkrecords.*;
import picard.sam.markduplicates.util.ReadEnds;
import picard.sam.util.PhysicalLocation;

import java.util.Map;

/** Java bean containing the raw fields of MarkDuplicatesSparkRecord for use with Spark Dataset Encoders. */
public class MarkDuplicatesSparkSqlRecord implements PhysicalLocation {
    // key fields
    private Long firstReadKeyValue;
    private Long secondReadKeyValue;

    // value fields
    private Integer partitionIndex;
    private String name;
    private Integer type;
    private Boolean isRead1ReverseStrand;
    private Boolean isRead2ReverseStrand;
    private Short score;
    private Boolean wasFlipped;

    // PhysicalLocation fields
    private short readGroup;
    private short tile;
    private int x;
    private int y;
    private short libraryId;

    public MarkDuplicatesSparkSqlRecord() {
    }

    public MarkDuplicatesSparkSqlRecord(Long firstReadKeyValue, Long secondReadKeyValue, Integer partitionIndex, String name, Integer type, Boolean isRead1ReverseStrand, Boolean isRead2ReverseStrand, Short score, Boolean wasFlipped) {
        this.firstReadKeyValue = firstReadKeyValue;
        this.secondReadKeyValue = secondReadKeyValue;
        this.partitionIndex = partitionIndex;
        this.name = name;
        this.type = type;
        this.isRead1ReverseStrand = isRead1ReverseStrand;
        this.isRead2ReverseStrand = isRead2ReverseStrand;
        this.score = score;
        this.wasFlipped = wasFlipped;
        this.readGroup = -1;
        this.tile = -1;
        this.x = -1;
        this.y = -1;
        this.libraryId = -1;
    }

    public MarkDuplicatesSparkSqlRecord(Long firstReadKeyValue, Long secondReadKeyValue, Integer partitionIndex, String name, Integer type, Boolean isRead1ReverseStrand, Boolean isRead2ReverseStrand, Short score, Boolean wasFlipped, short readGroup, short tile, int x, int y, short libraryId) {
        this.firstReadKeyValue = firstReadKeyValue;
        this.secondReadKeyValue = secondReadKeyValue;
        this.partitionIndex = partitionIndex;
        this.name = name;
        this.type = type;
        this.isRead1ReverseStrand = isRead1ReverseStrand;
        this.isRead2ReverseStrand = isRead2ReverseStrand;
        this.score = score;
        this.wasFlipped = wasFlipped;
        this.readGroup = readGroup;
        this.tile = tile;
        this.x = x;
        this.y = y;
        this.libraryId = libraryId;
    }

    public static MarkDuplicatesSparkSqlRecord newEmptyFragment(GATKRead read, SAMFileHeader header, Map<String, Byte> headerLibraryMap) {
        ReadsKey readsKey = ReadsKey.getKeyForFragment(ReadUtils.getStrandedUnclippedStart(read),
                read.isReverseStrand(),
                ReadUtils.getReferenceIndex(read, header),
                headerLibraryMap.get(MarkDuplicatesSparkUtils.getLibraryForRead(read, header, LibraryIdGenerator.UNKNOWN_LIBRARY)));
        ReadsKeySparkSqlRecord readsKeySparkSqlRecord = ReadsKeySparkSqlRecord.fromReadsKey(readsKey);
        return new MarkDuplicatesSparkSqlRecord(
                readsKeySparkSqlRecord.getFirstReadKeyValue(),
                readsKeySparkSqlRecord.getSecondReadKeyValue(),
                0,
                null,
                MarkDuplicatesSparkRecord.Type.EMPTY_FRAGMENT.ordinal(),
                read.isReverseStrand(),
                null,
                null,
                null);
    }

    public static MarkDuplicatesSparkSqlRecord newFragment(final GATKRead first, final SAMFileHeader header, int partitionIndex, MarkDuplicatesScoringStrategy scoringStrategy, Map<String, Byte> headerLibraryMap) {
        ReadsKey readsKey = ReadsKey.getKeyForFragment(ReadUtils.getStrandedUnclippedStart(first),
                first.isReverseStrand(),
                (short)ReadUtils.getReferenceIndex(first, header),
                headerLibraryMap.get(MarkDuplicatesSparkUtils.getLibraryForRead(first, header, LibraryIdGenerator.UNKNOWN_LIBRARY)));
        ReadsKeySparkSqlRecord readsKeySparkSqlRecord = ReadsKeySparkSqlRecord.fromReadsKey(readsKey);
        return new MarkDuplicatesSparkSqlRecord(
                readsKeySparkSqlRecord.getFirstReadKeyValue(),
                readsKeySparkSqlRecord.getSecondReadKeyValue(),
                partitionIndex,
                first.getName(),
                MarkDuplicatesSparkRecord.Type.FRAGMENT.ordinal(),
                first.isReverseStrand(),
                null,
                scoringStrategy.score(first),
                null);
    }

    public static MarkDuplicatesSparkSqlRecord newPair(GATKRead read1, GATKRead read2, SAMFileHeader header, int partitionIndex, MarkDuplicatesScoringStrategy scoringStrategy, Map<String, Byte> headerLibraryMap) {
        final String name1 = read1.getName();
        final String name2 = read2.getName();
        Utils.validate(name1.equals(name2), () -> "Paired reads have different names\n" + name1 + "\n" + name2);

        GATKRead first = read1;
        GATKRead second;

        final int read1UnclippedStart = ReadUtils.getStrandedUnclippedStart(read1);
        final int read2UnclippedStart = ReadUtils.getStrandedUnclippedStart(read2);

        int read1ReferenceIndex = ReadUtils.getReferenceIndex(read1,header);
        int read2ReferenceIndex = ReadUtils.getReferenceIndex(read2,header);

        if( read1ReferenceIndex != read2ReferenceIndex ? read1ReferenceIndex < read2ReferenceIndex : read1UnclippedStart <= read2UnclippedStart ){
            first = read1;
            second = read2;
        } else {
            first = read2;
            second = read1;
        }

        // if the two read ends are in the same position, pointing in opposite directions,
        // the orientation is undefined and the procedure above
        // will depend on the order of the reads in the file.
        // To avoid this, and match picard's behavior in this case, we ensure the orientation will be FR:
        if (read1ReferenceIndex == read2ReferenceIndex &&
                read1UnclippedStart == read2UnclippedStart &&
                first.isReverseStrand() && !second.isReverseStrand()) {
            // setting to FR for consistencies sake. (which involves flipping) if both reads had the same unclipped start
            GATKRead tmp = first;
            first = second;
            second = tmp;
        }

        ReadsKey readsKey = ReadsKey.getKeyForPair(header, first, second, headerLibraryMap);
        ReadsKeySparkSqlRecord readsKeySparkSqlRecord = ReadsKeySparkSqlRecord.fromReadsKey(readsKey);
        return new MarkDuplicatesSparkSqlRecord(
                readsKeySparkSqlRecord.getFirstReadKeyValue(),
                readsKeySparkSqlRecord.getSecondReadKeyValue(),
                partitionIndex,
                first.getName(),
                MarkDuplicatesSparkRecord.Type.PAIR.ordinal(),
                first.isReverseStrand(),
                second.isReverseStrand(),
                (short)(scoringStrategy.score(read1) + scoringStrategy.score(read2)),
                second.isFirstOfPair());
    }

    public static MarkDuplicatesSparkSqlRecord getPassthrough(GATKRead read, int partitionIndex) {
        ReadsKey readsKey = ReadsKey.hashKeyForPassthroughRead(read);
        ReadsKeySparkSqlRecord readsKeySparkSqlRecord = ReadsKeySparkSqlRecord.fromReadsKey(readsKey);
        return new MarkDuplicatesSparkSqlRecord(
                readsKeySparkSqlRecord.getFirstReadKeyValue(),
                readsKeySparkSqlRecord.getSecondReadKeyValue(),
                partitionIndex,
                read.getName(),
                MarkDuplicatesSparkRecord.Type.PASSTHROUGH.ordinal(),
                null,
                null,
                null,
                null);
    }

    public ReadsKeySparkSqlRecord key() {
        return new ReadsKeySparkSqlRecord(firstReadKeyValue, secondReadKeyValue);
    }

    public Long getFirstReadKeyValue() {
        return firstReadKeyValue;
    }

    public void setFirstReadKeyValue(Long firstReadKeyValue) {
        this.firstReadKeyValue = firstReadKeyValue;
    }

    public Long getSecondReadKeyValue() {
        return secondReadKeyValue;
    }

    public void setSecondReadKeyValue(Long secondReadKeyValue) {
        this.secondReadKeyValue = secondReadKeyValue;
    }

    public Integer getPartitionIndex() {
        return partitionIndex;
    }

    public void setPartitionIndex(Integer partitionIndex) {
        this.partitionIndex = partitionIndex;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public Integer getType() {
        return type;
    }

    public void setType(Integer type) {
        this.type = type;
    }

    public Boolean getRead1ReverseStrand() {
        return isRead1ReverseStrand;
    }

    public void setRead1ReverseStrand(Boolean read1ReverseStrand) {
        isRead1ReverseStrand = read1ReverseStrand;
    }

    public Boolean getRead2ReverseStrand() {
        return isRead2ReverseStrand;
    }

    public void setRead2ReverseStrand(Boolean read2ReverseStrand) {
        isRead2ReverseStrand = read2ReverseStrand;
    }

    public Short getScore() {
        return score;
    }

    public void setScore(Short score) {
        this.score = score;
    }

    public Boolean getWasFlipped() {
        return wasFlipped;
    }

    public void setWasFlipped(Boolean wasFlipped) {
        this.wasFlipped = wasFlipped;
    }

    // PhysicalLocation

    public short getReadGroup() {
        return readGroup;
    }

    public void setReadGroup(short readGroup) {
        this.readGroup = readGroup;
    }

    public short getTile() {
        return tile;
    }

    public void setTile(short tile) {
        this.tile = tile;
    }

    public int getX() {
        return x;
    }

    public void setX(int x) {
        this.x = x;
    }

    public int getY() {
        return y;
    }

    public void setY(int y) {
        this.y = y;
    }

    public short getLibraryId() {
        return libraryId;
    }

    public void setLibraryId(short libraryId) {
        this.libraryId = libraryId;
    }

    // For pair comparison

    /**
     * Returns the pair orientation suitable for optical duplicates,
     * which always goes by the first then the second end for the strands.
     * This is based on code in MarkDuplicatesGATK and ReadEnds.getOrientationByte.
     * Returns one of {@link ReadEnds#RR}, {@link ReadEnds#RF}, {@link ReadEnds#FR}, {@link ReadEnds#FF}
     */
    public byte getOrientationForOpticalDuplicates() {
        return getOrientation(true);
    }

    private byte getOrientation(final boolean optical) {
        if (isRead1ReverseStrand && isRead2ReverseStrand) {
            return ReadEnds.RR;
        }
        if (isRead1ReverseStrand) {
            return (optical && wasFlipped)? ReadEnds.FR : ReadEnds.RF; //at this point we know for sure isRead2ReverseStrand is false
        }
        if (isRead2ReverseStrand) {
            return (optical && wasFlipped)? ReadEnds.RF :ReadEnds.FR; //at this point we know for sure isRead1ReverseStrand is false
        }
        return ReadEnds.FF;  //at this point we know for sure isRead1ReverseStrand is false and isRead2ReverseStrand is false
    }

}
