package org.broadinstitute.hellbender.utils.read.markduplicates.sparkrecords;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.cmdline.argumentcollections.MarkDuplicatesSparkArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.transforms.markduplicates.MarkDuplicatesSparkUtils;
import org.broadinstitute.hellbender.utils.read.FlowBasedReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.read.markduplicates.LibraryIdGenerator;
import org.broadinstitute.hellbender.utils.read.markduplicates.MarkDuplicatesScoringStrategy;
import org.broadinstitute.hellbender.utils.read.markduplicates.ReadsKey;
import picard.sam.markduplicates.util.ReadEnds;

import java.util.Map;

/**
 * Class representing a single read fragment at a particular start location without a mapped mate.
 *
 * This class holds onto as little information as possible in an attempt to prevent excessive serialization of
 * during the processing step of MarkDuplicatesSpark
 */
public class FlowModeFragment extends TransientFieldPhysicalLocation {
    private static final long serialVersionUID = 0L;
    public static final String FLOW_DUPLICATE_SCORE_ATTR_NAME = "FlowDuplicateScore";

    protected transient ReadsKey key;
    protected int end = FlowBasedReadUtils.FLOW_BASED_INSIGNIFICANT_END;

    private final boolean R1R;

    protected final short score;

    public FlowModeFragment(final GATKRead first, final SAMFileHeader header, int partitionIndex, MarkDuplicatesScoringStrategy scoringStrategy, Map<String, Byte> headerLibraryMap, final MarkDuplicatesSparkArgumentCollection mdArgs) {
        super(partitionIndex, first.getName());

        this.R1R = first.isReverseStrand();

        int start = FlowBasedReadUtils.getStrandedUnclippedStartForFlow(first, header, mdArgs);
        if ( mdArgs.FLOW_END_LOCATION_SIGNIFICANT ) {
            this.end = FlowBasedReadUtils.getStrandedUnclippedEndForFlow(first, header, mdArgs);
        }
        this.key = ReadsKey.getKeyForFragment(start,
                isRead1ReverseStrand(),
                (short)ReadUtils.getReferenceIndex(first, header),
                headerLibraryMap.get(MarkDuplicatesSparkUtils.getLibraryForRead(first, header, LibraryIdGenerator.UNKNOWN_LIBRARY)));

        this.score = (this.end != FlowBasedReadUtils.FLOW_BASED_INSIGNIFICANT_END)
                ? ((mdArgs.FLOW_QUALITY_SUM_STRATEGY && FlowBasedReadUtils.isFlow(first)) ? computeFlowDuplicateScore(first, start, end) : scoringStrategy.score(first))
                : -1;
    }

    // compute fragment score using a flow-based specific method - cache in transient attribute
    private short computeFlowDuplicateScore(GATKRead rec, int start, int end) {

        Short storedScore = (Short)rec.getTransientAttribute(FLOW_DUPLICATE_SCORE_ATTR_NAME);
        if ( storedScore == null ) {
            short score = 0;

            score += (short) Math.min(getFlowSumOfBaseQualities(rec, start, end), Short.MAX_VALUE / 2);

            storedScore = score;
            rec.setTransientAttribute(FLOW_DUPLICATE_SCORE_ATTR_NAME, storedScore);
        }

        return storedScore;
    }

    /**
     * A quality summing scoring strategy used for flow based reads.
     *
     * The method walks on the bases of the read, in the synthesis direction. For each base, the effective
     * quality value is defined as the value on the first base on the hmer to which the base belongs to. The score
     * is defined to be the sum of all effective values above a given threshold.
     */
    private int getFlowSumOfBaseQualities(GATKRead rec, int start, int end) {
        int score = 0;

        if ( rec.isReverseStrand() ) {
            int     tmp = start;
            start = end;
            end = tmp;
        }

        // access qualities and bases
        byte[]      quals = rec.getBaseQualitiesNoCopy();
        byte[]      bases = rec.getBasesNoCopy();

        // establish range of bases/quals to work on
        int         startingOffset = Math.max(0, start - rec.getUnclippedStart());
        int         endOffset = Math.max(0, rec.getUnclippedEnd() - end);

        // loop on bases, extract qual related to homopolymer from start of homopolymer
        byte        lastBase = 0;
        byte        effectiveQual = 0;
        for ( int i = startingOffset ; i < bases.length - endOffset ; i++ ) {
            byte        base = bases[i];
            if ( base != lastBase ) {
                effectiveQual = quals[i];
            }
            if ( effectiveQual >= 15 ) {
                score += effectiveQual;
            }
            lastBase = base;
        }

        return score;
    }

    @Override
    public Type getType() {
      return Type.FRAGMENT;
    }

    @Override
    // NOTE: This is transient and thus may not exist if the object gets serialized
    public ReadsKey key() {
        return key;
    }

    @Override
    public short getScore() {
      return score;
    }

    @Override
    public boolean isRead1ReverseStrand() {
      return R1R;
    }

    @Override
    public byte getOrientationForPCRDuplicates() {
        return (R1R)? ReadEnds.R : ReadEnds.F;
    }

    @Override
    public String toString() {
        return "flow mode fragment: " + name;
    }

    public int getEnd() {
        return end;
    }
}
