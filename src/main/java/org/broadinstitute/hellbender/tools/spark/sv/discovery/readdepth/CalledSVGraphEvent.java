package org.broadinstitute.hellbender.tools.spark.sv.discovery.readdepth;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Arrays;

/**
 * Structural variants called using an SVGraph
 */
public final class CalledSVGraphEvent {

    public enum Type {
        DEL,
        DUP,
        INV,
        DUP_INV,    //Duplicated inversion
        UR          //Unresolved
    }

    public static String CONTIG_COLUMN_STRING = "CHR";
    public static String START_COLUMN_STRING = "POS";
    public static String END_COLUMN_STRING = "END";
    public static String TYPE_COLUMN_STRING = "TYPE";
    public static String SIZE_COLUMN_STRING = "SIZE";
    public static String ID_COLUMN_STRING = "ID";
    public static String GENOTYPE_QUAL_COLUMN_STRING = "GQ";
    public static String REF_QUAL_COLUMN_STRING = "RQ";
    public static String HOMOZYGOUS_COLUMN_STRING = "HOM";
    public static String LOW_MAPPABILITY_COLUMN_STRING = "LOMAP";
    public static String DUPLICATE_COLUMN_STRING = "DUP";
    
    private final int groupId;
    private final int pathId;
    private final Type type;
    private final SVInterval interval;
    private final boolean resolved; //False if solution not found
    private double refQuality;
    private double genotypeQuality;
    private final boolean homozygous;
    private final SVGraphEdgeEvidence evidence;
    private boolean isLowMappability;
    private boolean isDuplicate;

    public CalledSVGraphEvent(final Type type, final SVInterval interval,
                              final int groupId, final int pathId,
                              final boolean resolved, final double refQuality,
                              final boolean homozygous, final SVGraphEdgeEvidence evidence) {
        Utils.nonNull(type, "Type cannot be null");
        Utils.nonNull(interval, "Interval cannot be null");
        this.type = type;
        this.interval = interval;
        this.groupId = groupId;
        this.pathId = pathId;
        this.resolved = resolved;
        this.refQuality = refQuality;
        this.homozygous = homozygous;
        this.genotypeQuality = -1;
        this.evidence = evidence;
        this.isLowMappability = false;
        this.isDuplicate = false;
    }

    public void setLowMappability(final boolean isLowMappability) {
        this.isLowMappability = isLowMappability;
    }

    public void setDuplicate(final boolean isDuplicate) {
        this.isDuplicate = isDuplicate;
    }

    public SVGraphEdgeEvidence getEvidence() {
        return evidence;
    }

    public boolean isResolved() {
        return resolved;
    }

    public int getGroupId() {
        return groupId;
    }

    public int getPathId() {
        return pathId;
    }

    public Type getType() {
        return type;
    }

    public SVInterval getInterval() {
        return interval;
    }

    public double getRefQuality() {
        return refQuality;
    }
    public void setRefQuality(final double refQuality) {
        this.refQuality = refQuality;
    }
    public double getGenotypeQuality() {
        return genotypeQuality;
    }
    public void setGenotypeQuality(final double genotypeQuality) {
        this.genotypeQuality = genotypeQuality;
    }

    public boolean isHomozygous() { return homozygous; }

    public static String bedHeader() {
        return String.join("\t", Arrays.asList(CONTIG_COLUMN_STRING, START_COLUMN_STRING, END_COLUMN_STRING,
                TYPE_COLUMN_STRING, ID_COLUMN_STRING, REF_QUAL_COLUMN_STRING, GENOTYPE_QUAL_COLUMN_STRING,
                HOMOZYGOUS_COLUMN_STRING, LOW_MAPPABILITY_COLUMN_STRING, DUPLICATE_COLUMN_STRING, SIZE_COLUMN_STRING,
                CalledSVGraphGenotype.GRAPH_ID_COLUMN, CalledSVGraphGenotype.GENOTYPE_ID_COLUMN));
    }

    public String bedString(final SAMSequenceDictionary dictionary) {
        Utils.nonNull(dictionary, "Dictionary cannot be null");
        final SimpleInterval simpleInterval = SVIntervalUtils.convertToSimpleInterval(interval, dictionary);
        return simpleInterval.getContig() + "\t" + simpleInterval.getStart() + "\t" + simpleInterval.getEnd() + "\t" +
                getType().toString() + "\t" + evidence.getId() + "\t" + getRefQuality() + "\t" + getGenotypeQuality() + "\t" + (homozygous ? "1" : "0") + "\t" +
                (isLowMappability ? "1" : "0") + "\t" + (isDuplicate ? "1" : "0") + "\t" + getInterval().getLength() + "\t" + getGroupId() + "\t" + getPathId();
    }
}
