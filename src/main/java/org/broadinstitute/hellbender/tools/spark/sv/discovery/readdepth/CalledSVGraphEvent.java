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
    public static String PROBABILITY_COLUMN_STRING = "P";
    public static String HOMOZYGOUS_COLUMN_STRING = "HOM";
    
    private final int groupId;
    private final int pathId;
    private final Type type;
    private final SVInterval interval;
    private final boolean resolved; //False if solution not found
    private final double probability;
    private final boolean homozygous;

    public CalledSVGraphEvent(final Type type, final SVInterval interval,
                              final int groupId, final int pathId,
                              final boolean resolved, final double probability,
                              final boolean homozygous) {
        Utils.nonNull(type, "Type cannot be null");
        Utils.nonNull(interval, "Interval cannot be null");
        this.type = type;
        this.interval = interval;
        this.groupId = groupId;
        this.pathId = pathId;
        this.resolved = resolved;
        this.probability = probability;
        this.homozygous = homozygous;
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

    public double getProbability() {
        return probability;
    }

    public boolean isHomozygous() { return homozygous; }

    public static String bedHeader() {
        return String.join("\t", Arrays.asList(CONTIG_COLUMN_STRING, START_COLUMN_STRING, END_COLUMN_STRING, TYPE_COLUMN_STRING, PROBABILITY_COLUMN_STRING, HOMOZYGOUS_COLUMN_STRING, SIZE_COLUMN_STRING, CalledSVGraphGenotype.GRAPH_ID_COLUMN, CalledSVGraphGenotype.GENOTYPE_ID_COLUMN));
    }

    public String bedString(final SAMSequenceDictionary dictionary) {
        Utils.nonNull(dictionary, "Dictionary cannot be null");
        final SimpleInterval simpleInterval = SVIntervalUtils.convertToSimpleInterval(interval, dictionary);
        return simpleInterval.getContig() + "\t" + simpleInterval.getStart() + "\t" + simpleInterval.getEnd() + "\t" +
                getType().toString() + "\t" + getProbability() + "\t" + (homozygous ? "1" : "0") + "\t" +
                getInterval().getLength() + "\t" + getGroupId() + "\t" + getPathId();
    }
}
