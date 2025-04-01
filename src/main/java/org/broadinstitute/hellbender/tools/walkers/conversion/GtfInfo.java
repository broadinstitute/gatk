package org.broadinstitute.hellbender.tools.walkers.conversion;

import htsjdk.samtools.util.Interval;

/**
 * A class that represents information extracted from a feature in a GTF file.
 * The {@code GtfInfo} object encapsulates details about a specific gene or transcript,
 * including its type, gene name, and interval.
 *
 * <p>The {@code GtfInfo.Type} enum specifies whether the features is a
 * {@code GENE} or a {@code TRANSCRIPT}.</p>
 *
 * <p>The interval specifies the feature's contig (chromosome), start position,
 * and end position, providing the precise location of the gene or transcript
 * on the genome.</p>
 *
 * <p>Example usage:</p>
 * <pre>
 *     Interval interval = new Interval("chr1", 1000, 2000);
 *     GtfInfo gtfInfo = new GtfInfo(interval, GtfInfo.Type.GENE, "MAPK1");
 * </pre>
 **/

public class GtfInfo {

    public enum Type {
        GENE,
        TRANSCRIPT
    }

    private final Type type;
    private final String geneName;
    private final Interval interval;

    public GtfInfo(Interval interval, Type type, String geneName) {
        this.interval = interval;
        this.type = type;
        this.geneName = geneName;
    }

    public Type getType() {
        return type;
    }

    public String getGeneName() {
        return geneName;
    }

    public Interval getInterval() {
        return interval;
    }

    public Integer getStart() {
        return interval.getStart();
    }

    public Integer getEnd() {
        return interval.getEnd();
    }

    @Override
    public String toString() {
        return "GtfInfo{ " + "type = " + type + " geneName = " + geneName + "interval = " + interval;
    }

}
