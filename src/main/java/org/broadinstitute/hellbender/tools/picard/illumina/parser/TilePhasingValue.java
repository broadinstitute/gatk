package org.broadinstitute.hellbender.tools.picard.illumina.parser;

/**
 * Captures information about a phasing value - Which read it corresponds to, which phasing type and a median value
 *
 *  @author jgentry
 */
public class TilePhasingValue {
    private final TileTemplateRead tileTemplateRead;
    private final float phasingValue;
    private final float prePhasingValue;

    public TilePhasingValue(final TileTemplateRead tileTemplateRead, final float phasingValue, final float prePhasingValue) {
        this.tileTemplateRead = tileTemplateRead;
        this.phasingValue = phasingValue;
        this.prePhasingValue = prePhasingValue;
    }

    public TileTemplateRead getTileTemplateRead() {
        return tileTemplateRead;
    }

    public float getPhasingValue() {
        return phasingValue;
    }

    public float getPrePhasingValue() {
        return prePhasingValue;
    }
}
