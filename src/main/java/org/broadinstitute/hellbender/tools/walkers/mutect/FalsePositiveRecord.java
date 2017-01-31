package org.broadinstitute.hellbender.tools.walkers.mutect;

/**
 * Created by tsato on 2/14/17.
 */
public class FalsePositiveRecord {
    public static final String ID_COLUMN_NAME = "id";
    public static final String SNP_COLUMN_NAME = "snp";
    public static final String INDEL_COLUMN_NAME = "indel";
    public static final String SNP_FPR_COLUMN_NAME = "snp_FPR";
    public static final String INDEL_FPR_COLUMN_NAME = "indel_FPR";
    public static final String TARGET_TERRITORY_COLUMN_NAME = "target_territory";

    private String id;
    private long snpFalsePositives;
    private long indelFalsePositives;
    private long targetTerritory;

    public FalsePositiveRecord(final String id, final long snpFalsePositives, final long indelFalsePositives, final long targetTerritory) {
        this.id = id;
        this.snpFalsePositives = snpFalsePositives;
        this.indelFalsePositives = indelFalsePositives;
        this.targetTerritory = targetTerritory;
    }

    public String getId(){
        return id;
    }

    public long getSnpFalsePositives(){
        return snpFalsePositives;
    }

    public long getIndelFalsePositives(){
        return indelFalsePositives;
    }

    public long getTargetTerritory() {
        return targetTerritory;
    }

    // report FPR in units of per megabases
    public double getSnpFalsePositiveRate(){
        return snpFalsePositives / targetTerritory * 1000000.0;
    }

    public double getIndelFalsePositiveRate(){
        return indelFalsePositives / targetTerritory * 1000000.0;
    }

}
