package org.broadinstitute.hellbender.tools.spark.sv;

/**
 * Created by valentin on 6/1/17.
 */
public class AlignmentPenalties {

    public double unmappedMate = -4.0 * Math.log(10);
    public double unmappedPair = -8.0 * Math.log(10);
    public double softClip = -3.0 * Math.log(10) * 4;
    public double inversion = -4.0 * Math.log(10);
    public double indelStart = -4.5 * Math.log(10);
    public double indelExtend = Math.log(0.1);

}
