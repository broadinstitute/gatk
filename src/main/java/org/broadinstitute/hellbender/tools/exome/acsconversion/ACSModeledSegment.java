package org.broadinstitute.hellbender.tools.exome.acsconversion;

import org.broadinstitute.hellbender.tools.exome.ModeledSegment;
import org.broadinstitute.hellbender.utils.SimpleInterval;


/**
 * Class to represent AllelicCapSeg (R version) segments.  This class should only be used for conversion operations.
 */
public class ACSModeledSegment extends ModeledSegment {

    final double f;
    final double tau;
    final double sigmaTau;
    final double muMinor;
    final double sigmaMinor;
    final double muMajor;
    final double sigmaMajor;
    final int segLabelCNLOH;
    final int hetCount;

    public ACSModeledSegment(final SimpleInterval interval, final String call, final long targetCount, final double segmentMeanInLog2CR, final int hetCount, final double f, final double sigmaTau, final double muMinor, final double sigmaMinor, final double muMajor, final double sigmaMajor, final int segLabelCNLOH) {
        super(interval, call, targetCount, segmentMeanInLog2CR);
        this.f = f;
        this.tau = getSegmentMeanInCRSpace() * 2;
        this.sigmaTau = sigmaTau;
        this.muMinor = muMinor;
        this.sigmaMinor = sigmaMinor;
        this.muMajor = muMajor;
        this.sigmaMajor = sigmaMajor;
        this.segLabelCNLOH = segLabelCNLOH;
        this.hetCount = hetCount;
    }

    public double getF() {
        return f;
    }

    public double getTau() {
        return tau;
    }

    public double getSigmaTau() {
        return sigmaTau;
    }

    public double getMuMinor() {
        return muMinor;
    }

    public double getSigmaMinor() {
        return sigmaMinor;
    }

    public double getMuMajor() {
        return muMajor;
    }

    public double getSigmaMajor() {
        return sigmaMajor;
    }

    public int getSegLabelCNLOH() {
        return segLabelCNLOH;
    }
    public int getHetCount() {
        return hetCount;
    }
}