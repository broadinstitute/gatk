package org.broadinstitute.hellbender.tools.exome.conversion.acsconversion;

import org.broadinstitute.hellbender.tools.exome.ModeledSegment;
import org.broadinstitute.hellbender.utils.SimpleInterval;


/**
 * Class to represent AllelicCapSeg (R version) segments.  This class should only be used for conversion operations.
 *
 * f - minor allelic fraction
 * tau - total copy ratio * 2
 * sigmaTau - confidence interval around tau.
 * muMinor - minor allelic copy number estimate: tau * f
 * muMajor - major allelic copy number estimate: tau * (1-f)
 * sigmaMinor and Major -- equal to each other.  Estimated error around muMinor and muMajor based on the credible
 *  interval of f from ACNV.
 * segLabelCNLOH -- always set to 2
 * hetCount -- number of hets in a segment
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