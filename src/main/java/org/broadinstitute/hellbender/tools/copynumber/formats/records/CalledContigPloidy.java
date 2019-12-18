package org.broadinstitute.hellbender.tools.copynumber.formats.records;

import org.broadinstitute.hellbender.utils.Utils;

import java.util.Objects;

public class CalledContigPloidy {
    private String contig;
    private int ploidy;
    private double ploidyQuality;

    public CalledContigPloidy(final String contig, final int ploidy, final double ploidyQuality) {
        Utils.nonNull(contig, "Ploidy call contig cannot be null");
        Utils.validateArg(ploidy >= 0, "Negative ploidy is not allowed.");
        Utils.validateArg(ploidyQuality >= 0, "Negative ploidy quality is not allowed.");
        Utils.validateArg(Double.isFinite(ploidy), "Non-finite ploidy is not allowed.");
        Utils.validateArg(Double.isFinite(ploidyQuality), "Non-finite ploidy quality is not allowed.");
        this.contig = contig;
        this.ploidy = ploidy;
        this.ploidyQuality = ploidyQuality;
    }

    public String getContig() {
        return contig;
    }

    public int getPloidy() {
        return ploidy;
    }

    public double getPloidyQuality() {
        return ploidyQuality;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof CalledContigPloidy)) return false;
        CalledContigPloidy that = (CalledContigPloidy) o;
        return ploidy == that.ploidy &&
                Double.compare(that.ploidyQuality, ploidyQuality) == 0 &&
                contig.equals(that.contig);
    }

    @Override
    public int hashCode() {
        return Objects.hash(contig, ploidy, ploidyQuality);
    }

    @Override
    public String toString() {
        return "CalledContigPloidy{" +
                "contig='" + contig + '\'' +
                ", ploidy=" + ploidy +
                ", ploidyQuality=" + ploidyQuality +
                '}';
    }
}
