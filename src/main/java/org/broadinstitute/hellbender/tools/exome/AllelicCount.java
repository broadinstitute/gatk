package org.broadinstitute.hellbender.tools.exome;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.tools.exome.allelefraction.MinorAlleleFractionCache;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

/**
 * Reference and alternate allele counts at a SNP site specified by an interval.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class AllelicCount implements Locatable {

    /* these are mandatory */
    private final SimpleInterval interval;
    private final int refReadCount, altReadCount;

    /* these are extra metadata and can be null */
    private final Nucleotide refNucleotide, altNucleotide;
    private final Integer readDepth;
    private final Double hetLogOdds;

    /**
     * construct the allelic count object
     * @param interval the genomic interval (it is assumed to have the same begin and end points)
     * @param refReadCount ref nucleotide count in the pileup
     * @param altReadCount alt nucleotide count in the pileup
     * @param refNucleotide ref nucleotide
     * @param altNucleotide alt nucleotide
     * @param hetLogOdds the log odds of the pileup corresponding to a heterozygous site
     */
    public AllelicCount(final SimpleInterval interval, final int refReadCount, final int altReadCount,
                        final Nucleotide refNucleotide, final Nucleotide altNucleotide,
                        final Integer readDepth, final Double hetLogOdds) {

        this.refReadCount = ParamUtils.isPositiveOrZero(refReadCount, "Can't construct AllelicCount with negative read counts.");
        this.altReadCount = ParamUtils.isPositiveOrZero(altReadCount, "Can't construct AllelicCount with negative read counts.");
        ParamUtils.isPositive(altReadCount + refReadCount, "Can't construct AllelicCount with zero total counts.");

        this.interval = Utils.nonNull(interval, "Can't construct AllelicCount with null interval.");

        /* these can be null */
        this.refNucleotide = refNucleotide;
        this.altNucleotide = altNucleotide;
        this.readDepth = readDepth;
        this.hetLogOdds = hetLogOdds;
    }

    public AllelicCount(final SimpleInterval interval, final int refReadCount, final int altReadCount) {
        this(interval, refReadCount, altReadCount, null, null, null, null);
    }

    public AllelicCount(final SimpleInterval interval, final int refReadCount, final int altReadCount,
                        final Nucleotide refNucleotide, final Nucleotide altNucleotide,
                        final Integer readDepth) {
        this(interval, refReadCount, altReadCount, refNucleotide, altNucleotide, readDepth, null);
    }

    @Override
    public String getContig() { return interval.getContig(); }

    @Override
    public int getStart() { return interval.getStart(); }

    @Override
    public int getEnd() { return interval.getEnd(); }

    public SimpleInterval getInterval() { return interval; }

    public int getRefReadCount() { return refReadCount; }

    public int getAltReadCount() { return altReadCount; }

    public Nucleotide getRefNucleotide() {
        return Utils.nonNull(refNucleotide, "The Ref nucleotide information is not available.");
    }

    public Nucleotide getAltNucleotide() {
        return Utils.nonNull(altNucleotide, "The Alt nucleotide information is not available.");
    }

    public Integer getReadDepth() {
        return Utils.nonNull(readDepth, "The read depth information is not available.");
    }

    public Double getHetLogOdds() {
        return Utils.nonNull(hetLogOdds, "The log odds of heterozygosity is not available.");
    }

    /**
     * Returns the maximum likelihood estimate of the alternate-allele fraction.
     * @return      alternate-allele fraction
     */
    public double estimateAltAlleleFraction() {
        return (double) altReadCount / (refReadCount + altReadCount);
    }

    /**
     * Returns the maximum-likelihood estimate of the minor-allele fraction given a specified allelic bias.
     * This likelihood is unimodal on [0,1/2] so numerical max-finding is straightforward.  See docs/CNVs/CNV-methods.pdf.
     *
     * @param allelicBias   allelic bias to use in estimate of minor allele fraction
     * @return      maximum likelihood estimate of the minor allele fraction
     */
    public double estimateMinorAlleleFraction(final double allelicBias) {
        ParamUtils.isPositiveOrZero(allelicBias, "Allelic bias must be non-negative.");
        return MinorAlleleFractionCache.get(altReadCount, refReadCount, allelicBias);
    }

    /**
     * Returns the maximum-likelihood estimate of the minor-allele fraction, assuming no allelic bias.
     *
     * With no allelic bias, the likelihood of a alt reads and r ref reads is proportional to f^a(1-f)^r + f^r(1-f)^a,
     * where f in [0,1/2] is the minor allele fraction and NOT the alt allele fraction.  The two terms derive from the
     * a priori equally likely cases that the alt and ref alleles are the minor allele, respectively.
     * This likelihood is unimodal on [0,1/2] so numerical max-finding is straightforward.
     *
     * @return      maximum likelihood estimate of the minor allele fraction
     */
    public double estimateMinorAlleleFraction() {
        return MinorAlleleFractionCache.get(altReadCount, refReadCount, 1.);
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (!(o instanceof AllelicCount)) {
            return false;
        }

        final AllelicCount count = (AllelicCount) o;
        return interval.equals(count.interval)
                && refReadCount == count.refReadCount && altReadCount == count.altReadCount
                && (!((refNucleotide != null) && (count.refNucleotide != null)) || (refNucleotide == count.refNucleotide))
                && (!((altNucleotide != null) && (count.altNucleotide != null)) || (altNucleotide == count.altNucleotide));
    }

    @Override
    public int hashCode() {
        int result = interval.hashCode();
        result = 31 * result + refReadCount;
        result = 31 * result + altReadCount;
        return result;
    }
}
