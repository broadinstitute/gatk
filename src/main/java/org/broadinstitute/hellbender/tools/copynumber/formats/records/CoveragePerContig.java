package org.broadinstitute.hellbender.tools.copynumber.formats.records;

import org.broadinstitute.hellbender.tools.copynumber.DetermineGermlineContigPloidy;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.LinkedHashMap;

/**
 * Represents total coverage over each contig in an ordered set associated with a named sample.
 * Should only be used to write temporary files in {@link DetermineGermlineContigPloidy}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class CoveragePerContig {
    private final String sampleName;

    private final LinkedHashMap<String, Integer> coveragePerContig;

    public CoveragePerContig(final String sampleName,
                             final LinkedHashMap<String, Integer> coveragePerContig) {
        this.sampleName = Utils.nonEmpty(sampleName);
        this.coveragePerContig = Utils.nonNull(coveragePerContig);
    }

    public String getSampleName() {
        return sampleName;
    }

    public int getCoverage(final String contig) {
        return coveragePerContig.get(contig);
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (o == null || getClass() != o.getClass()) {
            return false;
        }

        final CoveragePerContig that = (CoveragePerContig) o;
        return sampleName.equals(that.sampleName) &&
                coveragePerContig.equals(that.coveragePerContig);
    }

    @Override
    public int hashCode() {
        int result = sampleName.hashCode();
        result = 31 * result + coveragePerContig.hashCode();
        return result;
    }
}
