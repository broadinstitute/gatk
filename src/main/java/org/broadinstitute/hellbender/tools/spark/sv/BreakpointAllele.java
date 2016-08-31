package org.broadinstitute.hellbender.tools.spark.sv;

import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.List;
import java.util.Objects;

import static org.broadinstitute.hellbender.tools.spark.sv.BreakpointAllele.InversionType.INV_3_TO_5;
import static org.broadinstitute.hellbender.tools.spark.sv.BreakpointAllele.InversionType.INV_5_TO_3;

/**
 * This class represents the allele of an SV breakpoint (a novel adjacency between two genomic locations
 */
class BreakpointAllele {
    final SimpleInterval leftAlignedLeftBreakpoint;
    final SimpleInterval leftAlignedRightBreakpoint;
    final String insertedSequence;
    final String homology;
    final InversionType inversionType;

    // not included in equals and hashCode so as not to break grouping by breakpoint allele if the mappings are different
    final List<String> insertionMappings;

    public BreakpointAllele(final SimpleInterval leftAlignedLeftBreakpoint, final SimpleInterval leftAlignedRightBreakpoint, final String insertedSequence, final String homology, final boolean fiveToThree, final boolean threeToFive, final List<String> insertionMappings) {
        this.leftAlignedLeftBreakpoint = leftAlignedLeftBreakpoint;
        this.leftAlignedRightBreakpoint = leftAlignedRightBreakpoint;
        this.insertedSequence = insertedSequence;
        this.homology = homology;
        this.inversionType = getInversionType(fiveToThree, threeToFive);
        this.insertionMappings = insertionMappings;
    }

    public boolean isInversion() {
        return leftAlignedLeftBreakpoint.getContig().equals(leftAlignedRightBreakpoint.getContig()) && (getInversionType() == INV_3_TO_5 || getInversionType() == INV_5_TO_3);
    }

    public enum InversionType{
        INV_3_TO_5, INV_5_TO_3, INV_NONE
    }

    public InversionType getInversionType() {
        return inversionType;
    }

    private InversionType getInversionType(final boolean fiveToThree, final boolean threeToFive){
        if(!fiveToThree && threeToFive){
            return InversionType.INV_3_TO_5;
        }else if(fiveToThree && !threeToFive){
            return InversionType.INV_5_TO_3;
        }
        return InversionType.INV_NONE;
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        final BreakpointAllele that = (BreakpointAllele) o;
        return inversionType == that.inversionType &&
                Objects.equals(leftAlignedLeftBreakpoint, that.leftAlignedLeftBreakpoint) &&
                Objects.equals(leftAlignedRightBreakpoint, that.leftAlignedRightBreakpoint) &&
                Objects.equals(insertedSequence, that.insertedSequence) &&
                Objects.equals(homology, that.homology);
    }

    @Override
    public int hashCode() {
        return Objects.hash(leftAlignedLeftBreakpoint, leftAlignedRightBreakpoint, insertedSequence, homology, inversionType);
    }
}
