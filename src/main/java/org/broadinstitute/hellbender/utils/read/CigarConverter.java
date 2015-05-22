package org.broadinstitute.hellbender.utils.read;

import com.google.api.services.genomics.model.CigarUnit;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

import java.util.*;

public final class CigarConverter {

    private static final List<String> CIGAR_UNIT_OPERATORS = Arrays.asList("ALIGNMENT_MATCH", "INSERT", "DELETE", "SKIP", "CLIP_SOFT", "CLIP_HARD", "PAD", "SEQUENCE_MATCH", "SEQUENCE_MISMATCH");
    private static final List<CigarOperator> SAM_CIGAR_ELEMENT_OPERATORS = Arrays.asList(CigarOperator.M, CigarOperator.I, CigarOperator.D, CigarOperator.N, CigarOperator.S, CigarOperator.H, CigarOperator.P, CigarOperator.EQ, CigarOperator.X);
    private static final Map<String, CigarOperator> cigarUnitOperatorToCigarElementOperatorTable;
    private static final Map<CigarOperator, String> cigarElementOperatorToCigarUnitOperatorTable;

    static {
        cigarUnitOperatorToCigarElementOperatorTable = new HashMap<>(CIGAR_UNIT_OPERATORS.size() * 2);
        cigarElementOperatorToCigarUnitOperatorTable = new HashMap<>(CIGAR_UNIT_OPERATORS.size() * 2);

        for ( int i = 0; i < CIGAR_UNIT_OPERATORS.size(); ++i ) {
            cigarUnitOperatorToCigarElementOperatorTable.put(CIGAR_UNIT_OPERATORS.get(i), SAM_CIGAR_ELEMENT_OPERATORS.get(i));
            cigarElementOperatorToCigarUnitOperatorTable.put(SAM_CIGAR_ELEMENT_OPERATORS.get(i), CIGAR_UNIT_OPERATORS.get(i));
        }
    }

    public static Cigar convertCigarUnitListToSAMCigar( final List<CigarUnit> cigarUnits ) {
        final List<CigarElement> cigarElements = new ArrayList<>(cigarUnits.size());
        for ( CigarUnit cigarUnit : cigarUnits ) {
            cigarElements.add(convertCigarUnitToSAMCigarElement(cigarUnit));
        }
        return new Cigar(cigarElements);
    }

    public static List<CigarUnit> convertSAMCigarToCigarUnitList( final Cigar cigar ) {
        final List<CigarUnit> cigarUnits = new ArrayList<>(cigar.numCigarElements());
        for ( CigarElement cigarElement : cigar.getCigarElements() ) {
            cigarUnits.add(convertSAMCigarElementToCigarUnit(cigarElement));
        }
        return cigarUnits;
    }

    public static CigarElement convertCigarUnitToSAMCigarElement( final CigarUnit cigarUnit ) {
        if ( cigarUnit.getOperationLength() == null || cigarUnit.getOperation() == null ) {
            throw new IllegalArgumentException("Cannot convert a CigarUnit with a missing operation or operation length to a CigarElement");
        }

        return new CigarElement(cigarUnit.getOperationLength().intValue(),
                cigarUnitOperatorToCigarElementOperatorTable.get(cigarUnit.getOperation()));
    }

    public static CigarUnit convertSAMCigarElementToCigarUnit( final CigarElement cigarElement ) {
        final CigarUnit cigarUnit = new CigarUnit();

        cigarUnit.setOperationLength((long)cigarElement.getLength());
        cigarUnit.setOperation(cigarElementOperatorToCigarUnitOperatorTable.get(cigarElement.getOperator()));
        return cigarUnit;
    }
}