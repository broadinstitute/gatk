package org.broadinstitute.hellbender.utils.read;

import com.google.api.services.genomics.model.CigarUnit;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import org.broadinstitute.hellbender.exceptions.GATKException;

import java.util.*;

/**
 * Utilities for bidirectional conversion between SAM Cigar and GA4GH-style List&lt;CigarUnit&gt;.
 *
 * Can convert either entire Cigars or individual Cigar elements.
 *
 * Conversion from List&lt;CigarUnit&gt; to SAM Cigar can be lossy, since there is no equivalent to
 * CigarUnit's referenceSequence field in SAM Cigars. This field is silently dropped during conversion,
 * if present.
 */
public final class CigarConversionUtils {

    private CigarConversionUtils() {}

    /**
     * Valid cigar operators for CigarUnits as defined in the GA4GH spec. The order of these operators matters,
     * and must match the ordering in {@link #SAM_CIGAR_ELEMENT_OPERATORS} below.
     */
    private static final List<String> CIGAR_UNIT_OPERATORS =
            Arrays.asList("ALIGNMENT_MATCH",
                          "INSERT",
                          "DELETE",
                          "SKIP",
                          "CLIP_SOFT",
                          "CLIP_HARD",
                          "PAD",
                          "SEQUENCE_MATCH",
                          "SEQUENCE_MISMATCH");

    /**
     * Valid cigar operators for SAM Cigars as defined in the SAM spec. The order of these operators matters,
     * and must match the ordering in {@link #CIGAR_UNIT_OPERATORS} above.
     */
    private static final List<CigarOperator> SAM_CIGAR_ELEMENT_OPERATORS =
            Arrays.asList(CigarOperator.M,
                          CigarOperator.I,
                          CigarOperator.D,
                          CigarOperator.N,
                          CigarOperator.S,
                          CigarOperator.H,
                          CigarOperator.P,
                          CigarOperator.EQ,
                          CigarOperator.X);

    /**
     * Lookup table from GA4GH CigarUnit operator -> equivalent SAM CigarOperator
     */
    private static final Map<String, CigarOperator> ga4ghToSAMOperatorTable;

    /**
     * Lookup table from SAM CigarOperator -> equivalent GA4GH CigarUnit operator
     */
    private static final Map<CigarOperator, String> samToGA4GHOperatorTable;

    static {
        // Populate cigar operator lookup tables

        Map<String, CigarOperator> ga4ghToSAMOperatorTableLocal = new HashMap<>(CIGAR_UNIT_OPERATORS.size() * 2);
        Map<CigarOperator, String> samToGA4GHOperatorTableLocal = new HashMap<>(CIGAR_UNIT_OPERATORS.size() * 2);

        for ( int i = 0; i < CIGAR_UNIT_OPERATORS.size(); ++i ) {
            ga4ghToSAMOperatorTableLocal.put(CIGAR_UNIT_OPERATORS.get(i), SAM_CIGAR_ELEMENT_OPERATORS.get(i));
            samToGA4GHOperatorTableLocal.put(SAM_CIGAR_ELEMENT_OPERATORS.get(i), CIGAR_UNIT_OPERATORS.get(i));
        }

        ga4ghToSAMOperatorTable = Collections.unmodifiableMap(ga4ghToSAMOperatorTableLocal);
        samToGA4GHOperatorTable = Collections.unmodifiableMap(samToGA4GHOperatorTableLocal);
    }

    /**
     * Converts a GA4GH-style List&lt;CigarUnit&gt; to a SAM Cigar.
     *
     * This conversion can be lossy, since there is no equivalent to CigarUnit's referenceSequence field
     * in SAM Cigars. This field is silently dropped during conversion, if present.
     *
     * @param cigarUnits List&lt;CigarUnit&gt; to convert to SAM Cigar.
     * @return equivalent SAM Cigar
     */
    public static Cigar convertCigarUnitListToSAMCigar( final List<CigarUnit> cigarUnits ) {
        if ( cigarUnits == null ) {
            throw new IllegalArgumentException("Cannot convert null CigarUnit list");
        }

        final List<CigarElement> cigarElements = new ArrayList<>(cigarUnits.size());
        for ( CigarUnit cigarUnit : cigarUnits ) {
            cigarElements.add(convertCigarUnitToSAMCigarElement(cigarUnit));
        }
        return new Cigar(cigarElements);
    }

    /**
     * Converts a SAM Cigar to a GA4GH-style List&lt;CigarUnit&gt;
     *
     * @param cigar SAM Cigar to convert
     * @return equivalent List&lt;CigarUnit&gt;
     */
    public static List<CigarUnit> convertSAMCigarToCigarUnitList( final Cigar cigar ) {
        if ( cigar == null ) {
            throw new IllegalArgumentException("Cannot convert null Cigar");
        }

        final List<CigarUnit> cigarUnits = new ArrayList<>(cigar.numCigarElements());
        for ( CigarElement cigarElement : cigar.getCigarElements() ) {
            cigarUnits.add(convertSAMCigarElementToCigarUnit(cigarElement));
        }
        return cigarUnits;
    }

    /**
     * Converts an individual GA4GH-style CigarUnit to an equivalent SAM CigarElement
     *
     * This conversion can be lossy, since there is no equivalent to CigarUnit's referenceSequence field
     * in SAM CigarElement. This field is silently dropped during conversion, if present.
     *
     * @param cigarUnit CigarUnit to convert to SAM CigarElement
     * @return equivalent SAM CigarElement
     */
    public static CigarElement convertCigarUnitToSAMCigarElement( final CigarUnit cigarUnit ) {
        if ( cigarUnit == null || cigarUnit.getOperationLength() == null || cigarUnit.getOperation() == null ) {
            throw new IllegalArgumentException("Cannot convert a null CigarUnit or a CigarUnit with a missing operation or operation length to a CigarElement");
        }

        final CigarOperator convertedOperation = ga4ghToSAMOperatorTable.get(cigarUnit.getOperation());
        if ( convertedOperation == null ) {
            throw new GATKException("Unable to convert CigarUnit to SAM CigarElement: no SAM-equivalent cigar operator found for CigarUnit operation " + cigarUnit.getOperation());
        }

        if ( cigarUnit.getOperationLength() < 0 || cigarUnit.getOperationLength() > Integer.MAX_VALUE ) {
            throw new GATKException("Unable to convert CigarUnit to SAM CigarElement: CigarUnit operation length " + cigarUnit.getOperationLength() + " must be between 0 and Integer.MAX_VALUE, inclusive");
        }

        return new CigarElement(cigarUnit.getOperationLength().intValue(), convertedOperation);
    }

    /**
     * Converts an individual SAM CigarElement to an equivalent GA4GH-style CigarUnit
     *
     * @param cigarElement SAM CigarElement to convert to CigarUnit
     * @return equivalent CigarUnit
     */
    public static CigarUnit convertSAMCigarElementToCigarUnit( final CigarElement cigarElement ) {
        if ( cigarElement == null ) {
            throw new IllegalArgumentException("Cannot convert a null CigarElement");
        }

        final String convertedOperation = samToGA4GHOperatorTable.get(cigarElement.getOperator());
        if ( convertedOperation == null ) {
            throw new GATKException("Unable to convert CigarElement to CigarUnit: unknown operator in CigarElement: " + cigarElement.getOperator());
        }

        if ( cigarElement.getLength() < 0 ) {
            throw new GATKException("Unable to convert CigarElement to CigarUnit: negative operation length in CigarElement: " + cigarElement.getLength());
        }

        final CigarUnit cigarUnit = new CigarUnit();
        cigarUnit.setOperationLength((long)cigarElement.getLength());
        cigarUnit.setOperation(convertedOperation);

        return cigarUnit;
    }
}
