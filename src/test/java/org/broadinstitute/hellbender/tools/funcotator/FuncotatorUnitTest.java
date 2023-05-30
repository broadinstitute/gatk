package org.broadinstitute.hellbender.tools.funcotator;

import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.testng.annotations.Test;

import java.util.LinkedHashSet;
import java.util.Set;

/**
 * Unit tests for the {@link org.broadinstitute.hellbender.tools.funcotator}.
 * Created by jonn on 8/22/17.
 */
public class FuncotatorUnitTest extends GATKBaseTest {

    /**
     * Unit test to check if our checkIfAlreadyAnnotated function works correctly.
     * In this case, we are using a header which does not have the
     * {@link FuncotatorConstants#FUNCOTATOR_VERSION_VCF_HEADERLINE_KEY}
     * in the vcf header, and thus, this function should throw an error because the given
     * vcf file has already been annotated and we do not need to annotate again.
     *
     */
    @Test(expectedExceptions = UserException.BadInput.class)
    public void testCheckWhenAlreadyAnnotated() {
        final Set<VCFHeaderLine> vcfHeaderLines = VCFHeader.makeHeaderVersionLineSet(VCFHeader.DEFAULT_VCF_VERSION);
        vcfHeaderLines.add(new VCFHeaderLine(
                FuncotatorConstants.FUNCOTATOR_VERSION_VCF_HEADERLINE_KEY,
                "content not relevant for this test"));
        final VCFHeader vcfHeader = new VCFHeader(vcfHeaderLines);
        Funcotator.checkIfAlreadyAnnotated(vcfHeader, null);
    }

    /**
     * Unit test to check if our checkIfAlreadyAnnotated function works correctly.
     * In this case, we are using the data provider which has the
     * {@link FuncotatorConstants#FUNCOTATOR_VERSION_VCF_HEADERLINE_KEY} header line
     * in the given vcf header, and thus, this function should
     * not throw an error because the given vcf file has not been annotated yet
     * and this file should be annotated.
     *
     */
    @Test
    public void testCheckWhenNotAlreadyAnnotated() {
        final VCFHeader vcfHeader = new VCFHeader(VCFHeader.makeHeaderVersionLineSet(VCFHeader.DEFAULT_VCF_VERSION));
        Funcotator.checkIfAlreadyAnnotated(vcfHeader, null);
        //this should not throw an exception.
    }




    //==================================================================================================================
    // Static Variables:

    //==================================================================================================================
    // Helper Methods:

    //==================================================================================================================
    // Data Providers:

    //==================================================================================================================
    // Tests:
}
