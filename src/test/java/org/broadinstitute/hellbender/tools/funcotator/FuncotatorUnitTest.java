package org.broadinstitute.hellbender.tools.funcotator;

import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

/**
 * Unit tests for the {@link org.broadinstitute.hellbender.tools.funcotator}.
 * Created by jonn on 8/22/17.
 */
public class FuncotatorUnitTest extends GATKBaseTest {

    /**
     * Data provider for the unit test for checkIfAlreadyAnnotated. In this case, we have a
     * "Funcotator" column in the vcf header, which indicates
     * that this vcf has already been annotated and thus the program
     * should not run on this input.
     *
     */
    @DataProvider
    Object[][] provideDataForAnnotationCheckWrong() {
        final Set<VCFHeaderLine> headerLines = new LinkedHashSet<>();

        headerLines.add(new VCFHeaderLine("fileformat", "fileformat"));
        headerLines.add(new VCFHeaderLine("FILTER", "FILTER"));
        headerLines.add(new VCFHeaderLine("FORMAT", "FORMAT"));
        headerLines.add(new VCFHeaderLine("Funcotator", "Funcotator"));
        headerLines.add(new VCFHeaderLine("GATKCommandLine", "GATKCommandLine"));
        headerLines.add(new VCFHeaderLine("INFO", "INFO"));

        return new Object[][] {
                { headerLines }
        };
    }

    /**
     * Data provider for the unit test for checkIfAlreadyAnnotated. In this case, we have a
     * "Funcotator Version" column in the vcf header, which indicates
     * that this vcf has already been annotated and thus the program
     * should not run on this input.
     *
     */
    @DataProvider
    Object[][] provideDataForAnnotationCheckWrongSecond() {
        final Set<VCFHeaderLine> headerLines = new LinkedHashSet<>();

        headerLines.add(new VCFHeaderLine("fileformat", "fileformat"));
        headerLines.add(new VCFHeaderLine("FILTER", "FILTER"));
        headerLines.add(new VCFHeaderLine("FORMAT", "FORMAT"));
        headerLines.add(new VCFHeaderLine("Funcotator Version", "Funcotator Version"));
        headerLines.add(new VCFHeaderLine("GATKCommandLine", "GATKCommandLine"));
        headerLines.add(new VCFHeaderLine("INFO", "INFO"));

        return new Object[][] {
                { headerLines }
        };
    }


    /**
     * Data provider for the unit test of our annotation function when the
     * given vcf file has not already been annotated. In this case, we see that
     * there is no longer a column labelled "Funcotator", and all other columns
     * are typically of a vcf file which has not been annotated yet. Thus,
     * the function should not throw an error and we should annotate this file.
     *
     */
    @DataProvider
    Object[][] provideDataForAnnotationCheckRight() {
        final Set<VCFHeaderLine> headerLines = new LinkedHashSet<>(); //= new List<>();
        //return Arrays.asList(
        headerLines.add(new VCFHeaderLine("fileformat", "VCF"));
        headerLines.add(new VCFHeaderLine("FILTER", "FILTER"));
        headerLines.add(new VCFHeaderLine("FORMAT", "FORMAT"));
        headerLines.add(new VCFHeaderLine("GATKCommandLine", "GATKCommandLine"));
        headerLines.add(new VCFHeaderLine("INFO", "INFO"));

        return new Object[][] {
                { headerLines }
        };
    }


    /**
     * Unit test to check if our checkIfAlreadyAnnotated function works correctly.
     * In this case, we are using the data provider which has an incorrect field name
     * in the vcf header, and thus, this function should throw an error because the given
     * vcf file has already been annotated and we do not need to annotate again.
     *
     */
    @Test(dataProvider = "provideDataForAnnotationCheckWrong", expectedExceptions = UserException.BadInput.class)
    public void testCheckWhenAlreadyAnnotated(final Set<VCFHeaderLine> vcfHeaderLines) {
        final VCFHeader vcfHeader = new VCFHeader(vcfHeaderLines);
        Funcotator.checkIfAlreadyAnnotated(vcfHeader, null);
    }

    /**
     * Unit test to check if our checkIfAlreadyAnnotated function works correctly.
     * In this case, we are using the data provider which has the second incorrect field name
     * in the vcf header, and thus, this function should throw an error because the given
     * vcf file has already been annotated and we do not need to annotate again.
     *
     */
    @Test(dataProvider = "provideDataForAnnotationCheckWrongSecond", expectedExceptions = UserException.BadInput.class)
    public void testCheckWhenAlreadyAnnotatedSecond(final Set<VCFHeaderLine> vcfHeaderLines) {
        final VCFHeader vcfHeader = new VCFHeader(vcfHeaderLines);
        Funcotator.checkIfAlreadyAnnotated(vcfHeader, null);
    }

    /**
     * Unit test to check if our checkIfAlreadyAnnotated function works correctly.
     * In this case, we are using the data provider which does not have any incorrect
     * field names in the given vcf header, and thus, this function should
     * not throw an error because the given vcf file has not been annotated yet
     * and this file should be annotated.
     *
     */
    @Test(dataProvider = "provideDataForAnnotationCheckRight")
    public void testCheckWhenNotAlreadyAnnotated(final Set<VCFHeaderLine> vcfHeaderLines) {
        final VCFHeader vcfHeader = new VCFHeader(vcfHeaderLines);
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
