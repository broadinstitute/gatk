package org.broadinstitute.hellbender.utils;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.testutils.FuncotatorReferenceTestUtils;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotationFactory;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.annotations.AfterClass;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class TranscriptUtilsUnitTest extends GATKBaseTest {

    //==================================================================================================================
    // Private Static Members:

    private static final ReferenceDataSource refDataSourceHg19Ch3;

    private static final List<AutoCloseable> autoCloseableList = new ArrayList<>();
    static {
        refDataSourceHg19Ch3 = ReferenceDataSource.of( IOUtils.getPath(FuncotatorReferenceTestUtils.retrieveHg19Chr3Ref()) );
        autoCloseableList.add( refDataSourceHg19Ch3 );
    }

    //==================================================================================================================
    // Setup and Breakdown Methods with Hooks:

    @BeforeClass
    public static void setupBeforeTests() {
        System.out.println("Setting up before tests...");
    }

    @AfterClass
    public static void cleanupAfterTests() {
        System.out.println("Cleaning up after tests...");

        for ( final AutoCloseable cl : autoCloseableList ) {

            try {
                cl.close();
            }
            catch ( final Exception ex ) {
                throw new GATKException("Could not close " + cl.toString(), ex);
            }
        }
    }

    //==================================================================================================================
    // Private Members:

    //==================================================================================================================
    // Helper Methods:

    //==================================================================================================================
    // Data Providers:

    @DataProvider(name = "provideForTestExtractTrascriptFromReference")
    Object[][] provideForTestExtractTrascriptFromReference() {

        final SimpleInterval cntn4Interval = new SimpleInterval("chr3", 2140497, 3099645);
        final ReferenceContext refContextCntn4 = new ReferenceContext( refDataSourceHg19Ch3, cntn4Interval);

        return new Object[][] {
                { refContextCntn4, List.of(new SimpleInterval("chr3", 2140497, 2140497)), "A" },
                { refContextCntn4, List.of(new SimpleInterval("chr3", 2140497, 2140498)), "AG" },
                { refContextCntn4, List.of(new SimpleInterval("chr3", 2140497, 2140499)), "AGC" },
                { refContextCntn4, List.of(new SimpleInterval("chr3", 2140497, 2140500)), "AGCC" },
                { refContextCntn4, List.of(new SimpleInterval("chr3", 2140497, 2140501)), "AGCCG" },
                { refContextCntn4, List.of(new SimpleInterval("chr3", 2140498, 2140501)), "GCCG" },
                { refContextCntn4, List.of(new SimpleInterval("chr3", 2140499, 2140501)), "CCG" },
                { refContextCntn4, List.of(new SimpleInterval("chr3", 2140500, 2140501)), "CG" },
                { refContextCntn4, List.of(new SimpleInterval("chr3", 2140501, 2140501)), "G" },
        };
    }

    //==================================================================================================================
    // Tests:

    @Test(dataProvider = "provideForTestExtractTrascriptFromReference")
    public void TestExtractTrascriptFromReference(final ReferenceContext refContext, final List<SimpleInterval> exons, final String expectedTranscript) {
        final String transcript = TranscriptUtils.extractTrascriptFromReference(refContext, exons);
        Assert.assertEquals(transcript, expectedTranscript);
    }

}
