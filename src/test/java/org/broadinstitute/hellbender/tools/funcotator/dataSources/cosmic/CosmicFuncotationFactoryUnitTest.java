package org.broadinstitute.hellbender.tools.funcotator.dataSources.cosmic;

import htsjdk.tribble.Feature;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.funcotator.Funcotation;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.TableFuncotation;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotationBuilder;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.testutils.FuncotatorReferenceTestUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.nio.file.Path;
import java.sql.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import static org.broadinstitute.hellbender.tools.funcotator.FuncotatorTestConstants.COSMIC_TEST_DB;
import static org.broadinstitute.hellbender.tools.funcotator.FuncotatorTestConstants.PIK3CA_POSITION;

/**
 * Class for running unit tests on {@link CosmicFuncotationFactory}.
 * Created by jonn on 12/17/17.
 */
public class CosmicFuncotationFactoryUnitTest extends GATKBaseTest {

    //==================================================================================================================
    // Private Static Members:

    private static final Path PATH_TO_TEST_DB = IOUtils.getPath(COSMIC_TEST_DB);

    private static final String PROTEIN_CHANGE_QUERY = "SELECT \"Mutation AA\" FROM Cosmic WHERE \"Mutation AA\" != \"\";";
    private static final String GENOMIC_POSITION_QUERY = "SELECT \"Mutation genome position\" FROM Cosmic WHERE \"Mutation genome position\" != \"\";";

    private static final VariantContext defaultVariantContext;
    private static final ReferenceContext defaultReferenceContext;

    private static final ReferenceDataSource PIK3CA_REF_DATA_SOURCE = ReferenceDataSource.of( new File(FuncotatorReferenceTestUtils.retrieveHg19Chr3Ref()).toPath() );

    static {
        final VariantContextBuilder variantContextBuilder = new VariantContextBuilder(
                FuncotatorReferenceTestUtils.retrieveHg19Chr3Ref(),
                PIK3CA_POSITION.getContig(),
                PIK3CA_POSITION.getStart(),
                PIK3CA_POSITION.getStart(),
                Arrays.asList(Allele.create("T", true), Allele.create("G"))
        );
        defaultVariantContext = variantContextBuilder.make();

        defaultReferenceContext = new ReferenceContext(
                PIK3CA_REF_DATA_SOURCE,
                PIK3CA_POSITION
        );
    }

    //==================================================================================================================
    // Private Members:

    //==================================================================================================================
    // Helper Methods:

    private VariantContext createVariantContext(final String contig,
                                                final int start,
                                                final int end,
                                                final String refString,
                                                final String altString) {

        final Allele refAllele = Allele.create(refString, true);
        final Allele altAllele = Allele.create(altString);

        final VariantContextBuilder variantContextBuilder = new VariantContextBuilder(
                FuncotatorReferenceTestUtils.retrieveHg19Chr3Ref(),
                contig,
                start,
                end,
                Arrays.asList(refAllele, altAllele)
        );

        return variantContextBuilder.make();
    }

    private Object[] helpProvideForTestCreateFuncotations(final String contig,
                                                          final int start,
                                                          final int end,
                                                          final String refAlleleString,
                                                          final String altAlleleString,
                                                          final List<GencodeFuncotation> gencodeFuncotationList,
                                                          final List<Funcotation> expected) {
        return new Object[] {
                createVariantContext(contig, start, end, refAlleleString, altAlleleString),
                new ReferenceContext( PIK3CA_REF_DATA_SOURCE, new SimpleInterval(contig, start, end)),
                Collections.emptyList(),
                gencodeFuncotationList,
                expected
        };
    }

    //==================================================================================================================
    // Data Providers:

    @DataProvider
    private Object[][] provideDataForTestProteinPositionRegex() {

        final ArrayList<String> proteinPositions = new ArrayList<>();
        final ArrayList<String> genomePositions = new ArrayList<>();

        // Get all protein positions in the test data:
        try {
            Class.forName("org.sqlite.JDBC");
            try ( final Connection dbConnection = DriverManager.getConnection("jdbc:sqlite:" + PATH_TO_TEST_DB.toUri().toString())) {
                try ( final Statement statement = dbConnection.createStatement() ) {

                    try ( final ResultSet resultSet = statement.executeQuery(PROTEIN_CHANGE_QUERY) ) {
                        while ( resultSet.next() ) {
                            proteinPositions.add( resultSet.getString("Mutation AA") );
                        }
                    }

                    try ( final ResultSet resultSet = statement.executeQuery(GENOMIC_POSITION_QUERY) ) {
                        while ( resultSet.next() ) {
                            genomePositions.add( resultSet.getString("Mutation genome position") );
                        }
                    }
                }
            }
        }
        catch (final SQLException ex) {
            throw new GATKException("Unable to read from database: " + PATH_TO_TEST_DB.toUri().toString());
        }
        catch (final ClassNotFoundException ex) {
            throw new GATKException("Cannot load SQLite Java Package!", ex);
        }

        return new Object[][] {
                { CosmicFuncotationFactory.PROTEIN_POSITION_REGEX, proteinPositions, 424 },
                { CosmicFuncotationFactory.GENOME_POSITION_REGEX, genomePositions, 346 }
        };
    }

    @DataProvider
    private Object[][] provideForTestCreateFuncotations() {
        return new Object[][] {
                // Trivial Case: No gencode funcotations:
                helpProvideForTestCreateFuncotations("chr3", 178916617, 178916617, "C", "A",
                        Collections.emptyList(),
                        Collections.singletonList(
                                TableFuncotation.create(Collections.singletonList("Cosmic_overlapping_mutations"), Collections.singletonList(""), Allele.create("A"), "Cosmic", null)
                        )
                ),
                // Trivial Case: Position is not in the database:
                helpProvideForTestCreateFuncotations("chr3", 178916617, 178916617, "C", "A",
                        Collections.singletonList(
                                new GencodeFuncotationBuilder().setHugoSymbol("PIK3CA").setChromosome("chr3").setStart(178916617).setEnd(178916617).setProteinChange("p.P2T").build()
                        ),
                        Collections.singletonList(
                                TableFuncotation.create(Collections.singletonList("Cosmic_overlapping_mutations"), Collections.singletonList(""), Allele.create("A"), "Cosmic", null)
                        )
                ),
                // Protein position match:
                // NOTE: Intentional bad values in the variant position:
                helpProvideForTestCreateFuncotations("chr3", 1, 1, "G", "A",
                        Collections.singletonList(
                                new GencodeFuncotationBuilder().setHugoSymbol("PIK3CA").setChromosome("chr3").setStart(178936091).setEnd(178936091).setProteinChange("p.E545K").build()
                        ),
                        Collections.singletonList(
                                TableFuncotation.create(Collections.singletonList("Cosmic_overlapping_mutations"), Collections.singletonList("p.E545K(2)"), Allele.create("A"), "Cosmic", null)
                        )
                ),
                // Genome position match:
                helpProvideForTestCreateFuncotations("chr3", 178936091, 178936091, "G", "A",
                        Collections.singletonList(
                                new GencodeFuncotationBuilder().setHugoSymbol("PIK3CA").setChromosome("chr3").setStart(178936091).setEnd(178936091).setProteinChange("p.E999K").build()
                        ),
                        Collections.singletonList(
                                TableFuncotation.create(Collections.singletonList("Cosmic_overlapping_mutations"), Collections.singletonList("p.E545K(1)"), Allele.create("A"), "Cosmic", null)
                        )
                ),
                // Protein and Genome position match:
                helpProvideForTestCreateFuncotations("chr3", 178936091, 178936091, "G", "A",
                        Collections.singletonList(
                                new GencodeFuncotationBuilder().setHugoSymbol("PIK3CA").setChromosome("chr3").setStart(178936091).setEnd(178936091).setProteinChange("p.E545K").build()
                        ),
                        Collections.singletonList(
                                TableFuncotation.create(Collections.singletonList("Cosmic_overlapping_mutations"), Collections.singletonList("p.E545K(2)"), Allele.create("A"), "Cosmic", null)
                        )
                ),
                // Huge Protein Match - all in the file:
                helpProvideForTestCreateFuncotations("chr3", 1, 1, "G", "A",
                        Collections.singletonList(
                                new GencodeFuncotationBuilder().setHugoSymbol("PIK3CA").setChromosome("chr3").setStart(178936091).setEnd(178936091).setProteinChange("p.E1_K3455").build()
                        ),
                        Collections.singletonList(
                                TableFuncotation.create(Collections.singletonList("Cosmic_overlapping_mutations"), Collections.singletonList("p.E545K(2)|p.E542K(2)|p.H1047R(2)|p.N345K(1)"), Allele.create("A"), "Cosmic", null)
                        )
                ),
                // Huge Genome Position Match - all in the file:
                helpProvideForTestCreateFuncotations("chr3", 178921553, 178952085,
                        String.join("", Collections.nCopies(178952085 - 178921553 + 1, "G")),
                        String.join("", Collections.nCopies(178952085 - 178921553 + 1, "A")),
                        Collections.singletonList(
                                new GencodeFuncotationBuilder().setHugoSymbol("PIK3CA").setChromosome("chr3").setStart(178921553).setEnd(178952085).setProteinChange("p.E1T").build()
                        ),
                        Collections.singletonList(
                                TableFuncotation.create(Collections.singletonList("Cosmic_overlapping_mutations"), Collections.singletonList("p.E542K(2)|p.H1047R(1)|p.E545K(1)|p.N345K(1)"), Allele.create(String.join("", Collections.nCopies(178952085 - 178921553 + 1, "A"))), "Cosmic", null)
                        )
                ),
        };
    }

    //==================================================================================================================
    // Tests:

    @Test(dataProvider = "provideDataForTestProteinPositionRegex")
    public void testPositionRegex(final Pattern regex, final List<String> dbPositions, final int expectedNumResults) {

        int numMatches = 0;
        for ( final String pos : dbPositions ) {
            final Matcher matcher = regex.matcher(pos);
            if ( matcher.matches() ) {
                ++numMatches;
            }
            else {
                System.out.println("");
            }
        }

        Assert.assertEquals(numMatches, expectedNumResults);
    }

    @Test
    public void testGetName() {
        final CosmicFuncotationFactory cosmicFuncotationFactory = new CosmicFuncotationFactory(PATH_TO_TEST_DB);
        Assert.assertEquals(cosmicFuncotationFactory.getName(), "Cosmic");
    }

    @Test
    public void testGetSupportedFuncotationFields() {
        final CosmicFuncotationFactory cosmicFuncotationFactory = new CosmicFuncotationFactory(PATH_TO_TEST_DB);

        final LinkedHashSet<String> expectedFieldNames = new LinkedHashSet<>(1);
        expectedFieldNames.add("Cosmic_overlapping_mutations");

        Assert.assertEquals(cosmicFuncotationFactory.getSupportedFuncotationFields(), expectedFieldNames);
    }

    @Test(expectedExceptions = GATKException.class)
    public void testCreateFuncotationsNoGencodeInput() {
        final CosmicFuncotationFactory cosmicFuncotationFactory = new CosmicFuncotationFactory(PATH_TO_TEST_DB);

        cosmicFuncotationFactory.createFuncotationsOnVariant(defaultVariantContext, defaultReferenceContext, Collections.emptyList());
    }

    @Test(dataProvider = "provideForTestCreateFuncotations")
    public void testCreateFuncotations(final VariantContext variant,
                                       final ReferenceContext referenceContext,
                                       final List<Feature> featureList,
                                       final List<GencodeFuncotation> gencodeFuncotations,
                                       final List<Funcotation> expected) {

        final CosmicFuncotationFactory cosmicFuncotationFactory = new CosmicFuncotationFactory(PATH_TO_TEST_DB);

        Assert.assertEquals(
            cosmicFuncotationFactory.createFuncotationsOnVariant(
                variant,
                referenceContext,
                featureList,
                gencodeFuncotations
            ),
            expected
        );
    }
}
