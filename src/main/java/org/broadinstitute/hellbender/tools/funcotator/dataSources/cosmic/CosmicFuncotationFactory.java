package org.broadinstitute.hellbender.tools.funcotator.dataSources.cosmic;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.tribble.Feature;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.funcotator.DataSourceFuncotationFactory;
import org.broadinstitute.hellbender.tools.funcotator.Funcotation;
import org.broadinstitute.hellbender.tools.funcotator.FuncotatorArgumentDefinitions;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.TableFuncotation;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.sqlite.SQLiteConfig;

import java.nio.file.Path;
import java.sql.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;


/**
 * Factory for creating {@link Funcotation}s by handling a SQLite database containing information from COSMIC.
 * The raw datasource (http://cancer.sanger.ac.uk/cosmic/download - CosmicCompleteTargetedScreensMutantExport.tsv.gz)
 * must be unzipped and preprocessed with the script `createSqliteCosmicDb.sh`.
 *
 *
 * This is a high-level object that interfaces with the internals of {@link org.broadinstitute.hellbender.tools.funcotator.Funcotator}.
 * Created by jonn on 12/16/17.
 */
public class CosmicFuncotationFactory extends DataSourceFuncotationFactory {

    //==================================================================================================================
    // Public Static Members:

    //==================================================================================================================
    // Private Static Members:

    /** Standard Logger.  */
    protected static final Logger logger = LogManager.getLogger(CosmicFuncotationFactory.class);

    /**
     * Regular expression for matching a genome position in a Cosmic record.
     */
    @VisibleForTesting
    static final Pattern GENOME_POSITION_REGEX =  Pattern.compile("(\\d+):(\\d+)-(\\d+)");

    /**
     * Regular expression for matching a protein position in a Cosmic record.
     */
    @VisibleForTesting
    static final Pattern PROTEIN_POSITION_REGEX =  Pattern.compile("[pP]\\.\\(?[A-Z](\\d+)\\)?[A-Z]?_?(?:[A-Z]?(\\d+)?.*)?");

    /**
     * Placeholder contig name for protein positions.
     */
    private static final String PROTEIN_CONTIG = "P";

    /** Name of the Cosmic table in the DB. */
    private static final String TABLE_NAME = "Cosmic";

    /**
     * The name of the column containing protein position in the DB.
     */
    private static final String PROTEIN_POSITION_COLUMN_NAME = "Mutation AA";

    /**
     * The name of the column containing genome position in the DB.
     */
    private static final String GENOME_POSITION_COLUMN_NAME = "Mutation genome position";

    /**
     * The name of the column containing gene name in the DB.
     */
    private static final String GENE_NAME_COLUMN = "Gene name";

    /**
     * Fields to ignore when returning results from this {@link CosmicFuncotationFactory}.
     */
    private static final HashSet<String> IGNORE_FIELDS = new HashSet<>(Arrays.asList(PROTEIN_POSITION_COLUMN_NAME, GENOME_POSITION_COLUMN_NAME, GENE_NAME_COLUMN));


    /** Query to get the field names from the DB */
    private static final String FIELD_NAME_QUERY = "SELECT * FROM " + TABLE_NAME + " LIMIT 1;";

    /**
     * Template for results query for matching genes in the database.
     */
    private static final String RESULT_QUERY_TEMPLATE = "SELECT * FROM " + TABLE_NAME + " WHERE \""
            + GENE_NAME_COLUMN + "\" == ";

    //==================================================================================================================
    // Private Members:

    /**
     * The name of this {@link CosmicFuncotationFactory}.
     */
    private final String name = "Cosmic";

    /**
     * The path to the Sql
     */
    private final Path pathToCosmicDb;

    /**
     * The connection to the SQLite database for this {@link CosmicFuncotationFactory}.
     */
    private final Connection dbConnection;

    /**
     * The ordered set of fields that this {@link CosmicFuncotationFactory} supports.
     */
    private final LinkedHashSet<String> supportedFields;

    //==================================================================================================================
    // Constructors:

    public CosmicFuncotationFactory(final Path pathToCosmicDb) {
        this(pathToCosmicDb, new LinkedHashMap<>(), DEFAULT_VERSION_STRING);
    }

    public CosmicFuncotationFactory(final Path pathToCosmicDb,
                                    final LinkedHashMap<String, String> annotationOverridesMap,
                                    final String version) {
        this.pathToCosmicDb = pathToCosmicDb;

        this.version = version;

        // Connect to the DB:
        try {
            Class.forName("org.sqlite.JDBC");

            // Set our configuration options for the connection here:
            final SQLiteConfig config = new SQLiteConfig();
            // We only want to read from the DB:
            config.setReadOnly(true);

            logger.debug("Connecting to SQLite database at: " + pathToCosmicDb.toUri().toString());
            dbConnection = DriverManager.getConnection("jdbc:sqlite:" + pathToCosmicDb.toUri().toString(), config.toProperties());
            logger.debug("Connected to SQLite database!");
        }
        catch (final SQLException ex) {
            throw new UserException("Unable to open SQLite DB for COSMIC at: " + pathToCosmicDb.toUri().toString(), ex);
        }
        catch (final ClassNotFoundException ex) {
            throw new UserException("Cannot load SQLite Java Package!", ex);
        }

        // Get the supported fields:
        supportedFields = new LinkedHashSet<>(1);
        supportedFields.add(name + "_overlapping_mutations");

        // Initialize our annotation overrides:
        initializeAnnotationOverrides(annotationOverridesMap);
    }

    //==================================================================================================================
    // Override Methods:

    @Override
    public void close() {
        if (dbConnection != null) {
            try {
                dbConnection.close();
            }
            catch (final SQLException ex) {
                throw new GATKException("Unable to close the connection to DB: " + pathToCosmicDb.toUri().toString(), ex);
            }
        }
    }

    @Override
    public String getName() {
        return name;
    }

    @Override
    public LinkedHashSet<String> getSupportedFuncotationFields() {
        return supportedFields;
    }

    @Override
    /**
     * {@inheritDoc}
     * This method should never be called on a {@link CosmicFuncotationFactory} - knowledge of the applied
     * {@link GencodeFuncotation}s is required to create a {@link Funcotation} from here.
     */
    public List<Funcotation> createFuncotations(final VariantContext variant, final ReferenceContext referenceContext, final List<Feature> featureList) {
        // TODO: this should be allowed, but with a warning to the user that only annotations with good Genome Positions can be used.
        throw new GATKException(this.getClass().getName() + " requires a set of GencodeFuncotations in order to createFuncotations!  This method should never be called on a " + this.getClass().getName());
    }

    @Override
    public List<Funcotation> createFuncotations(final VariantContext variant, final ReferenceContext referenceContext, final List<Feature> featureList, final List<GencodeFuncotation> gencodeFuncotations) {

        final List<Funcotation> outputFuncotations = new ArrayList<>();

        // Keep count of our overlapping mutations here:
        int numOverlappingMutations = 0;

        // If we have gencodeFuncotations we go through them and get the gene name
        // Then query our DB for matches on the gene name.
        // Then grab Genome position / Protein position and see if we overlap.
        // If any do, we create our CosmicFuncotation
        for (  final GencodeFuncotation gencodeFuncotation : gencodeFuncotations ) {
            final String geneName = gencodeFuncotation.getHugoSymbol();

            final SimpleInterval genomePosition = new SimpleInterval(variant.getContig(), variant.getStart(), variant.getEnd());

            final SimpleInterval proteinPosition;
            if ( gencodeFuncotation.getProteinChange() != null ) {
                proteinPosition = parseProteinString(gencodeFuncotation.getProteinChange());
            }
            else {
                proteinPosition = null;
            }

            try {
                try ( final Statement statement = dbConnection.createStatement() ) {
                    try ( final ResultSet resultSet = statement.executeQuery(RESULT_QUERY_TEMPLATE + "\"" + geneName + "\";") ) {
                        // iterate through our results:
                        while ( resultSet.next() ) {
                            // Get the genome position and protein position:
                            final SimpleInterval cosmicGenomePosition = getGenomePositionFromResults(resultSet);
                            final SimpleInterval cosmicProteinPosition = getProteinPositionFromResults(resultSet);

                            // Try to match on genome position first:
                            if ( cosmicGenomePosition != null ) {
                                // If we overlap the records, we update the counter:
                                if ( genomePosition.overlaps(cosmicGenomePosition) ) {
                                    ++numOverlappingMutations;
                                    continue;
                                }
                            }

                            // Now try to match on protein position:
                            if ( proteinPosition != null ) {
                                // If we overlap the records, we update the counter:
                                if ( proteinPosition.overlaps(cosmicProteinPosition) ) {
                                    ++numOverlappingMutations;
                                }
                            }
                            // NOTE: We can't annotate if the protein position and the genome position are null.
                        }
                    }
                }
            }
            catch (final SQLException ex) {
                throw new GATKException("Unable to query the database for geneName: " + geneName, ex);
            }
        }

        // Add our tally for all alternate alleles in this variant:
        for ( final Allele altAllele : variant.getAlternateAlleles() ) {
            outputFuncotations.add(
                    new TableFuncotation(
                            new ArrayList<>(supportedFields),
                            new ArrayList<>(Collections.singletonList(String.valueOf(numOverlappingMutations))),
                            altAllele,
                            name
                    )
            );
        }

        setOverrideValuesInFuncotations(outputFuncotations);

        return outputFuncotations;
    }

    @Override
    public FuncotatorArgumentDefinitions.DataSourceType getType() {
        return FuncotatorArgumentDefinitions.DataSourceType.COSMIC;
    }

    //==================================================================================================================
    // Static Methods:

    //==================================================================================================================
    // Instance Methods:

    /**
     * Get the genome position of the current record in the given {@link ResultSet}.
     * @param resultSet The results of a query on the database with a current row (must not be {@code null}).
     * @return A {@link SimpleInterval} represnting the genome position of the current record in the given {@link ResultSet}; or {@code null}.
     */
    private final SimpleInterval getGenomePositionFromResults(final ResultSet resultSet) {
        Utils.nonNull(resultSet);

        try {
            final String rawPosition = resultSet.getString(GENOME_POSITION_COLUMN_NAME);
            final Matcher matcher = GENOME_POSITION_REGEX.matcher(rawPosition);
            if ( matcher.matches() ) {
                // We have a position, so we should parse it:
                final String rawContig =  matcher.group(1);
                final String contig;
                if ( rawContig.startsWith("chr") ) {
                    contig = rawContig;
                }
                else {
                    contig = "chr" + rawContig;
                }
                final int start = Integer.valueOf(matcher.group(2));
                final int end = Integer.valueOf(matcher.group(3));

                return new SimpleInterval(contig, start, end);
            }
        }
        catch (final SQLException ex) {
            throw new GATKException("Cannot get Genome Position from column: " + GENOME_POSITION_COLUMN_NAME, ex);
        }

        return null;
    }

    /**
     * Pulls a protein change / protein position out of the current record in the given {@link ResultSet}.
     * @param resultSet The results of a query on the database with a current row (must not be {@code null}).
     * @return A {@link SimpleInterval} representing the extents of the protein position in the current record in the given {@link ResultSet} or {@code null}.
     */
    private final SimpleInterval getProteinPositionFromResults(final ResultSet resultSet) {
        Utils.nonNull(resultSet);

        try {
            final String rawPosition = resultSet.getString(PROTEIN_POSITION_COLUMN_NAME);
            return parseProteinString(rawPosition);
        }
        catch (final SQLException ex) {
            throw new GATKException("Cannot get Protein Position from column: " + GENOME_POSITION_COLUMN_NAME, ex);
        }
    }

    /**
     * Parse a {@link SimpleInterval} from a protein position / protein change.
     * @param proteinPositionString A {@link String} representing a protein position / protein change.
     * @return A {@link SimpleInterval} representing the extents of the given {@code proteinPositionString} or {@code null}.
     */
    private final SimpleInterval parseProteinString(final String proteinPositionString) {
        Utils.nonNull(proteinPositionString);

        final Matcher matcher = PROTEIN_POSITION_REGEX.matcher(proteinPositionString);
        if ( matcher.matches() ) {
            // We have a position, so we should parse it.
            // If we have an end position, we should use it:
            if ( matcher.group(2) != null ) {
                return new SimpleInterval(PROTEIN_CONTIG, Integer.valueOf(matcher.group(1)), Integer.valueOf(matcher.group(2)));
            }
            else {
                return new SimpleInterval(PROTEIN_CONTIG, Integer.valueOf(matcher.group(1)), Integer.valueOf(matcher.group(1)));
            }
        }
        return null;
    }

    //==================================================================================================================
    // Helper Data Types:

}
