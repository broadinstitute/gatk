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
import java.util.stream.Collectors;


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
    protected Class<? extends Feature> getAnnotationFeatureClass() {
        // Returning Feature.class here implies that this class doesn't care about what features it gets.
        return Feature.class;
    }

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
    protected List<Funcotation> createDefaultFuncotationsOnVariant( final VariantContext variant, final ReferenceContext referenceContext ) {

        final List<Funcotation> funcotationList = new ArrayList<>();

        // Add our tally for all alternate alleles in this variant:
        for ( final Allele altAllele : variant.getAlternateAlleles() ) {
            funcotationList.add(
                    TableFuncotation.create(
                            new ArrayList<>(supportedFields),
                            new ArrayList<>(Collections.singletonList("")),
                            altAllele,
                            name, null
                    )
            );
        }

        return funcotationList;
    }

    @Override
    /**
     * {@inheritDoc}
     * This method should never be called on a {@link CosmicFuncotationFactory} - knowledge of the applied
     * {@link GencodeFuncotation}s is required to create a {@link Funcotation} from here.
     */
    protected List<Funcotation> createFuncotationsOnVariant(final VariantContext variant, final ReferenceContext referenceContext, final List<Feature> featureList) {
        // TODO: this should be allowed, but with a warning to the user that only annotations with good Genome Positions can be used.  Issue: https://github.com/broadinstitute/gatk/issues/4594
        throw new GATKException(this.getClass().getName() + " requires a set of GencodeFuncotations in order to createFuncotationsOnVariant!  This method should never be called on a " + this.getClass().getName());
    }

    @Override
    protected List<Funcotation> createFuncotationsOnVariant(final VariantContext variant, final ReferenceContext referenceContext, final List<Feature> featureList, final List<GencodeFuncotation> gencodeFuncotations) {

        final List<Funcotation> outputFuncotations = new ArrayList<>();

        // Keep count of each overlapping mutation here:
        final Map<String, Integer> proteinChangeCounts = new LinkedHashMap<>();

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

                            // Get the genome position:
                            final SimpleInterval cosmicGenomePosition = getGenomePositionFromResults(resultSet);

                            // Try to match on genome position first:
                            if ( cosmicGenomePosition != null ) {
                                if ( genomePosition.overlaps(cosmicGenomePosition) ) {
                                    // If we overlap the records, we get the protein change and add it to the map:
                                    updateProteinChangeCountMap(proteinChangeCounts, resultSet);
                                    continue;
                                }
                            }

                            // Get the protein position:
                            final SimpleInterval cosmicProteinPosition = getProteinPositionFromResults(resultSet);

                            // Now try to match on protein position:
                            if ( proteinPosition != null ) {
                                // If we overlap the records, we update the counter:
                                if ( proteinPosition.overlaps(cosmicProteinPosition) ) {
                                    updateProteinChangeCountMap(proteinChangeCounts, resultSet);
                                }
                            }
                            // NOTE: We can't annotate if the protein position is null.
                        }
                    }
                }
            }
            catch (final SQLException ex) {
                throw new GATKException("Unable to query the database for geneName: " + geneName, ex);
            }
        }

        // Add our counts to all alternate alleles in this variant:
        for ( final Allele altAllele : variant.getAlternateAlleles() ) {
            outputFuncotations.add(
                    TableFuncotation.create(
                            new ArrayList<>(supportedFields),
                            Collections.singletonList(proteinChangeCounts.entrySet().stream()
                                    .map(entry -> entry.getKey() + '('+ entry.getValue() + ')')
                                    .collect(Collectors.joining("|"))),
                            altAllele,
                            name, null
                    )
            );
        }

        return outputFuncotations;
    }

    private void updateProteinChangeCountMap(final Map<String, Integer> proteinChangeCounts, final ResultSet resultSet) {
        final String proteinChange = getProteinChangeStringFromResults(resultSet);
        if ( !proteinChange.isEmpty() ) {
            final int count = proteinChangeCounts.getOrDefault(proteinChange, 0);
            proteinChangeCounts.put(proteinChange, count + 1);
        }
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

                try {
                    return new SimpleInterval(contig, start, end);
                }
                catch (final IllegalArgumentException ex) {
                    // If we have poorly bounded genomic positions, we need to warn the user and move on.
                    // These may occur occasionally in the data.
                    logger.warn("Warning - unable to parse genome position string due to invalid position information.  Ignoring potential COSMIC match with genome position: " + rawPosition);
                    return null;
                }
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
    private SimpleInterval getProteinPositionFromResults(final ResultSet resultSet) {
        return parseProteinString( getProteinChangeStringFromResults(resultSet) );
    }

    /**
     * Pulls a protein change string out of the current record in the given {@link ResultSet}.
     * @param resultSet The results of a query on the database with a current row (must not be {@code null}).
     * @return A {@link String} representing the protein change as found in the current record of the given {@link ResultSet}.  Will not be {@code null}.
     */
    private String getProteinChangeStringFromResults(final ResultSet resultSet) {
        Utils.nonNull(resultSet);

        try {
            final String proteinChangeString = resultSet.getString(PROTEIN_POSITION_COLUMN_NAME);
            return proteinChangeString == null ? "" : proteinChangeString;
        }
        catch (final SQLException ex) {
            throw new GATKException("Cannot get protein change from column: " + GENOME_POSITION_COLUMN_NAME, ex);
        }
    }

    /**
     * Parse a {@link SimpleInterval} from a protein position / protein change.
     * @param proteinPositionString A {@link String} representing a protein position / protein change.
     * @return A {@link SimpleInterval} representing the extents of the given {@code proteinPositionString} or {@code null}.
     */
    private SimpleInterval parseProteinString(final String proteinPositionString) {
        Utils.nonNull(proteinPositionString);

        final Matcher matcher = PROTEIN_POSITION_REGEX.matcher(proteinPositionString);
        if ( matcher.matches() ) {
            try {

                // We have a position, so we should parse it.
                // If we have an end position, we should use it:
                if ( matcher.group(2) != null ) {
                    return new SimpleInterval(PROTEIN_CONTIG, Integer.valueOf(matcher.group(1)), Integer.valueOf(matcher.group(2)));
                }
                else {
                    return new SimpleInterval(PROTEIN_CONTIG, Integer.valueOf(matcher.group(1)), Integer.valueOf(matcher.group(1)));
                }
            }
            catch (final IllegalArgumentException ex) {
                // If we have poorly bounded protein changes, we need to warn the user and move on.
                // These occur occasionally in the data.
                logger.warn("Warning - unable to parse protein string due to invalid position information.  Ignoring potential COSMIC match with protein sequence: " + proteinPositionString);
                return null;
            }
        }
        return null;
    }

    /**
     * Print the given {@link ResultSet} to stdout.
     * @param resultSet The {@link ResultSet} to print.
     */
    private final void printResultSet(final ResultSet resultSet) {
        try {
            final ResultSetMetaData metadata    = resultSet.getMetaData();
            final int               columnCount = metadata.getColumnCount();
            for (int i = 1; i <= columnCount; i++) {
                System.out.print(metadata.getColumnName(i) + ", ");
            }
            System.out.println();
            while ( resultSet.next() ) {
                for ( int i = 1; i <= columnCount; i++ ) {
                    System.out.print(resultSet.getString(i) + ", ");
                }
                System.out.println();
            }
        }
        catch (final SQLException ex) {
            throw new GATKException("Cannot print resultSet!");
        }
    }

    //==================================================================================================================
    // Helper Data Types:

}
