package org.broadinstitute.hellbender.tools.funcotator.dataSources.XSV;

import htsjdk.tribble.Feature;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.funcotator.DataSourceFuncotationFactory;
import org.broadinstitute.hellbender.tools.funcotator.Funcotation;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation;
import org.broadinstitute.hellbender.utils.nio.PathLineIterator;

import java.nio.file.Path;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Factory for creating {@link XSVFuncotation}s by handling `Separated Value` files with arbitrary delimiters
 * (e.g. CSV/TSV files) which contain data that use a simple key (i.e. {@link XsvDataKeyType}).
 *
 * This is a high-level object that interfaces with the internals of {@link org.broadinstitute.hellbender.tools.funcotator.Funcotator}.
 * Created by jonn on 11/28/17.
 */
public class SimpleKeyXsvFuncotationFactory extends DataSourceFuncotationFactory {

    //==================================================================================================================
    // Private Static Members:

    //==================================================================================================================
    // Private Members:

    /**
     * The name of the data source associated with this {@link SimpleKeyXsvFuncotationFactory}.
     */
    private final String name;

    /**
     * Delimiter used in this XSV file.
     */
    private final String delimiter;

    /**
     * Path to the given XSV file.
     */
    private final Path xsvFilePath;

    /**
     * The column (0-indexed) containing the key for this XSV file.
     * The key is used to annotate a given {@link VariantContext} with the data XSV file with the other funcotations.
     */
    private final int keyColumn;

    /**
     * The type of key used in this {@link SimpleKeyXsvFuncotationFactory}.
     */
    private final XsvDataKeyType keyType;

    /**
     * The names of the columns containing values that will be added to the resulting {@link XSVFuncotation}.
     */
    private final List<String> annotationColumnNames;

    /**
     * Map containing the annotations that we have to
     */
    private final Map<String, List<String>> annotationMap;

    //==================================================================================================================
    // Constructors:

    public SimpleKeyXsvFuncotationFactory(final String name, final Path filePath, final String delim, final int keyColumn, final XsvDataKeyType keyType) {
        this(name, filePath, delim, keyColumn, keyType, 0);
    }

    public SimpleKeyXsvFuncotationFactory(final String name,
                                          final Path filePath,
                                          final String delim,
                                          final int keyColumn,
                                          final XsvDataKeyType keyType,
                                          final int numHeaderLinesToIgnore) {
        this.name = name;

        delimiter = delim;
        xsvFilePath = filePath;

        this.keyColumn = keyColumn;
        this.keyType = keyType;

        // Initialize our annotations map:
        annotationMap = new HashMap<>();

        // Create our iterator:
        try ( final PathLineIterator pathLineIterator = new PathLineIterator(xsvFilePath) ) {

            // Get a line iterator for our lines:
            final Iterator<String> it = pathLineIterator.iterator();

            // Get our column names:
            annotationColumnNames = createColumnNames( it, numHeaderLinesToIgnore );

            // Populate our annotation map:
            populateAnnotationMap(it);
        }
    }

    //==================================================================================================================
    // Override Methods:

    @Override
    public String getName() {
        return name;
    }

    @Override
    public LinkedHashSet<String> getSupportedFuncotationFields() {
        return new LinkedHashSet<>(annotationColumnNames.stream().map(f -> getName() + "_" + f).collect(Collectors.toList()));
    }

    @Override
    /**
     * {@inheritDoc}
     * This method should never be called on a {@link SimpleKeyXsvFuncotationFactory} - knowledge of the applied
     * {@link GencodeFuncotation}s is required to create an {@link XSVFuncotation} from here.
     *
     */
    public List<Funcotation> createFuncotations(final VariantContext variant, final ReferenceContext referenceContext, final List<Feature> featureList) {
        throw new GATKException("This method should never be called on a " + this.getClass().getName());
    }

    @Override
    /**
     * {@inheritDoc}
     * For each {@link org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation}, the Transcript ID or Gene Name (Hugo Symbol)
     * is checked for a match against the key of any annotation in {@link SimpleKeyXsvFuncotationFactory#annotationMap}.
     * If a match is found, an {@link XSVFuncotation} is added to the list to be returned.
     */
    public List<Funcotation> createFuncotations(final VariantContext variant,
                                                final ReferenceContext referenceContext,
                                                final List<Feature> featureList,
                                                final List<GencodeFuncotation> gencodeFuncotations) {
        final List<Funcotation> outputFuncotations = new ArrayList<>();

        // If we have gencodeFuncotations we go through them and check for the correct Gene Name / TranscriptID.
        // If any match, we create our xsvFuncotation for this variant:
        for (  final GencodeFuncotation gencodeFuncotation : gencodeFuncotations ) {

            // Get our key from the gencode funcotation:
            final String key;
            if ( keyType == XsvDataKeyType.GENE_NAME ) {
                key = gencodeFuncotation.getHugoSymbol();
            }
            else {
                key = gencodeFuncotation.getAnnotationTranscript();
            }

            // Get our annotations:
            final List<String> annotations = annotationMap.get( key );
            if ( annotations != null ) {
                // Add our annotations to the list:
                outputFuncotations.add( new XSVFuncotation(annotationColumnNames, annotations) );
            }
        }

        setOverrideValuesInFuncotations(outputFuncotations);

        return outputFuncotations;
    }

    //==================================================================================================================
    // Static Methods:

    //==================================================================================================================
    // Instance Methods:

    /**
     * Creates the annotation column names from the given iterator.
     * @param lineIterator An iterator at the start of an XSV file from which to get the header columns.
     * @param numHeaderLinesToIgnore The number of lines at the start of the file to ignore before beginning parsing.
     */
    private List<String> createColumnNames(final Iterator<String> lineIterator, final int numHeaderLinesToIgnore) {
        // Ignore the leading lines that we were told to ignore:
        for ( int i = 0; i < numHeaderLinesToIgnore ; ++i ) {
            lineIterator.next();
        }

        // We're at the header, so we need to initialize the header columns:
        final List<String> annotationColumnNames = new ArrayList<>( Arrays.asList(lineIterator.next().split(delimiter)) );

        // If the number of columns is < 2, we don't have any data (because we don't add in the column containing
        // the key).  This is an error:
        if ( annotationColumnNames.size() < 2 ) {
            throw new UserException.MalformedFile("Data Source is badly formatted (" + xsvFilePath.toUri().toString() + ") - contains too few columns (" + annotationColumnNames.size() + ")!");
        }

        // Pull out the column containing the key so it doesn't appear in our data:
        annotationColumnNames.remove(keyColumn);

        return annotationColumnNames;
    }

    /**
     * Populates {@link SimpleKeyXsvFuncotationFactory#annotationMap} with data from the given iterator.
     * Assumes that {@link SimpleKeyXsvFuncotationFactory#annotationColumnNames} is populated.
     * @param it An {@link Iterator} of {@link String} starting at the first data line in the file to parse.
     */
    private void populateAnnotationMap(final Iterator<String> it) {
        // Parse the rest of the data:
        while ( it.hasNext() ) {
            final List<String> dataRow = new ArrayList<>( Arrays.asList(it.next().split(delimiter)) );

            // Make sure we have the same number of columns:
            if ( dataRow.size() != annotationColumnNames.size() ) {
                throw new UserException.MalformedFile("Data Source is badly formatted (" + xsvFilePath.toUri().toString() + ") - row does not contain the same number of columns as header (" + dataRow.size() + " != " + annotationColumnNames.size() + ")!");
            }

            // Remove the key column:
            final String rowKey = dataRow.remove(keyColumn);

            // Store this in our map:
            annotationMap.put(rowKey, dataRow);
        }
    }

    //==================================================================================================================
    // Helper Data Types:

    public enum XsvDataKeyType {
        /**
         * The key specified is a Gene Name which will be used to match and annotate a {@link VariantContext}.
         * All given {@link VariantContext}s matching the Gene Name will be annotated with that row's data.
         */
        GENE_NAME,

        /**
         * The key specified is a Transcript ID which will be used to match and annotate a {@link VariantContext}.
         * All given {@link VariantContext}s matching the Transcript ID will be annotated with that row's data.
         */
        TRANSCRIPT_ID
    }

    // TODO: Use `Files.readAttributes` to get size/last modified time / etc. to approximate equality test.
}
