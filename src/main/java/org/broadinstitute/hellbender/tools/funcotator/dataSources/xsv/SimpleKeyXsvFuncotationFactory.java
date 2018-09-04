package org.broadinstitute.hellbender.tools.funcotator.dataSources.xsv;

import htsjdk.tribble.Feature;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.funcotator.DataSourceFuncotationFactory;
import org.broadinstitute.hellbender.tools.funcotator.Funcotation;
import org.broadinstitute.hellbender.tools.funcotator.FuncotatorArgumentDefinitions;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.TableFuncotation;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.nio.PathLineIterator;

import java.nio.file.Path;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Factory for creating {@link TableFuncotation}s by handling `Separated Value` files with arbitrary delimiters
 * (e.g. CSV/TSV files) which contain data that use a simple key (i.e. {@link XsvDataKeyType}).
 *
 * This is a high-level object that interfaces with the internals of {@link org.broadinstitute.hellbender.tools.funcotator.Funcotator}.
 * Created by jonn on 11/28/17.
 */
public class SimpleKeyXsvFuncotationFactory extends DataSourceFuncotationFactory {

    //==================================================================================================================
    // Public Static Members:

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
    private final Path xsvInputPath;

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
     * The number of lines at the top of the file to ignore before starting to parse the input file.
     */
    private final int numHeaderLinesToIgnore;

    /**
     * The names of the columns containing values that will be added to the resulting {@link TableFuncotation}.
     */
    private final List<String> annotationColumnNames;

    /**
     * A {@link List} of empty {@link String}s to use when we must create an annotation but the features do not match
     * the query.
     */
    private final List<String> emptyAnnotationList;

    /**
     * Map containing the annotations that we have to
     */
    private final Map<String, List<String>> annotationMap;

    //==================================================================================================================
    // Constructors:

    public SimpleKeyXsvFuncotationFactory(final String name, final Path filePath, final String version, final String delim, final int keyColumn, final XsvDataKeyType keyType) {
        this(name, filePath, version, delim, keyColumn, keyType, new LinkedHashMap<>(), 0, false);
    }

    public SimpleKeyXsvFuncotationFactory(final String name, final Path filePath, final String version, final String delim, final int keyColumn, final XsvDataKeyType keyType, final LinkedHashMap<String, String> annotationOverrides) {
        this(name, filePath, version, delim, keyColumn, keyType, annotationOverrides, 0, false);
    }

    public SimpleKeyXsvFuncotationFactory(final String name, final Path filePath, final String version, final String delim, final int keyColumn, final XsvDataKeyType keyType, final LinkedHashMap<String, String> annotationOverrides, final int numHeaderLinesToIgnore) {
        this(name, filePath, version, delim, keyColumn, keyType, annotationOverrides, numHeaderLinesToIgnore, false);
    }

    public SimpleKeyXsvFuncotationFactory(final String name,
                                          final Path filePath,
                                          final String version,
                                          final String delim,
                                          final int keyColumn,
                                          final XsvDataKeyType keyType,
                                          final LinkedHashMap<String, String> annotationOverrides,
                                          final int numHeaderLinesToIgnore,
                                          final boolean permissiveColumns ) {
        this.name = name;

        delimiter = delim;
        xsvInputPath = filePath;

        this.version = version;
        this.keyColumn = keyColumn;
        this.keyType = keyType;

        annotationOverrideMap = annotationOverrides;

        this.numHeaderLinesToIgnore = numHeaderLinesToIgnore;

        // Initialize our annotations map:
        annotationMap = new HashMap<>();

        // Create our iterator:
        try ( final PathLineIterator pathLineIterator = new PathLineIterator(xsvInputPath) ) {

            // Get a line iterator for our lines:
            final Iterator<String> it = pathLineIterator.iterator();

            // Get our column names:
            annotationColumnNames = createColumnNames(it, numHeaderLinesToIgnore);

            // Populate our empty annotation list:
            emptyAnnotationList = new ArrayList<>(annotationColumnNames.size());
            for ( final String s : annotationColumnNames ) {
                emptyAnnotationList.add("");
            }

            // Populate our annotation map:
            populateAnnotationMap(it, permissiveColumns);
        }

        // Initialize overrides / defaults:
        initializeAnnotationOverrides( annotationOverrides );
    }

    //==================================================================================================================
    // Override Methods:

    @Override
    protected Class<? extends Feature> getAnnotationFeatureClass() {
        // Returning Feature.class here implies that this class doesn't care about what features it gets.
        return Feature.class;
    }

    @Override
    public String getName() {
        return name;
    }

    @Override
    public LinkedHashSet<String> getSupportedFuncotationFields() {
        return new LinkedHashSet<>(annotationColumnNames);
    }

    @Override
    protected List<Funcotation> createDefaultFuncotationsOnVariant( final VariantContext variant, final ReferenceContext referenceContext ) {
        return createDefaultFuncotationsOnVariantHelper(variant, referenceContext, Collections.emptySet());
    }

    @Override
    /**
     * {@inheritDoc}
     * This method should never be called on a {@link SimpleKeyXsvFuncotationFactory} - knowledge of the applied
     * {@link GencodeFuncotation}s is required to create an {@link TableFuncotation} from here.
     */
    protected List<Funcotation> createFuncotationsOnVariant(final VariantContext variant, final ReferenceContext referenceContext, final List<Feature> featureList) {
        throw new GATKException(this.getClass().getName() + " requires a set of GencodeFuncotations in order to createFuncotationsOnVariant!  This method should never be called on a " + this.getClass().getName());
    }

    @Override
    /**
     * {@inheritDoc}
     * For each {@link org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation}, the Transcript ID or Gene Name (Hugo Symbol)
     * is checked for a match against the key of any annotation in {@link SimpleKeyXsvFuncotationFactory#annotationMap}.
     * If a match is found, an {@link TableFuncotation} is added to the list to be returned.
     */
    protected List<Funcotation> createFuncotationsOnVariant(final VariantContext variant,
                                                         final ReferenceContext referenceContext,
                                                         final List<Feature> featureList,
                                                         final List<GencodeFuncotation> gencodeFuncotations) {
        final List<Funcotation> outputFuncotations = new ArrayList<>();

        // If we have gencodeFuncotations we go through them and check for the correct Gene Name / TranscriptID.
        // If any match, we create our xsvFuncotation for this variant:
        for ( final GencodeFuncotation gencodeFuncotation : gencodeFuncotations ) {

            // Create a set to put our annotated Alternate alleles in.
            // We'll use this to determine if the alt allele has been annotated.
            final Set<Allele> annotatedAltAlleles = new HashSet<>(variant.getAlternateAlleles().size());

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
                // Create 1 annotation for each alt allele and add our annotations to the list:
                for ( final Allele altAllele : variant.getAlternateAlleles() ) {
                    outputFuncotations.add(TableFuncotation.create(annotationColumnNames, annotations, altAllele, name, null));
                    annotatedAltAlleles.add(altAllele);
                }
            }

            // If we didn't add funcotations for an allele, we should add in blank funcotations to that allele for
            // each field that can be produced by this SimpleKeyXsvFuncotationFactory:
            if ( annotatedAltAlleles.size() != variant.getAlternateAlleles().size() ) {
                outputFuncotations.addAll(createDefaultFuncotationsOnVariantHelper(variant, referenceContext, annotatedAltAlleles));
            }
        }

        return outputFuncotations;
    }

    @Override
    public FuncotatorArgumentDefinitions.DataSourceType getType() {
        return FuncotatorArgumentDefinitions.DataSourceType.SIMPLE_XSV;
    }

    //==================================================================================================================
    // Static Methods:

    //==================================================================================================================
    // Instance Methods:

    private List<Funcotation> createDefaultFuncotationsOnVariantHelper( final VariantContext variant, final ReferenceContext referenceContext, final Set<Allele> annotatedAltAlleles  ) {

        final List<Funcotation> funcotationList = new ArrayList<>();

        final List<Allele> alternateAlleles = variant.getAlternateAlleles();

        for ( final Allele altAllele : alternateAlleles ) {
            if ( !annotatedAltAlleles.contains(altAllele) ) {
                funcotationList.add(TableFuncotation.create(annotationColumnNames, emptyAnnotationList, altAllele, name, null));
            }
        }

        return funcotationList;
    }

    /**
     * Creates the annotation column names from the given iterator.
     * @param lineIterator An iterator at the start of an XSV file from which to get the header columns.
     * @param numHeaderLinesToIgnore The number of lines at the start of the file to ignore before beginning parsing.
     */
    private List<String> createColumnNames(final Iterator<String> lineIterator,
                                            final int numHeaderLinesToIgnore) {
        // Ignore the leading lines that we were told to ignore:
        for ( int i = 0; i < numHeaderLinesToIgnore ; ++i ) {
            lineIterator.next();
        }

        // We're at the header, so we need to initialize the header columns,
        // And fix the column headers to not contain any spaces:
        final List<String> annotationColumnNames =
                Utils.split(lineIterator.next(), delimiter).stream()
                        .map(n -> n.replaceAll("^\\s+", "").replaceAll("\\s+$", ""))
                        .map(n -> getName() + "_" + n.replaceAll(" ", "_"))
                        .collect(Collectors.toCollection(ArrayList::new));

        // If the number of columns is < 2, we don't have any data (because we don't add in the column containing
        // the key).  This is an error:
        if ( annotationColumnNames.size() < 2 ) {
            throw new UserException.MalformedFile("Data Source is badly formatted (" + xsvInputPath.toUri().toString() + ") - contains too few columns (" + annotationColumnNames.size() + ")!  Is the delimiter specified incorrectly?");
        }

        // Pull out the column containing the key so it doesn't appear in our data:
        annotationColumnNames.remove(keyColumn);

        return annotationColumnNames;
    }

    /**
     * Populates {@link SimpleKeyXsvFuncotationFactory#annotationMap} with data from the given iterator.
     * Assumes that {@link SimpleKeyXsvFuncotationFactory#annotationColumnNames} is populated.
     * @param it An {@link Iterator} of {@link String} starting at the first data line in the file to parse.
     * @param permissiveColumnNumbers A flag which if true indicates to allow mismatches between the number of columns in the header row and data rows.
     *                                In the event of a mismatch, the data row with either be padded to the number of columns in the header
     *                                or it will be truncated to match the number of columns in the header.
     */
    private void populateAnnotationMap(final Iterator<String> it,
                                       final boolean permissiveColumnNumbers) {

        boolean emptyLineFlag = false;

        int dataRowNum = numHeaderLinesToIgnore + 2;

        // Parse the rest of the data:
        while ( it.hasNext() ) {

            final String rawRow = it.next();

            // Check for an empty line.
            // The only permissible place is at the end of the file.
            if ( rawRow.length() == 0 ) {
                if ( !emptyLineFlag ) {
                    emptyLineFlag = true;
                    continue;
                }
            }
            if ( emptyLineFlag ) {
                    throw new UserException.MalformedFile("File contains an empty line (" + dataRowNum + ").  All lines must have data.");
            }

            // The use of Utils.split is NECESSARY here because delimiter could potentially evaluate
            // as a regular expression (as in the case of '|').
            final List<String> dataRow = Utils.split(rawRow, delimiter);

            // Remove the key column:
            final String rowKey = dataRow.remove(keyColumn);

            // Make sure we have the same number of columns:
            if ( (dataRow.size() != annotationColumnNames.size()) ) {
                if ( !permissiveColumnNumbers ) {
                    throw new UserException.MalformedFile("Data Source is badly formatted (" + xsvInputPath.toUri().toString() + ") - row " + dataRowNum + " does not contain the same number of columns as header (" + dataRow.size() + " != " + annotationColumnNames.size() + ")!");
                }
                else if ( dataRow.size() > annotationColumnNames.size() ) {
                    dataRow.remove(dataRow.size()-1);
                }
                else {
                    while ( dataRow.size() < annotationColumnNames.size() ) {
                        dataRow.add( "" );
                    }
                }
            }

            // Store this in our map:
            annotationMap.put(rowKey, dataRow);

            // Increment our row counter:
            ++dataRowNum;
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
}
