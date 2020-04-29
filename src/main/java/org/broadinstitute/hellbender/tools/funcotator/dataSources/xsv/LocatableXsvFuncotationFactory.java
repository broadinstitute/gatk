package org.broadinstitute.hellbender.tools.funcotator.dataSources.xsv;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.tribble.Feature;
import htsjdk.tribble.readers.AsciiLineReader;
import htsjdk.tribble.readers.AsciiLineReaderIterator;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.funcotator.DataSourceFuncotationFactory;
import org.broadinstitute.hellbender.tools.funcotator.Funcotation;
import org.broadinstitute.hellbender.tools.funcotator.FuncotatorArgumentDefinitions;
import org.broadinstitute.hellbender.tools.funcotator.FuncotatorUtils;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.TableFuncotation;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.codecs.xsvLocatableTable.XsvLocatableTableCodec;
import org.broadinstitute.hellbender.utils.codecs.xsvLocatableTable.XsvTableFeature;

import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

/**
 * Factory for creating {@link TableFuncotation}s by handling `Separated Value` files with arbitrary delimiters
 * (e.g. CSV/TSV files) which contain data that are locatable (i.e. {@link org.broadinstitute.hellbender.utils.codecs.xsvLocatableTable.XsvTableFeature}).
 *
 * This is a high-level object that interfaces with the internals of {@link org.broadinstitute.hellbender.tools.funcotator.Funcotator}.
 * Created by jonn on 12/6/17.
 */
public class LocatableXsvFuncotationFactory extends DataSourceFuncotationFactory {

    //==================================================================================================================
    // Public Static Members:

    //==================================================================================================================
    // Private Static Members:

    @VisibleForTesting
    static String DEFAULT_NAME = "LocatableXsv";

    //==================================================================================================================
    // Private Members:

    /**
     * Name of this {@link LocatableXsvFuncotationFactory}.
     */
    private final String name;

    /**
     * {@link LinkedHashSet} of the names of all fields supported by this {@link LocatableXsvFuncotationFactory}.
     * Set by {@link #setSupportedFuncotationFields(Path)}.
     */
    private LinkedHashSet<String> supportedFieldNames = null;

    /**
     * {@link List} of the names of all fields supported by this {@link LocatableXsvFuncotationFactory}.
     * Set by {@link #setSupportedFuncotationFields(Path)}.
     */
    private List<String> supportedFieldNameList = null;

    /**
     * {@link List} of empty {@link String}s of the same length as {@link #supportedFieldNames}.
     * Cached for faster output.
     * Set by {@link #setSupportedFuncotationFields(Path)}.
     */
    private List<String> emptyFieldList = null;

    //==================================================================================================================
    // Constructors:

    /**
     * Create a {@link LocatableXsvFuncotationFactory}.
     * @param name A {@link String} containing the name of this {@link LocatableXsvFuncotationFactory}.
     * @param version  The version {@link String} of the backing data source from which {@link Funcotation}s will be made.
     * @param annotationOverridesMap A {@link LinkedHashMap<String,String>} containing user-specified overrides for specific {@link Funcotation}s.
     * @param mainSourceFileAsFeatureInput The backing {@link FeatureInput} for this {@link LocatableXsvFuncotationFactory}, from which all {@link Funcotation}s will be created.
     */
    public LocatableXsvFuncotationFactory(final String name, final String version, final LinkedHashMap<String, String> annotationOverridesMap,
                                          final FeatureInput<? extends Feature> mainSourceFileAsFeatureInput){
        this(name, version, annotationOverridesMap, mainSourceFileAsFeatureInput, false);
    }

    /**
     * Create a {@link LocatableXsvFuncotationFactory}.
     * @param name A {@link String} containing the name of this {@link LocatableXsvFuncotationFactory}.
     * @param version  The version {@link String} of the backing data source from which {@link Funcotation}s will be made.
     * @param annotationOverridesMap A {@link LinkedHashMap<String,String>} containing user-specified overrides for specific {@link Funcotation}s.
     * @param mainSourceFileAsFeatureInput The backing {@link FeatureInput} for this {@link LocatableXsvFuncotationFactory}, from which all {@link Funcotation}s will be created.
     * @param isDataSourceB37 If {@code true}, indicates that the data source behind this {@link LocatableXsvFuncotationFactory} contains B37 data.
     */
    public LocatableXsvFuncotationFactory(final String name, final String version, final LinkedHashMap<String, String> annotationOverridesMap,
                                          final FeatureInput<? extends Feature> mainSourceFileAsFeatureInput,
                                          final boolean isDataSourceB37){
        this(name, version, annotationOverridesMap, mainSourceFileAsFeatureInput, isDataSourceB37, FuncotatorUtils.DEFAULT_MIN_NUM_BASES_FOR_VALID_SEGMENT);
    }

    /**
     * Create a {@link LocatableXsvFuncotationFactory}.
     * @param name A {@link String} containing the name of this {@link LocatableXsvFuncotationFactory}.
     * @param version  The version {@link String} of the backing data source from which {@link Funcotation}s will be made.
     * @param annotationOverridesMap A {@link LinkedHashMap<String,String>} containing user-specified overrides for specific {@link Funcotation}s.
     * @param mainSourceFileAsFeatureInput The backing {@link FeatureInput} for this {@link LocatableXsvFuncotationFactory}, from which all {@link Funcotation}s will be created.
     * @param isDataSourceB37 If {@code true}, indicates that the data source behind this {@link LocatableXsvFuncotationFactory} contains B37 data.
     * @param minBasesForValidSegment The minimum number of bases for a segment to be considered valid.
     */
    public LocatableXsvFuncotationFactory(final String name, final String version, final LinkedHashMap<String, String> annotationOverridesMap,
                                          final FeatureInput<? extends Feature> mainSourceFileAsFeatureInput,
                                          final boolean isDataSourceB37, final int minBasesForValidSegment){

        super(mainSourceFileAsFeatureInput, minBasesForValidSegment);

        this.name = name;
        this.version = version;

        this.annotationOverrideMap = new LinkedHashMap<>(annotationOverridesMap);
        this.dataSourceIsB37 = isDataSourceB37;
    }

    //==================================================================================================================
    // Override Methods:

    @Override
    public Class<? extends Feature> getAnnotationFeatureClass() {
        return XsvTableFeature.class;
    }

    @Override
    public String getName() {
        return name;
    }

    @Override
    public LinkedHashSet<String> getSupportedFuncotationFields() {

        if ( supportedFieldNames == null ) {
            throw new GATKException("Must set supportedFuncotationFields before querying for them!");
        }
        else {
            return new LinkedHashSet<>(supportedFieldNames);
        }
    }

    @Override
    protected List<Funcotation> createDefaultFuncotationsOnVariant(final VariantContext variant, final ReferenceContext referenceContext ) {
        return createDefaultFuncotationsOnVariantHelper(variant, referenceContext, Collections.emptySet());
    }

    @Override
    protected List<Funcotation> createFuncotationsOnVariant(final VariantContext variant, final ReferenceContext referenceContext, final List<Feature> featureList) {
        final List<Funcotation> outputFuncotations = new ArrayList<>();

        final List<Allele> alternateAlleles = variant.getAlternateAlleles();

        // Create a set to put our annotated Alternate alleles in.
        // We'll use this to determine if the alt allele has been annotated.
        final Set<Allele> annotatedAltAlleles = new HashSet<>(alternateAlleles.size());

        if ( !featureList.isEmpty() ) {
            for ( final Feature feature : featureList ) {

                // Get the kind of feature we want here:
                if ( feature != null ) {

                    // By this point we know the feature type is correct, so we cast it:
                    final XsvTableFeature tableFeature = (XsvTableFeature) feature;

                    // Make sure we are annotating the correct XSVTableFeature:
                    // NOTE: This is probably redundant to the name check in DataSourceFuncotationFactory
                    if ( tableFeature.getDataSourceName().equals(name) ) {

                        // Now we create one funcotation for each Alternate allele:
                        for ( final Allele altAllele : alternateAlleles ) {
                            outputFuncotations.add(TableFuncotation.create(tableFeature, altAllele, name, null));
                            annotatedAltAlleles.add(altAllele);
                        }
                    }

                    // TODO: Must break the loop now to prevent multiple entries messing up the number of fields in the funcotation (issue #4930 - https://github.com/broadinstitute/gatk/issues/4930)
                    break;
                }
            }
        }

        // If we didn't add funcotations for an allele, we should add in blank funcotations to that allele for each field that can be produced
        // by this LocatableXsvFuncotationFactory:
        if ( annotatedAltAlleles.size() != alternateAlleles.size() ) {
            outputFuncotations.addAll(createDefaultFuncotationsOnVariantHelper(variant, referenceContext, annotatedAltAlleles));
        }

        return outputFuncotations;
    }

    @Override
    protected List<Funcotation> createFuncotationsOnVariant(final VariantContext variant, final ReferenceContext referenceContext, final List<Feature> featureList, final List<GencodeFuncotation> gencodeFuncotations) {
        return createFuncotationsOnVariant(variant, referenceContext, featureList);
    }

    @Override
    public FuncotatorArgumentDefinitions.DataSourceType getType() {
        return FuncotatorArgumentDefinitions.DataSourceType.LOCATABLE_XSV;
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
                funcotationList.add(TableFuncotation.create(supportedFieldNameList, emptyFieldList, altAllele, name, null));
            }
        }

        return funcotationList;
    }

    /**
     * Set the field names that this {@link LocatableXsvFuncotationFactory} can create.
     * Does so by reading the headers of backing data files for this {@link LocatableXsvFuncotationFactory}.
     * @param inputDataFilePath {@link Path} to a backing data file from which annotations can be made for this {@link LocatableXsvFuncotationFactory}.  Must not be {@code null}.
     */
    public void setSupportedFuncotationFields(final Path inputDataFilePath) {

        Utils.nonNull(inputDataFilePath);

        if ( supportedFieldNames == null ) {
            synchronized ( this ) {
                if ( supportedFieldNames == null ) {

                    // Approximate / arbitrary starting size:
                    supportedFieldNames = new LinkedHashSet<>(10);

                    // Set up a codec here to read the config file.
                    // We have to call canDecode to set up the internal state of the XsvLocatableTableCodec:
                    final XsvLocatableTableCodec codec = new XsvLocatableTableCodec();
                    try {
                        if ( !codec.canDecode(mainSourceFileAsFeatureInput.getFeaturePath()) ) {
                            // This should never happen because we have already validated this config file by the time we
                            // reach here:
                            throw new GATKException.ShouldNeverReachHereException("Could not decode from data file: " + mainSourceFileAsFeatureInput.getFeaturePath());
                        }
                    }
                    catch ( final NullPointerException ex ) {
                        // This should never happen because we have already validated this config file by the time we
                        // reach here:
                        throw new GATKException.ShouldNeverReachHereException("Could not decode from data file!  Has not been set yet!");
                    }

                    // Get the info from our path:
                    final List<String> columnNames;
                    try (final InputStream fileInputStream = Files.newInputStream(inputDataFilePath)) {

                        final AsciiLineReaderIterator lineReaderIterator = new AsciiLineReaderIterator(AsciiLineReader.from(fileInputStream));
                        codec.readActualHeader(lineReaderIterator);
                        columnNames = codec.getHeaderWithoutLocationColumns();

                    } catch (final IOException ioe) {
                        throw new UserException.BadInput("Could not read header from data file: " + inputDataFilePath.toUri().toString(), ioe);
                    }

                    // Make sure we actually read the header:
                    if ( columnNames == null ) {
                        throw new UserException.MalformedFile("Could not decode from data file: " + inputDataFilePath.toUri().toString());
                    }

                    supportedFieldNames.addAll(columnNames);


                    // Initialize our field name lists:
                    initializeFieldNameLists();

                    // Adjust the manual annotations to make sure we don't try to annotate any fields we aren't
                    // responsible for:
                    annotationOverrideMap.entrySet().removeIf( e -> !supportedFieldNames.contains(e.getKey()) );
                }
            }
        }
    }

    /**
     * Initialize {@link #supportedFieldNameList} and {@link #emptyFieldList} given a populated {@link #supportedFieldNames}.
     */
    private void initializeFieldNameLists() {

        if ( supportedFieldNames == null ) {
            throw new GATKException("Must set supportedFuncotationFields before initializing field name lists!");
        }

        supportedFieldNameList = new ArrayList<>(supportedFieldNames);
        emptyFieldList = new ArrayList<>(supportedFieldNameList.size());
        for ( final String s : supportedFieldNameList) {
            emptyFieldList.add("");
        }
    }

    //==================================================================================================================
    // Helper Data Types:

}
