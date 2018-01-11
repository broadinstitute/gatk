package org.broadinstitute.hellbender.tools.funcotator.dataSources.xsv;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.tribble.Feature;
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
import org.broadinstitute.hellbender.utils.codecs.xsvLocatableTable.XsvLocatableTableCodec;
import org.broadinstitute.hellbender.utils.codecs.xsvLocatableTable.XsvTableFeature;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.util.*;
import java.util.stream.Collectors;

/**
 *  Factory for creating {@link TableFuncotation}s by handling `Separated Value` files with arbitrary delimiters
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
     * Names of all fields supported by this {@link LocatableXsvFuncotationFactory}.
     */
    private LinkedHashSet<String> supportedFieldNames = null;

    //==================================================================================================================
    // Constructors:

    public LocatableXsvFuncotationFactory(){
        this(DEFAULT_NAME);
    }

    public LocatableXsvFuncotationFactory(final String name){
        this.name = name;
    }

    //==================================================================================================================
    // Override Methods:

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
    public List<Funcotation> createFuncotations(final VariantContext variant, final ReferenceContext referenceContext, final List<Feature> featureList) {
        final List<Funcotation> outputFuncotations = new ArrayList<>();

        if ( !featureList.isEmpty() ) {
            for ( final Feature feature : featureList ) {
                // Get the kind of feature we want here:
                if ( (feature != null) && XsvTableFeature.class.isAssignableFrom(feature.getClass()) ) {
                    outputFuncotations.add( new TableFuncotation((XsvTableFeature) feature) );
                }
            }
        }

        return outputFuncotations;
    }

    @Override
    public List<Funcotation> createFuncotations(final VariantContext variant, final ReferenceContext referenceContext, final List<Feature> featureList, final List<GencodeFuncotation> gencodeFuncotations) {
        return createFuncotations(variant, referenceContext, featureList);
    }

    @Override
    public FuncotatorArgumentDefinitions.DataSourceType getType() {
        return FuncotatorArgumentDefinitions.DataSourceType.LOCATABLE_XSV;
    }

    //==================================================================================================================
    // Static Methods:

    //==================================================================================================================
    // Instance Methods:

    public void setSupportedFuncotationFields(final List<Path> inputDataFilePaths) {

        // Approximate starting size:
        supportedFieldNames = new LinkedHashSet<>(inputDataFilePaths.size() * 10);

        for ( final Path dataPath : inputDataFilePaths ) {

            // Get the associated config file:
            final Path configPath = XsvLocatableTableCodec.getConfigFilePath(dataPath);

            // Get the associated config properties:
            final Properties configProperties = XsvLocatableTableCodec.getAndValidateConfigFileContents(configPath);

            final int contigColumn   = Integer.valueOf(configProperties.getProperty(XsvLocatableTableCodec.CONFIG_FILE_CONTIG_COLUMN_KEY));
            final int startColumn    = Integer.valueOf(configProperties.getProperty(XsvLocatableTableCodec.CONFIG_FILE_START_COLUMN_KEY));
            final int endColumn      = Integer.valueOf(configProperties.getProperty(XsvLocatableTableCodec.CONFIG_FILE_END_COLUMN_KEY));
            final String delimiter      = configProperties.getProperty(XsvLocatableTableCodec.CONFIG_FILE_DELIMITER_KEY);
            final String dataSourceName = configProperties.getProperty(XsvLocatableTableCodec.CONFIG_FILE_DATA_SOURCE_NAME_KEY);

            // Create our index to remove list:
            final Set<Integer> indicesToRemove = new HashSet<>(Arrays.asList(contigColumn, startColumn, endColumn));

            // Get the raw header:
            List<String> header = null;
            try( final BufferedReader inputReader =
                         new BufferedReader(new InputStreamReader(Files.newInputStream(dataPath, StandardOpenOption.READ))) ) {

                String line = inputReader.readLine();
                while (line != null) {
                    if ( !line.startsWith(XsvLocatableTableCodec.COMMENT_DELIMITER) ) {
                        header = Utils.split(line, delimiter).stream()
                                .map(x -> dataSourceName + "_" + x)
                                .collect(Collectors.toCollection(ArrayList::new));
                        break;
                    }
                    line = inputReader.readLine();
                }
            }
            catch (final Exception ex) {
                throw new UserException.BadInput("Error while reading from input data file: " + dataPath.toUri().toString(), ex);
            }

            // Make sure we actually read the header:
            if ( header == null ) {
                throw new UserException.MalformedFile("Could not read header from data file: " + dataPath.toUri().toString());
            }

            // Add the header fields to the supportedFieldNames:
            for ( int i = 0 ; i < header.size() ; ++i ) {
                if ( !indicesToRemove.contains(i) ) {
                    supportedFieldNames.add(header.get(i));
                }
            }
        }
    }

    //==================================================================================================================
    // Helper Data Types:

}
