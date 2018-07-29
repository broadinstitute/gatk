package org.broadinstitute.hellbender.tools.funcotator.dataSources;

import com.google.common.collect.Sets;
import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.funcotator.Funcotation;
import org.broadinstitute.hellbender.tools.funcotator.FuncotatorUtils;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.xsv.LocatableXsvFuncotationFactory;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.xsv.SimpleKeyXsvFuncotationFactory;
import org.broadinstitute.hellbender.tools.funcotator.metadata.FuncotationMetadata;
import org.broadinstitute.hellbender.tools.funcotator.metadata.FuncotationMetadataUtils;
import org.broadinstitute.hellbender.tools.funcotator.vcfOutput.VcfOutputRenderer;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.codecs.xsvLocatableTable.XsvTableFeature;

import java.util.*;
import java.util.stream.Collectors;

/**
 * A {@link Funcotation} to hold data from simple tabular data.
 * This class will be used any time an annotation is created by reading data from any data type
 * that can be expressed as a row in a table / database (such as via {@link SimpleKeyXsvFuncotationFactory},
 * {@link LocatableXsvFuncotationFactory}, and
 * {@link org.broadinstitute.hellbender.tools.funcotator.dataSources.cosmic.CosmicFuncotationFactory}).
 * Created by jonn on 11/28/17.
 */
public class TableFuncotation implements Funcotation {

    //==================================================================================================================
    // Private Static Members:

    //==================================================================================================================
    // Private Members:

    /**
     * Name of the factory that created this {@link TableFuncotation}.
     */
    final private String dataSourceName;

    /**
     * Names of the fields in this {@link TableFuncotation}.
     */
    final private LinkedHashMap<String, String> fieldMap;

    /** The alternate {@link Allele} associated with this {@link TableFuncotation} */
    final private Allele altAllele;

    final private FuncotationMetadata metadata;

    //==================================================================================================================
    // Constructors:

    private TableFuncotation(final List<String> fieldNames, final List<String> fieldValues, final Allele altAllele, final String dataSourceName, final FuncotationMetadata metadata ) {
        if ( fieldNames.size() != fieldValues.size() ) {
            throw new UserException.BadInput("Field names and Field values are of different lengths!  This must not be!");
        }

        fieldMap = new LinkedHashMap<>(fieldNames.size());
        for ( int i = 0; i < fieldNames.size() ; ++i ) {
            fieldMap.put(fieldNames.get(i), fieldValues.get(i));
        }

        this.altAllele = altAllele;
        this.dataSourceName = dataSourceName;

        if (metadata == null) {
            this.metadata = FuncotationMetadataUtils.createWithUnknownAttributes(fieldNames);
        } else {
            this.metadata = metadata;
        }

        // Validate that the metadata is okay.
        final Set<String> metadataFieldNames = this.metadata.retrieveAllHeaderInfo().stream().map(f -> f.getID()).collect(Collectors.toSet());
        final HashSet<String> funcotationFieldNames = new HashSet<>(fieldNames);
        if (!metadataFieldNames.equals(funcotationFieldNames)) {
            throw new UserException.BadInput("Metadata was not valid for the given field names.  Unmatched fields: " +
                    Sets.symmetricDifference(metadataFieldNames, funcotationFieldNames).stream().collect(Collectors.joining(", ")));
        }
    }

    //==================================================================================================================
    // Override Methods:

    @Override
    public Allele getAltAllele() {
        return altAllele;
    }

    @Override
    public String getDataSourceName() {
        return dataSourceName;
    }

    @Override
    public void setFieldSerializationOverrideValue(final String fieldName, final String overrideValue) {
        if ( !fieldMap.containsKey(fieldName) ) {
            throw new GATKException("Attempted to override a field that is not contained in this TableFuncotation: "
                    + fieldName + " is not one of [" + String.join(",", fieldMap.keySet()) + "]");
        }
        fieldMap.put(fieldName, overrideValue);
    }

    @Override
    public String serializeToVcfString() {
        return fieldMap.values().stream()
                .map(f -> (f == null ? "" : FuncotatorUtils.sanitizeFuncotationForVcf(f)))
                .collect(Collectors.joining(VcfOutputRenderer.FIELD_DELIMITER));
    }

    @Override
    public LinkedHashSet<String> getFieldNames() {
        return new LinkedHashSet<>(fieldMap.keySet());
    }

    @Override
    public String getField(final String fieldName) {
        if ( fieldMap.containsKey(fieldName) ) {
            return fieldMap.get(fieldName);
        }
        else {
            throw new GATKException(this.getClass().getSimpleName() + ": Does not contain field: " + fieldName);
        }
    }

    @Override
    public boolean hasField(final String fieldName) {
        return fieldMap.containsKey(fieldName);
    }

    @Override
    public FuncotationMetadata getMetadata() {
        return this.metadata;
    }

    @Override
    public boolean equals(final Object o) {
        if ( this == o ) return true;
        if ( o == null || getClass() != o.getClass() ) return false;

        final TableFuncotation that = (TableFuncotation) o;

        if ( dataSourceName != null ? !dataSourceName.equals(that.dataSourceName) : that.dataSourceName != null )
            return false;
        if ( fieldMap != null ? !fieldMap.equals(that.fieldMap) : that.fieldMap != null ) return false;
        return altAllele != null ? altAllele.equals(that.altAllele) : that.altAllele == null;
    }

    @Override
    public int hashCode() {
        int result = dataSourceName != null ? dataSourceName.hashCode() : 0;
        result = 31 * result + (fieldMap != null ? fieldMap.hashCode() : 0);
        result = 31 * result + (altAllele != null ? altAllele.hashCode() : 0);
        return result;
    }

    @Override
    public String toString() {
        return "TableFuncotation{" +
                "dataSourceName='" + dataSourceName + '\'' +
                ", fieldMap={" + fieldMap.keySet().stream().map(k -> k + ":" + fieldMap.get(k)).collect(Collectors.joining(" , ")) + '}' +
                ", altAllele=" + altAllele +
                '}';
    }


    //==================================================================================================================
    // Static Methods:

    /**
     * Create a TableFuncotation with given name value pairs.
     *
     * @param fieldNames Names corresponding to the values in <code>fieldValues</code>.  Never {@code null}
     * @param fieldValues Values corresponding to the names in <code>fieldNames</code>.Never {@code null}
     * @param altAllele Alternate allele to use for all of the fields in this funcotation.  Never {@code null}
     * @param dataSourceName Datasource name to use for all of the fields in this funcotation.  Never {@code null}
     * @param metadata Metadata to use for the given fields.  Use {@code null} to get metadata with unknown description.
     * @return Never {@code null}
     */
    public static TableFuncotation create(final List<String> fieldNames, final List<String> fieldValues, final Allele altAllele, final String dataSourceName, final FuncotationMetadata metadata ) {
        Utils.nonNull(fieldNames);
        Utils.nonNull(fieldValues);
        Utils.nonNull(altAllele);
        Utils.nonNull(dataSourceName);
        return new TableFuncotation(fieldNames, fieldValues, altAllele, dataSourceName, metadata);
    }

    /**
     * Convenience method to create a table funcotation from a {@link LinkedHashSet}, which preserves order.
     *
     * See {@link TableFuncotation#create(List, List, Allele, String, FuncotationMetadata)}
     *
     * @param fieldNames See {@link TableFuncotation#create(List, List, Allele, String, FuncotationMetadata)}
     * @param fieldValues See {@link TableFuncotation#create(List, List, Allele, String, FuncotationMetadata)}
     * @param altAllele See {@link TableFuncotation#create(List, List, Allele, String, FuncotationMetadata)}
     * @param dataSourceName See {@link TableFuncotation#create(List, List, Allele, String, FuncotationMetadata)}
     * @param metadata See {@link TableFuncotation#create(List, List, Allele, String, FuncotationMetadata)}
     * @return See {@link TableFuncotation#create(List, List, Allele, String, FuncotationMetadata)}
     */
    public static TableFuncotation create(final LinkedHashSet<String> fieldNames, final List<String> fieldValues, final Allele altAllele, final String dataSourceName, final FuncotationMetadata metadata ) {
        return create(new ArrayList<>(fieldNames), fieldValues, altAllele, dataSourceName, metadata );
    }

    /**
     *  See {@link TableFuncotation#create(List, List, Allele, String, FuncotationMetadata)}
     *
     * @param data Map for field name to field value.  The field value will be converted to a String.  Never {@code null}
     * @param altAllele See {@link TableFuncotation#create(List, List, Allele, String, FuncotationMetadata)}
     * @param dataSourceName See {@link TableFuncotation#create(List, List, Allele, String, FuncotationMetadata)}
     * @param metadata See {@link TableFuncotation#create(List, List, Allele, String, FuncotationMetadata)}
     * @return See {@link TableFuncotation#create(List, List, Allele, String, FuncotationMetadata)}
     */
    public static TableFuncotation create(final Map<String, Object> data, final Allele altAllele, final String dataSourceName, final FuncotationMetadata metadata ) {
        final List<String> fieldNames = new ArrayList<>(data.keySet());
        final List<String> fieldValues = fieldNames.stream().map(f -> data.get(f).toString()).collect(Collectors.toList());
        return create(fieldNames, fieldValues, altAllele, dataSourceName, metadata);
    }

    /**
     * See {@link TableFuncotation#create(List, List, Allele, String, FuncotationMetadata)}
     *
     * @param xsvTableFeature {@link XsvTableFeature} with field names and values fully populated.  Never {@code null}
     * @param altAllele See {@link #create(List, List, Allele, String, FuncotationMetadata)}
     * @param dataSourceName See {@link #create(List, List, Allele, String, FuncotationMetadata)}
     * @param metadata See {@link #create(List, List, Allele, String, FuncotationMetadata)}
     * @return See {@link #create(List, List, Allele, String, FuncotationMetadata)}
     */
    public static TableFuncotation create(final XsvTableFeature xsvTableFeature, final Allele altAllele, final String dataSourceName, final FuncotationMetadata metadata ) {
        Utils.nonNull(xsvTableFeature);
        final List<String> fieldNames = xsvTableFeature.getHeaderWithoutLocationColumns();
        final List<String> fieldValues = xsvTableFeature.getValuesWithoutLocationColumns();
        return create(fieldNames, fieldValues, altAllele, dataSourceName, metadata);
    }

    //==================================================================================================================
    // Instance Methods:

    /**
     * Get the value in this {@link TableFuncotation} corresponding to the given key.
     * If the key is not in this {@link TableFuncotation}, returns {@code null}.
     * @param key Key to get from this {@link TableFuncotation}.
     * @return The value corresponding to the given key or {@code null}.
     */
    public String get(final String key) {
        return fieldMap.get(key);
    }

    /**
     * @return The {@link Set} of field names in this {@link TableFuncotation}.
     */
    public Set<String> keySet() {
        return fieldMap.keySet();
    }

    /**
     * @return The {@link Collection} of field values in this {@link TableFuncotation}.
     */
    public Collection<String> values() {
        return fieldMap.values();
    }

    /**
     * @return The number of key-value pairs in this {@link TableFuncotation}.
     */
    public int size() {
        return fieldMap.size();
    }

    //==================================================================================================================
    // Helper Data Types:
}
