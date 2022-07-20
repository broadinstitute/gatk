package org.broadinstitute.hellbender.tools.walkers.sv.GqRecalibrator;

import net.minidev.json.JSONArray;
import net.minidev.json.JSONObject;
import net.minidev.json.JSONValue;
import net.minidev.json.parser.ParseException;
import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.jetbrains.annotations.NotNull;

import java.io.*;
import java.math.BigDecimal;
import java.nio.charset.StandardCharsets;
import java.util.*;
import java.util.function.ToIntFunction;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Class to manage table with mixed columnar and matrix properties; with different primitive types. Supported types:
 *     boolean, int, long, float, double
 * <p>
 * This class is useful for cases when the exact properties are not necessarily fixed in advance and it is desirable to
 * build records dynamically from (name, value) pairs.
 * <p>
 * Conceptually the table is organized as numRows x numProperties x numColumns
 * Properties are stored as Object that can be cast to primitive arrays or array-of-arrays (matrix) in their original
 *     primitive type
 * Rows are extracted as float[] (suitable for machine learning), with one of two options
 *    normalize = false: translated to float type but otherwise unchanged from raw values
 *    normalize = true: shifted so the median is 0, and scaled so the standard deviation (over the central half of the
 *                      data) is 1.0. This is done on the fly when extracting rows.
 *                      The baseline and scale can be provided when adding a column (to provide consistency between
 *                      training and inference data) or computed automatically.
 */
class PropertiesTable implements Iterable<PropertiesTable.Property> {
    private static final String PROPERTY_NAMES_KEY = "propertyNames";
    private static final String PROPERTY_CLASSES_KEY = "propertyClasses";
    private static final String PROPERTY_BASELINE_KEY = "propertyBaseline";
    private static final String PROPERTY_SCALE_KEY = "propertyScale";
    private static final int DEFAULT_INITIAL_NUM_ALLOCATED_ROWS = 1000;
    private boolean allNumeric = true;

    enum PropertyClass {
        BooleanArrProperty, BooleanMatProperty,
        ByteArrProperty, ByteMatProperty,
        ShortArrProperty, ShortMatProperty,
        IntArrProperty, IntMatProperty,
        LongArrProperty, LongMatProperty,
        FloatArrProperty, FloatMatProperty,
        DoubleArrProperty, DoubleMatProperty,
        StringArrProperty, StringMatProperty,
        StringSetArrProperty, StringSetMatProperty,
    }

    //private final List<Property> properties = new ArrayList<>();
    private final Map<String, Property> properties = new HashMap<>();
    private final Map<String, List<String>> labelsEncoding = new HashMap<>();
    private List<String> orderedPropertyNames = null;
    private float[] propertiesRow = null;
    private int initialNumAllocatedRows;

    PropertiesTable(final int initialNumAllocatedRows) { this.initialNumAllocatedRows = initialNumAllocatedRows; }
    PropertiesTable() { this(DEFAULT_INITIAL_NUM_ALLOCATED_ROWS); }

    @NotNull @Override public Iterator<Property> iterator() {
        if(orderedPropertyNames == null) {
            setOrderedPropertyNames();
        }
        return orderedPropertyNames.stream().map(properties::get).iterator();
    }

    public List<String> getPropertyNames() {
        if(orderedPropertyNames == null) {
            setOrderedPropertyNames();
        }
        return orderedPropertyNames;
    }

    private void setOrderedPropertyNames() {
        orderedPropertyNames = properties.keySet().stream().sorted().collect(Collectors.toList());
    }

    public Property get(final String propertyName) {
        return properties.get(propertyName);
    }

    public Property getOrCreateProperty(final String propertyName, final PropertyClass propertyClass) {
        return getOrCreateProperty(propertyName, propertyClass, initialNumAllocatedRows);
    }

    public Property getOrCreateProperty(final String propertyName, final PropertyClass propertyClass,
                                        final int numAllocatedRows) {
        Property property = properties.getOrDefault(propertyName, null);
        if(property == null) {
            property = Property.create(propertyClass, propertyName, numAllocatedRows);
            if(!property.isNumeric()) {
                allNumeric = false;
            }
            properties.put(propertyName, property);
        } else if(property.getAllocatedRows() < numAllocatedRows) {
            property.setAllocatedRows(numAllocatedRows);
        }
        return property;
    }

    public void setBaselineAndScale(final String propertyName, final PropertyClass propertyClass,
                                    final double baseline, final double scale) {
        getOrCreateProperty(propertyName, propertyClass).setBaselineAndScale(baseline, scale);
    }

    public void addLabelEncoding(final String propertyName, final List<String> propertyLabels) {
        labelsEncoding.put(propertyName, propertyLabels.stream().sorted().collect(Collectors.toList()));
    }

    public void remove(final String propertyName) {
        properties.remove(propertyName);
    }

    public void insert(final Property property) {
        if(!property.isNumeric()) {
            allNumeric = false;
        }
        if (properties.containsKey(property.name)) {
            final Property oldProperty = properties.get(property.name);
            if(oldProperty.getNumRows() > 0) {
                throw new IllegalArgumentException("PropertiesTable already contains property: " + property.name);
            } else {
                // property exists but it's empty. If it has useful baseline and scale information, transfer those,
                // then overwrite with new property
                if(oldProperty.normalizationIsSet()) {
                    property.setBaselineAndScale(oldProperty.getBaseline(), oldProperty.getScale());
                }
                properties.put(property.name, property);
            }
        } else {
            properties.put(property.name, property);
        }
    }

    /**
     * "clear" existing data by setting numRows to 0. Property types and memory allocation is not altered so that any
     * subsequent data can be written to first row without continual reallocation
     */
    public void clearRows() {
        setNumRows(0);
    }

    public void setNumRows(final int numRows) {
        for(final Property property : properties.values()) {
            property.numRows = numRows;
        }
    }

    public void setNumAllocatedRows(final int numRows) {
        initialNumAllocatedRows = numRows;
        for(final Property property : properties.values()) {
            property.setAllocatedRows(numRows);
        }
    }

    public void append(final String propertyName, final boolean value) {
        getOrCreateProperty(propertyName, PropertyClass.BooleanArrProperty).append(value);
    }
    public void append(final String propertyName, final boolean[] value) {
        getOrCreateProperty(propertyName, PropertyClass.BooleanMatProperty).append(value);
    }
    public void append(final String propertyName, final byte value) {
        getOrCreateProperty(propertyName, PropertyClass.ByteArrProperty).append(value);
    }
    public void append(final String propertyName, final byte[] value) {
        getOrCreateProperty(propertyName, PropertyClass.ByteMatProperty).append(value);
    }
    public void append(final String propertyName, final short value) {
        getOrCreateProperty(propertyName, PropertyClass.ShortArrProperty).append(value);
    }
    public void append(final String propertyName, final short[] value) {
        getOrCreateProperty(propertyName, PropertyClass.ShortMatProperty).append(value);
    }
    public void append(final String propertyName, final int value) {
        getOrCreateProperty(propertyName, PropertyClass.IntArrProperty).append(value);
    }
    public void append(final String propertyName, final int[] value) {
        getOrCreateProperty(propertyName, PropertyClass.IntMatProperty).append(value);
    }
    public void append(final String propertyName, final long value) {
        getOrCreateProperty(propertyName, PropertyClass.LongArrProperty).append(value);
    }
    public void append(final String propertyName, final long[] value) {
        getOrCreateProperty(propertyName, PropertyClass.LongMatProperty).append(value);
    }
    public void append(final String propertyName, final float value) {
        getOrCreateProperty(propertyName, PropertyClass.FloatArrProperty).append(value);
    }
    public void append(final String propertyName, final float[] value) {
        getOrCreateProperty(propertyName, PropertyClass.FloatMatProperty).append(value);
    }
    public void append(final String propertyName, final double value) {
        getOrCreateProperty(propertyName, PropertyClass.DoubleArrProperty).append(value);
    }
    public void append(final String propertyName, final double[] value) {
        getOrCreateProperty(propertyName, PropertyClass.DoubleMatProperty).append(value);
    }
    public void append(final String propertyName, final String value) {
        getOrCreateProperty(propertyName, PropertyClass.StringArrProperty).append(value);
    }
    public void append(final String propertyName, final String[] value) {
        getOrCreateProperty(propertyName, PropertyClass.StringMatProperty).append(value);
    }
    public void append(final String propertyName, final Set<String> value) {
        getOrCreateProperty(propertyName, PropertyClass.StringSetArrProperty).append(value);
    }
    public void append(final String propertyName, final Set<String>[] value) {
        getOrCreateProperty(propertyName, PropertyClass.StringSetMatProperty).append(value);
    }

    public void set(final String propertyName, final boolean[] values) {
        getOrCreateProperty(propertyName, PropertyClass.BooleanArrProperty).set(values);
    }
    public void set(final String propertyName, final boolean[][] values) {
        getOrCreateProperty(propertyName, PropertyClass.BooleanMatProperty).set(values);
    }
    public void set(final String propertyName, final byte[] values) {
        getOrCreateProperty(propertyName, PropertyClass.ByteArrProperty).set(values);
    }
    public void set(final String propertyName, final byte[][] values) {
        getOrCreateProperty(propertyName, PropertyClass.ByteMatProperty).set(values);
    }
    public void set(final String propertyName, final short[] values) {
        getOrCreateProperty(propertyName, PropertyClass.ShortArrProperty).set(values);
    }
    public void set(final String propertyName, final short[][] values) {
        getOrCreateProperty(propertyName, PropertyClass.ShortMatProperty).set(values);
    }
    public void set(final String propertyName, final int[] values) {
        getOrCreateProperty(propertyName, PropertyClass.IntArrProperty).set(values);
    }
    public void set(final String propertyName, final int[][] values) {
        getOrCreateProperty(propertyName, PropertyClass.IntMatProperty).set(values);
    }
    public void set(final String propertyName, final long[] values) {
        getOrCreateProperty(propertyName, PropertyClass.LongArrProperty).set(values);
    }
    public void set(final String propertyName, final long[][] values) {
        getOrCreateProperty(propertyName, PropertyClass.LongMatProperty).set(values);
    }
    public void set(final String propertyName, final float[] values) {
        getOrCreateProperty(propertyName, PropertyClass.FloatArrProperty).set(values);
    }
    public void set(final String propertyName, final float[][] values) {
        getOrCreateProperty(propertyName, PropertyClass.FloatMatProperty).set(values);
    }
    public void set(final String propertyName, final double[] values) {
        getOrCreateProperty(propertyName, PropertyClass.DoubleArrProperty).set(values);
    }
    public void set(final String propertyName, final double[][] values) {
        getOrCreateProperty(propertyName, PropertyClass.DoubleMatProperty).set(values);
    }
    public void set(final String propertyName, final String[] values) {
        getOrCreateProperty(propertyName, PropertyClass.StringArrProperty).set(values);
    }
    public void set(final String propertyName, final String[][] values) {
        getOrCreateProperty(propertyName, PropertyClass.StringMatProperty).set(values);
    }
    public void set(final String propertyName, final Set<String>[] values) {
        getOrCreateProperty(propertyName, PropertyClass.StringSetArrProperty).set(values);
    }
    public void set(final String propertyName, final Set<String>[][] values) {
        getOrCreateProperty(propertyName, PropertyClass.StringSetMatProperty).set(values);
    }

    protected int getConsistentIntPropertyValue(ToIntFunction<Property> method, final String valuesLabel) {
        final int[] distinctValues = properties.values().stream()
                .mapToInt(method)
                .filter(c -> c >= 0)
                .distinct()
                .toArray();
        switch(distinctValues.length) {
            case 0: // no values
                return -1;
            case 1: // consistent value
                return distinctValues[0];
            default: // inconsistent values
                System.out.println(valuesLabel + ":");
                properties.forEach(
                    (propertyName, property) -> System.out.println("\t" + propertyName + ": " +
                                                                   method.applyAsInt(property))
                );
                throw new IllegalArgumentException("PropertiesTable contains inconsistent numbers of " + valuesLabel);
        }
    }

    @SuppressWarnings("UnusedReturnValue")
    public int getNumColumns() {
        return getConsistentIntPropertyValue(Property::getNumColumns, "columns");
    }

    public int getNumRows() {
        return getConsistentIntPropertyValue(Property::getNumRows, "rows");
    }

    public int getNumProperties() {
        if(allNumeric) {
            return properties.size();
        }
        return properties.values().stream()
                .mapToInt(
                        property -> {
                            final List<String> allLabels = labelsEncoding.containsKey(property.name) ?
                                    labelsEncoding.get(property.name) :
                                    property.getAllLabels();
                            return allLabels == null ?
                                    1 :
                                    allLabels.size() > 2 ?
                                            allLabels.size() :
                                            allLabels.size() == 2 ? 1 : 0;
                        }
                )
                .sum();
    }

    protected void trim() {
        // save memory by right-sizing the dynamically allocated arrays
        for(final Property property : properties.values()){
            property.setAllocatedRows(property.getNumRows());
        }
    }

    protected void oneHot() {
        // Can't edit properties while iterating, so form a copy first
        for(final Property property : new ArrayList<>(properties.values())) {
            final List<String> allLabels = labelsEncoding.containsKey(property.name) ?
                labelsEncoding.get(property.name) :
                property.getAllLabels();
            if(allLabels != null) {
                if(!labelsEncoding.containsKey(property.name)) {
                    addLabelEncoding(property.name, allLabels);
                }
                remove(property.name);
                property.oneHot(allLabels, this);
            }
        }
        allNumeric = true;
    }

    protected void setBaselineAndScales() {
        final List<String> setProperties = new ArrayList<>();
        final List<String> previouslySetProperties = new ArrayList<>();
        for(final Property property : properties.values()) {
            if(!property.normalizationIsSet()) {
                property.calculateBaselineAndScale();
                setProperties.add(property.name);
            } else {
                previouslySetProperties.add(property.name);
            }
        }
        if(setProperties.size() > 0 && previouslySetProperties.size() > 0) {
            // Should be that either the baseline and scales were all calculated previously, or none have been
            // calculated. If some have and some haven't, state has gotten messed up somehow.
            throw new IllegalArgumentException(
                "Some but not all properties have baseline and scale set:" +
                "\n\tpreviously set:" + String.join(", ", previouslySetProperties) +
                "\n\tpreviously unset:" + String.join(", ", setProperties)
            );
        }
    }

    /**
     * Once all data has been added to the table
     *   1) check to make sure number of rows and columns are consistent
     *   2) right-size memory (match number of allocated and actual rows)
     *   3) one-hot encode any label properties
     *   4) ensure that baseline and scale are appropriately set for every Property
     *   5) allocate float[] propertiesRow so that rows can be efficiently queried for inference
     */
    public void validateAndFinalize() {
        oneHot();
        setOrderedPropertyNames();
        getNumColumns();
        getNumRows();
        trim();
        setBaselineAndScales();
        propertiesRow = new float[getNumProperties()];
    }

    /**
     * Copy table data from specified row into prepared array buffer
     * @param outArray array buffer to copy into
     * @param outIndex offset of buffer to start copying
     * @param rowIndex requested row from table
     * @param columnIndex requested column from table (properties without columns will ignore this value)
     * @param normalize if true, normalize (stable z-score) data, if false, copy raw data
     * @return next outIndex for writing
     */
    public int copyPropertiesRow(final float[] outArray, int outIndex,
                                 final int rowIndex, final int columnIndex, final boolean normalize) {
        if(orderedPropertyNames == null) {
            setOrderedPropertyNames();
        }
        for(final String propertyName : orderedPropertyNames) {
            final Property property = properties.get(propertyName);
            outArray[outIndex] = property.getAsFloat(rowIndex, columnIndex, normalize);
            if(!Float.isFinite(outArray[outIndex])) {
                System.out.println(propertyName + " = " + outArray[outIndex] + " at row " + rowIndex + ", column " + columnIndex);
            }
            ++outIndex;
        }
        return outIndex;
    }

    public float[] getPropertiesRow(final int rowIndex, final int columnIndex, final boolean normalize) {
        copyPropertiesRow(propertiesRow, 0, rowIndex, columnIndex, normalize);
        return propertiesRow;
    }

    public Map<String, List<String>> getLabelsEncoding() {
        // don't want anyone messing with this, just want to expose the information
        return Collections.unmodifiableMap(labelsEncoding);
    }

    protected void saveNormalizationStats(final OutputStream outputStream) {
        final JSONArray propNames = new JSONArray();
        final JSONArray propClasses = new JSONArray();
        final JSONArray propBase = new JSONArray();
        final JSONArray propScale = new JSONArray();
        if(orderedPropertyNames == null) {
            setOrderedPropertyNames();
        }
        for (final String propertyName : orderedPropertyNames) {
            final Property property = properties.get(propertyName);
            propNames.add(propertyName);
            propClasses.add(property.getClass().getSimpleName());
            propBase.add(property.baseline);
            propScale.add(property.scale);
        }
        final JSONObject jsonObject = new JSONObject();
        jsonObject.put(PROPERTY_NAMES_KEY, propNames);
        jsonObject.put(PROPERTY_CLASSES_KEY, propClasses);
        jsonObject.put(PROPERTY_BASELINE_KEY, propBase);
        jsonObject.put(PROPERTY_SCALE_KEY, propScale);

        try {
            outputStream.write(jsonObject.toJSONString().getBytes(StandardCharsets.UTF_8));
            outputStream.write("\n".getBytes(StandardCharsets.UTF_8));
        } catch(IOException ioException) {
            throw new GATKException("Error saving data summary json", ioException);
        }
    }

    protected void saveLabelsEncoding(final OutputStream outputStream) {
        final JSONObject jsonObject = new JSONObject();
        for(final Map.Entry<String, List<String>> encodingEntry : labelsEncoding.entrySet()) {
            final JSONArray labels = new JSONArray();
            labels.addAll(encodingEntry.getValue());
            jsonObject.put(encodingEntry.getKey(), labels);
        }
        try {
            outputStream.write(jsonObject.toJSONString().getBytes(StandardCharsets.UTF_8));
            outputStream.write("\n".getBytes(StandardCharsets.UTF_8));
        } catch(IOException ioException) {
            throw new GATKException("Error saving data summary json", ioException);
        }

    }

    public void saveDataEncoding(final OutputStream outputStream) throws IOException {
        saveNormalizationStats(outputStream);
        saveLabelsEncoding(outputStream);
    }

    public void saveTable(final OutputStream outputStream) throws IOException {
        if(orderedPropertyNames == null) {
            setOrderedPropertyNames();
        }
        final JSONObject jsonObject = new JSONObject();
        for(final String propertyName : orderedPropertyNames) {
            final Property property = properties.get(propertyName);
            jsonObject.put(property.name, property.getJSON());
        }
        outputStream.write(jsonObject.toJSONString().getBytes());
    }

    public void save(final OutputStream outputStream) throws IOException {
        saveDataEncoding(outputStream);
        outputStream.write("\n".getBytes());
        saveTable(outputStream);
    }

    protected static long getNumTrue(final boolean[] values) {
        return getNumTrue(values, values.length);
    }
    protected static long getNumTrue(final boolean[] values, final int numRows) {
        long numTrue = 0;
        for(int i = 0; i < numRows; ++i) {
            if(values[i]) {
                ++numTrue;
            }
        }
        return numTrue;
    }

    protected static double getBaselineOrdered(final double[] orderedValues) {
        // get baseline as median of values
        return orderedValues.length == 0 ?
                0 :
                orderedValues.length % 2 == 1 ?
                        orderedValues[orderedValues.length / 2] :
                        (orderedValues[orderedValues.length / 2 - 1] + orderedValues[orderedValues.length / 2]) / 2.0;
    }

    protected static double getScaleOrdered(final double[] orderedValues, final double baseline) {
        // get scale as root-mean-square difference from baseline, over central half of data (to exclude outliers)
        switch(orderedValues.length) {
            case 0:
            case 1:
                return 1.0;
            default:
                final int q1 = orderedValues.length / 4;
                // this formula can have early problems with indices larger than 2^31
                // final int q3 = 3 * orderedValues.length / 4;
                // this is equivalent
                final int q3 = orderedValues.length - 1 - q1;
                final int start, stop;
                if(orderedValues[q1] == orderedValues[q3]) {
                    // the center of the distribution is effectively discrete, quantiles won't work
                    start = 0; stop = orderedValues.length;
                } else {
                    start = q1; stop = q3;
                }
                double scale = 0.0;
                for(int idx = start; idx < stop; ++idx) {
                    scale += (orderedValues[idx] - baseline) * (orderedValues[idx] - baseline);
                }
                return FastMath.max(FastMath.sqrt(scale / (1 + stop - start)), 1.0e-6);
        }
    }

    protected static double getDoubleFromJSON(final Object jsonObject) {
        if(jsonObject instanceof Double) {
            return (Double) jsonObject;
        } else if(jsonObject instanceof Float) {
            return (Float) jsonObject;
        } else if(jsonObject instanceof BigDecimal) {
            return ((BigDecimal)jsonObject).doubleValue();
        } else {
            throw new GATKException("Unknown conversion to double for " + jsonObject.getClass().getName());
        }
    }


    final String getNextLine(final InputStream inputStream) {
        try {
            byte[] bytes = new byte[1000];
            int ind = 0;
            for(int c = inputStream.read(); (byte)c != '\n' && c != -1; c = inputStream.read()) {
                bytes[ind] = (byte)c;
                ++ind;
                if(ind == bytes.length) {
                    bytes = Arrays.copyOf(bytes, bytes.length * 2);
                }
            }
            return new String(bytes, StandardCharsets.UTF_8);
        } catch (IOException ioException) {
            throw new GATKException("Unable to get next line from inputStream", ioException);
        }
    }

    final JSONObject getNextJSONObject(final InputStream inputStream) {
        final String nextLine = getNextLine(inputStream);
        try {
            return (JSONObject) JSONValue.parseWithException(nextLine);
        } catch (ParseException | ClassCastException parseException) {
            throw new GATKException("Unable to parse JSON from inputStream", parseException);
        }
    }

    protected void loadNormalizations(final InputStream inputStream) {
        final JSONObject jsonObject = getNextJSONObject(inputStream);
        final JSONArray propNames = ((JSONArray) jsonObject.get(PROPERTY_NAMES_KEY));
        final JSONArray propClasses = ((JSONArray) jsonObject.get(PROPERTY_CLASSES_KEY));
        final JSONArray propBase = ((JSONArray) jsonObject.get(PROPERTY_BASELINE_KEY));
        final JSONArray propScale = ((JSONArray) jsonObject.get(PROPERTY_SCALE_KEY));
        for (int idx = 0; idx < propNames.size(); ++idx) {
            setBaselineAndScale(
                    (String)propNames.get(idx), PropertyClass.valueOf((String)propClasses.get(idx)),
                    getDoubleFromJSON(propBase.get(idx)), getDoubleFromJSON(propScale.get(idx))
            );
        }
    }

    protected void loadEncodings(final InputStream inputStream) {
        final JSONObject jsonObject = getNextJSONObject(inputStream);
        for(final String propertyName : jsonObject.keySet()) {
            final JSONArray labels = (JSONArray) jsonObject.get(propertyName);
            addLabelEncoding(propertyName, labels.stream().map(label -> (String)label).collect(Collectors.toList()));
        }
    }

    public void loadDataEncoding(final InputStream inputStream) {
        loadNormalizations(inputStream);
        loadEncodings(inputStream);
    }

    static abstract class Property {
        public final String name;
        private Float baseline = null;
        private Float scale = null;
        protected int numRows;

        static final int ALLOCATION_GROWTH_SCALE = 2;

        Property(final String name) { this.name = name; }

        abstract Property set(final Object values);
        abstract public float getAsFloat(final int rowIndex, final int hyperIndex);
        abstract protected int getAllocatedRows();
        abstract public void setAllocatedRowsUnguarded(final int numRows);
        abstract protected void assignNextValue(final Object value);
        abstract protected double[] getValuesAsOrderedDoubles();

        public Float getBaseline() { return baseline; }
        public Float getScale() { return scale; }
        public boolean normalizationIsSet() {
            if(baseline == null || scale == null) {
                if(!(baseline == null && scale == null)) {
                    throw new IllegalArgumentException("baseline or scale assigned a null value");
                }
                return false;
            } else {
                return true;
            }
        }
        public boolean isNumeric() { return true; }
        public List<String> getAllLabels() { return null; }
        public void oneHot(final List<String> allLabels, final PropertiesTable propertiesTable) {}

        public static Property create(final PropertyClass propertyClass, final String name) {
            return create(propertyClass, name, DEFAULT_INITIAL_NUM_ALLOCATED_ROWS);
        }

        public static Property create(final PropertyClass propertyClass, final String name,
                                      final int numRows) {
            switch(propertyClass) {
                case BooleanArrProperty: return new BooleanArrProperty(name, numRows);
                case BooleanMatProperty: return new BooleanMatProperty(name, numRows);
                case ByteArrProperty: return new ByteArrProperty(name, numRows);
                case ByteMatProperty: return new ByteMatProperty(name, numRows);
                case ShortArrProperty: return new ShortArrProperty(name, numRows);
                case ShortMatProperty: return new ShortMatProperty(name, numRows);
                case IntArrProperty: return new IntArrProperty(name, numRows);
                case IntMatProperty: return new IntMatProperty(name, numRows);
                case LongArrProperty: return new LongArrProperty(name, numRows);
                case LongMatProperty: return new LongMatProperty(name, numRows);
                case FloatArrProperty: return new FloatArrProperty(name, numRows);
                case FloatMatProperty: return new FloatMatProperty(name, numRows);
                case DoubleArrProperty: return new DoubleArrProperty(name, numRows);
                case DoubleMatProperty: return new DoubleMatProperty(name, numRows);
                case StringArrProperty: return new StringArrProperty(name, numRows);
                case StringMatProperty: return new StringMatProperty(name, numRows);
                case StringSetArrProperty: return new StringSetArrProperty(name, numRows);
                case StringSetMatProperty: return new StringSetMatProperty(name, numRows);
                default: throw new IllegalArgumentException("Unable to create Property of type " + propertyClass);
            }
        }

        public int getNumColumns() {
            final int[] numColumns = getArrayOfDistinctNumColumns();
            switch(numColumns.length) {
                case 0: // no data in property or an Array property
                    return -1;
                case 1: // normal case for matrix property
                    return numColumns[0];
                default: // error, inconsistent numbers of columns
                    throw new IllegalArgumentException("Inconsistent number of columns in matrix property: " + name);
            }
        }
        protected int[] getArrayOfDistinctNumColumns() { return new int[0]; }

        public int getNumRows() { return numRows; }

        public void setAllocatedRows(final int numRows) {
            if(numRows <= getNumRows()) {
                if(numRows < getNumRows()) {
                    throw new IllegalArgumentException("Can't shrink allocated memory for property to less than used space.");
                }
            } else {
                setAllocatedRowsUnguarded(numRows);
            }
        }

        public void clearAllocatedRows() {
            setAllocatedRowsUnguarded(0);
        }

        public float getAsFloat(final int rowIndex, final int hyperIndex, final boolean normalize) {
            final float rawValue = getAsFloat(rowIndex, hyperIndex);
            return normalize ?
                    (rawValue - baseline) / scale :
                    rawValue;
        }
        public boolean getAsBool(final int rowIndex, final int hyperIndex) {
            throw new IllegalArgumentException(
                "getAsBool() not defined for Property of type " + this.getClass().getSimpleName()
            );
        }
        public int getAsInt(final int rowIndex, final int hyperIndex) {
            throw new IllegalArgumentException(
                "getAsInt() not defined for Property of type " + this.getClass().getSimpleName()
            );
        }
        public long getAsLong(final int rowIndex, final int hyperIndex) {
            throw new IllegalArgumentException(
                "getAsLong() not defined for Property of type " + this.getClass().getSimpleName()
            );
        }

        abstract public JSONArray getJSON();

        public void append(final Object value) {
            if(getNumRows() >= getAllocatedRows()) {
                setAllocatedRows(FastMath.max(DEFAULT_INITIAL_NUM_ALLOCATED_ROWS, getNumRows() * ALLOCATION_GROWTH_SCALE));
            }
            assignNextValue(value);
            ++numRows;
        }

        public void setBaselineAndScale(final double baseline, final double scale) {
            this.baseline = (float)baseline;
            this.scale = (float)scale;
        }

        protected void calculateBaselineAndScale() {
            final double[] orderedValues = getValuesAsOrderedDoubles();
            final double baseline = getBaselineOrdered(orderedValues);
            this.baseline = (float)baseline;
            scale = (float)getScaleOrdered(orderedValues, baseline);
        }
    }

    static public class BooleanArrProperty extends Property {
        public boolean[] values;

        BooleanArrProperty(final String name, final int numValues) { super(name); values = new boolean[numValues]; }
        BooleanArrProperty(final String name) { this(name, DEFAULT_INITIAL_NUM_ALLOCATED_ROWS); }

        @Override public Property set(Object values) { this.values = (boolean[]) values; return this; }
        @Override public JSONArray getJSON() {
            final JSONArray jsonArray = new JSONArray();
            IntStream.range(0, numRows).forEach(i -> jsonArray.add(values[i]));
            return jsonArray;
        }
        @Override public float getAsFloat(final int rowIndex, final int hyperIndex) {
            return values[rowIndex] ? 1F : 0F;
        }
        @Override public boolean getAsBool(final int rowIndex, final int hyperIndex) { return values[rowIndex]; }
        @Override public int getAllocatedRows() { return values.length; }
        @Override public void setAllocatedRowsUnguarded(final int numRows) {
            values = Arrays.copyOf(values, numRows);
        }
        @Override protected void assignNextValue(final Object value) { values[numRows] = (boolean)value; }
        @Override protected void calculateBaselineAndScale() {
            // special case for booleans
            final long numTrue = getNumTrue(values, numRows);
            final double baseline = numTrue / (double) numRows;
            final double scale = numTrue == 0 || numTrue == numRows ?
                    1.0 :
                    FastMath.sqrt(baseline * (1.0 - baseline));
            setBaselineAndScale(baseline, scale);
        }
        @Override protected double[] getValuesAsOrderedDoubles() {
            // special case for booleans
            throw new GATKException("Method not needed for booleans");
        }
    }

    static public class BooleanMatProperty extends Property {
        public boolean[][] values;

        BooleanMatProperty(final String name, final int numValues) { super(name); values = new boolean[numValues][]; }
        BooleanMatProperty(final String name) { this(name, DEFAULT_INITIAL_NUM_ALLOCATED_ROWS); }

        @Override public Property set(Object values) { this.values = (boolean[][]) values; return this; }
        @Override public JSONArray getJSON() {
            final JSONArray jsonArray = new JSONArray();
            IntStream.range(0, numRows).forEach(i -> jsonArray.add(values[i]));
            return jsonArray;
        }
        @Override public float getAsFloat(final int rowIndex, final int hyperIndex) {
            return values[rowIndex][hyperIndex] ? 1F : 0F;
        }
        @Override public boolean getAsBool(final int rowIndex, final int hyperIndex) {
            return values[rowIndex][hyperIndex];
        }
        @Override protected int[] getArrayOfDistinctNumColumns() {
            return IntStream.range(0, numRows).map(i -> values[i].length).distinct().toArray();
        }
        @Override public int getAllocatedRows() { return values.length; }
        @Override public void setAllocatedRowsUnguarded(final int numRows) {
            values = Arrays.copyOf(values, numRows);
        }
        @Override protected void assignNextValue(final Object value) { values[numRows] = (boolean[])value; }
        protected void assignColumnValue(final int column, final boolean value) {
            values[numRows][column] = value;
        }
        @Override protected void calculateBaselineAndScale() {
            // special case for booleans
            final long numTrue = IntStream.range(0, numRows).mapToLong(i -> getNumTrue(values[i])).sum();
            final double baseline = numTrue / (double) numRows;
            final double scale = numTrue == 0 || numTrue == numRows ?
                    1.0 :
                    FastMath.sqrt(baseline * (1.0 - baseline));
            setBaselineAndScale(baseline, scale);
        }
        @Override protected double[] getValuesAsOrderedDoubles() {
            // special case for booleans
            throw new GATKException("Method not needed for booleans");
        }
    }

    static public class ByteArrProperty extends Property {
        public byte[] values;

        ByteArrProperty(final String name, final int numValues) { super(name); values = new byte[numValues]; }
        ByteArrProperty(final String name) { this(name, DEFAULT_INITIAL_NUM_ALLOCATED_ROWS); }

        @Override public Property set(Object values) { this.values = (byte[]) values; return this; }
        @Override public JSONArray getJSON() {
            final JSONArray jsonArray = new JSONArray();
            IntStream.range(0, numRows).forEach(i -> jsonArray.add(values[i]));
            return jsonArray;
        }
        @Override public float getAsFloat(final int rowIndex, final int hyperIndex) { return values[rowIndex]; }
        @Override public int getAsInt(final int rowIndex, final int hyperIndex) { return values[rowIndex]; }
        @Override public int getAllocatedRows() { return values.length; }
        @Override public void setAllocatedRowsUnguarded(final int numRows) {
            values = Arrays.copyOf(values, numRows);
        }
        @Override protected void assignNextValue(final Object value) { values[numRows] = (byte)value; }
        @Override protected double[] getValuesAsOrderedDoubles() {
            return IntStream.range(0, numRows).mapToDouble(i -> (double)values[i]).sorted().toArray();
        }
    }

    static public class ByteMatProperty extends Property {
        public byte[][] values;

        ByteMatProperty(final String name, final int numValues) { super(name); values = new byte[numValues][]; }
        ByteMatProperty(final String name) { this(name, DEFAULT_INITIAL_NUM_ALLOCATED_ROWS); }

        @Override public Property set(Object values) { this.values = (byte[][]) values; return this; }
        @Override public JSONArray getJSON() {
            final JSONArray jsonArray = new JSONArray();
            Arrays.stream(values, 0, numRows).forEach(jsonArray::add);
            return jsonArray;
        }
        @Override public float getAsFloat(final int rowIndex, final int hyperIndex) {
            return values[rowIndex][hyperIndex];
        }
        @Override public int getAsInt(final int rowIndex, final int hyperIndex) {
            return values[rowIndex][hyperIndex];
        }
        @Override protected int[] getArrayOfDistinctNumColumns() {
            return IntStream.range(0, numRows).map(i -> values[i].length).distinct().toArray();
        }
        @Override public int getAllocatedRows() { return values.length; }
        @Override public void setAllocatedRowsUnguarded(final int numRows) {
            values = Arrays.copyOf(values, numRows);
        }
        @Override protected void assignNextValue(final Object value) { values[numRows] = (byte[])value; }
        @Override protected double[] getValuesAsOrderedDoubles() {
            return Arrays.stream(values, 0, numRows)
                    .flatMapToDouble(arr -> IntStream.range(0, arr.length).mapToDouble(i -> arr[i]))
                    .sorted()
                    .toArray();
        }
    }

    static public class ShortArrProperty extends Property {
        public short[] values;

        ShortArrProperty(final String name, final int numValues) { super(name); values = new short[numValues]; }
        ShortArrProperty(final String name) { this(name, DEFAULT_INITIAL_NUM_ALLOCATED_ROWS); }

        @Override public Property set(Object values) { this.values = (short[]) values; return this; }
        @Override public JSONArray getJSON() {
            final JSONArray jsonArray = new JSONArray();
            IntStream.range(0, numRows).forEach(i -> jsonArray.add(values[i]));
            return jsonArray;
        }
        @Override public float getAsFloat(final int rowIndex, final int hyperIndex) { return values[rowIndex]; }
        @Override public int getAsInt(final int rowIndex, final int hyperIndex) { return values[rowIndex]; }
        @Override public int getAllocatedRows() { return values.length; }
        @Override public void setAllocatedRowsUnguarded(final int numRows) {
            values = Arrays.copyOf(values, numRows);
        }
        @Override protected void assignNextValue(final Object value) { values[numRows] = (short)value; }
        @Override protected double[] getValuesAsOrderedDoubles() {
            return IntStream.range(0, numRows).mapToDouble(i -> (double)values[i]).sorted().toArray();
        }
    }

    static public class ShortMatProperty extends Property {
        public short[][] values;

        ShortMatProperty(final String name, final int numValues) { super(name); values = new short[numValues][]; }
        ShortMatProperty(final String name) { this(name, DEFAULT_INITIAL_NUM_ALLOCATED_ROWS); }

        @Override public Property set(Object values) { this.values = (short[][]) values; return this; }
        @Override public JSONArray getJSON() {
            final JSONArray jsonArray = new JSONArray();
            Arrays.stream(values, 0, numRows).forEach(jsonArray::add);
            return jsonArray;
        }
        @Override public float getAsFloat(final int rowIndex, final int hyperIndex) {
            return values[rowIndex][hyperIndex];
        }
        @Override public int getAsInt(final int rowIndex, final int hyperIndex) {
            return values[rowIndex][hyperIndex];
        }
        @Override protected int[] getArrayOfDistinctNumColumns() {
            return IntStream.range(0, numRows).map(i -> values[i].length).distinct().toArray();
        }
        @Override public int getAllocatedRows() { return values.length; }
        @Override public void setAllocatedRowsUnguarded(final int numRows) {
            values = Arrays.copyOf(values, numRows);
        }
        @Override protected void assignNextValue(final Object value) { values[numRows] = (short[])value; }
        @Override protected double[] getValuesAsOrderedDoubles() {
            return Arrays.stream(values, 0, numRows)
                    .flatMapToDouble(arr -> IntStream.range(0, arr.length).mapToDouble(i -> arr[i]))
                    .sorted()
                    .toArray();
        }
    }

    static public class IntArrProperty extends Property {
        public int[] values;

        IntArrProperty(final String name, final int numValues) { super(name); values = new int[numValues]; }
        IntArrProperty(final String name) { this(name, DEFAULT_INITIAL_NUM_ALLOCATED_ROWS); }

        @Override public Property set(Object values) { this.values = (int[]) values; return this; }
        @Override public JSONArray getJSON() {
            final JSONArray jsonArray = new JSONArray();
            Arrays.stream(values, 0, numRows).forEach(jsonArray::add);
            return jsonArray;
        }
        @Override public float getAsFloat(final int rowIndex, final int hyperIndex) { return values[rowIndex]; }
        @Override public int getAsInt(final int rowIndex, final int hyperIndex) { return values[rowIndex]; }
        @Override public int getAllocatedRows() { return values.length; }
        @Override public void setAllocatedRowsUnguarded(final int numRows) {
            values = Arrays.copyOf(values, numRows);
        }
        @Override protected void assignNextValue(final Object value) { values[numRows] = (int)value; }
        @Override protected double[] getValuesAsOrderedDoubles() {
            return Arrays.stream(values, 0, numRows).sorted().mapToDouble(x -> x).toArray();
        }
    }

    static public class IntMatProperty extends Property {
        public int[][] values;

        IntMatProperty(final String name, final int numValues) { super(name); values = new int[numValues][]; }
        IntMatProperty(final String name) { this(name, DEFAULT_INITIAL_NUM_ALLOCATED_ROWS); }

        @Override public Property set(Object values) { this.values = (int[][]) values; return this; }
        @Override public JSONArray getJSON() {
            final JSONArray jsonArray = new JSONArray();
            Arrays.stream(values, 0, numRows).forEach(jsonArray::add);
            return jsonArray;
        }
        @Override public float getAsFloat(final int rowIndex, final int hyperIndex) {
            return values[rowIndex][hyperIndex];
        }
        @Override public int getAsInt(final int rowIndex, final int hyperIndex) {
            return values[rowIndex][hyperIndex];
        }
        @Override protected int[] getArrayOfDistinctNumColumns() {
            return IntStream.range(0, numRows).map(i -> values[i].length).distinct().toArray();
        }
        @Override public int getAllocatedRows() { return values.length; }
        @Override public void setAllocatedRowsUnguarded(final int numRows) {
            values = Arrays.copyOf(values, numRows);
        }
        @Override protected void assignNextValue(final Object value) { values[numRows] = (int[])value; }
        @Override protected double[] getValuesAsOrderedDoubles() {
            return Arrays.stream(values, 0, numRows).flatMapToInt(Arrays::stream).sorted()
                .mapToDouble(x -> x).toArray();
        }
    }

    static public class LongArrProperty extends Property {
        public long[] values;

        LongArrProperty(final String name, final int numValues) { super(name); values = new long[numValues]; }
        LongArrProperty(final String name) { this(name, DEFAULT_INITIAL_NUM_ALLOCATED_ROWS); }

        @Override public Property set(Object values) { this.values = (long[]) values; return this; }
        @Override public JSONArray getJSON() {
            final JSONArray jsonArray = new JSONArray();
            Arrays.stream(values, 0, numRows).forEach(jsonArray::add);
            return jsonArray;
        }
        @Override public float getAsFloat(final int rowIndex, final int hyperIndex) { return values[rowIndex]; }
        @Override public long getAsLong(final int rowIndex, final int hyperIndex) { return values[rowIndex]; }
        @Override public int getAllocatedRows() { return values.length; }
        @Override public void setAllocatedRowsUnguarded(final int numRows) {
            values = Arrays.copyOf(values, numRows);
        }
        @Override protected void assignNextValue(final Object value) { values[numRows] = (long)value; }
        @Override protected double[] getValuesAsOrderedDoubles() {
            return Arrays.stream(values, 0, numRows).sorted().mapToDouble(x -> x).toArray();
        }
    }

    static public class LongMatProperty extends Property {
        public long[][] values;

        LongMatProperty(final String name, final int numValues) { super(name); values = new long[numValues][]; }
        LongMatProperty(final String name) { this(name, DEFAULT_INITIAL_NUM_ALLOCATED_ROWS); }

        @Override public Property set(Object values) { this.values = (long[][]) values; return this; }
        @Override public JSONArray getJSON() {
            final JSONArray jsonArray = new JSONArray();
            Arrays.stream(values, 0, numRows).forEach(jsonArray::add);
            return jsonArray;
        }
        @Override public float getAsFloat(final int rowIndex, final int hyperIndex) {
            return values[rowIndex][hyperIndex];
        }
        @Override public long getAsLong(final int rowIndex, final int hyperIndex) {
            return values[rowIndex][hyperIndex];
        }
        @Override protected int[] getArrayOfDistinctNumColumns() {
            return IntStream.range(0, numRows).map(i -> values[i].length).distinct().toArray();
        }
        @Override public int getAllocatedRows() { return values.length; }
        @Override public void setAllocatedRowsUnguarded(final int numRows) {
            values = Arrays.copyOf(values, numRows);
        }
        @Override protected void assignNextValue(final Object value) { values[numRows] = (long[])value; }
        @Override protected double[] getValuesAsOrderedDoubles() {
            return Arrays.stream(values, 0, numRows).flatMapToLong(Arrays::stream).sorted()
                .mapToDouble(x -> x).toArray();
        }
    }

    static public class FloatArrProperty extends Property {
        public float[] values;

        FloatArrProperty(final String name, final int numValues) { super(name); values = new float[numValues]; }
        FloatArrProperty(final String name) { this(name, DEFAULT_INITIAL_NUM_ALLOCATED_ROWS); }

        @Override public Property set(Object values) { this.values = (float[]) values; return this; }
        @Override public JSONArray getJSON() {
            final JSONArray jsonArray = new JSONArray();
            IntStream.range(0, numRows).forEach(i -> jsonArray.add(values[i]));
            return jsonArray;
        }
        @Override public float getAsFloat(final int rowIndex, final int hyperIndex) { return values[rowIndex]; }
        @Override public int getAllocatedRows() { return values.length; }
        @Override public void setAllocatedRowsUnguarded(final int numRows) {
            values = Arrays.copyOf(values, numRows);
        }
        @Override protected void assignNextValue(final Object value) { values[numRows] = (float)value; }
        @Override protected double[] getValuesAsOrderedDoubles() {
            return IntStream.range(0, numRows).mapToDouble(i -> (double)values[i]).sorted().toArray();
        }
    }

    static public class FloatMatProperty extends Property {
        public float[][] values;

        FloatMatProperty(final String name, final int numValues) { super(name); values = new float[numValues][]; }
        FloatMatProperty(final String name) { this(name, DEFAULT_INITIAL_NUM_ALLOCATED_ROWS); }

        @Override public Property set(Object values) { this.values = (float[][]) values; return this; }
        @Override public JSONArray getJSON() {
            final JSONArray jsonArray = new JSONArray();
            IntStream.range(0, numRows).forEach(i -> jsonArray.add(values[i]));
            return jsonArray;
        }
        @Override public float getAsFloat(final int rowIndex, final int hyperIndex) {
            return values[rowIndex][hyperIndex];
        }
        @Override protected int[] getArrayOfDistinctNumColumns() {
            return IntStream.range(0, numRows).map(i -> values[i].length).distinct().toArray();
        }
        @Override public int getAllocatedRows() { return values.length; }
        @Override public void setAllocatedRowsUnguarded(final int numRows) {
            values = Arrays.copyOf(values, numRows);
        }
        @Override protected void assignNextValue(final Object value) { values[numRows] = (float[])value; }
        @Override protected double[] getValuesAsOrderedDoubles() {
            return Arrays.stream(values, 0, numRows)
                .flatMapToDouble(arr -> IntStream.range(0, arr.length).mapToDouble(i -> arr[i]))
                .sorted()
                .toArray();
        }
    }

    static public class DoubleArrProperty extends Property {
        public double[] values;

        DoubleArrProperty(final String name, final int numValues) { super(name); values = new double[numValues]; }
        DoubleArrProperty(final String name) { this(name, DEFAULT_INITIAL_NUM_ALLOCATED_ROWS); }

        @Override public Property set(Object values) { this.values = (double[]) values; return this; }
        @Override public JSONArray getJSON() {
            final JSONArray jsonArray = new JSONArray();
            Arrays.stream(values, 0, numRows).forEach(jsonArray::add);
            return jsonArray;
        }
        public double getAsDouble(final int rowIndex) { return values[rowIndex]; }
        @Override public float getAsFloat(final int rowIndex, final int hyperIndex) { return (float)values[rowIndex]; }
        @Override public int getAllocatedRows() { return values.length; }
        @Override public void setAllocatedRowsUnguarded(final int numRows) {
            values = Arrays.copyOf(values, numRows);
        }
        @Override protected void assignNextValue(final Object value) { values[numRows] = (double)value; }
        @Override protected double[] getValuesAsOrderedDoubles() {
            return Arrays.stream(values, 0, numRows).sorted().toArray();
        }
    }

    static public class DoubleMatProperty extends Property {
        public double[][] values;

        DoubleMatProperty(final String name, final int numValues) { super(name); values = new double[numValues][]; }
        DoubleMatProperty(final String name) { this(name, DEFAULT_INITIAL_NUM_ALLOCATED_ROWS); }

        @Override public Property set(Object values) { this.values = (double[][]) values; return this; }
        @Override public JSONArray getJSON() {
            final JSONArray jsonArray = new JSONArray();
            Arrays.stream(values, 0, numRows).forEach(jsonArray::add);
            return jsonArray;
        }
        @Override public float getAsFloat(final int rowIndex, final int hyperIndex) {
            return (float)values[rowIndex][hyperIndex];
        }
        @Override protected int[] getArrayOfDistinctNumColumns() {
            return IntStream.range(0, numRows).map(i -> values[i].length).distinct().toArray();
        }
        @Override public int getAllocatedRows() { return values.length; }
        @Override public void setAllocatedRowsUnguarded(final int numRows) {
            values = Arrays.copyOf(values, numRows);
        }
        @Override protected void assignNextValue(final Object value) { values[numRows] = (double[])value; }
        @Override protected double[] getValuesAsOrderedDoubles() {
            return Arrays.stream(values, 0, numRows).flatMapToDouble(Arrays::stream).sorted().toArray();
        }
    }


    static public class StringArrProperty extends Property {
        private final List<String> indexToString = new ArrayList<>();
        private final Map<String, Integer> stringToIndex = new HashMap<>();
        private int[] ordinalEncoding;

        StringArrProperty(final String name, final int numValues) {
            super(name);
            ordinalEncoding = new int[numValues];
        }
        StringArrProperty(final String name) { this(name, DEFAULT_INITIAL_NUM_ALLOCATED_ROWS); }

        private int encode(final String value) {
            final Integer index = stringToIndex.getOrDefault(value, null);
            if(index == null) {
                final int newIndex = indexToString.size();
                indexToString.add(value);
                stringToIndex.put(value, newIndex);
                return newIndex;
            } else {
                return index;
            }
        }

        @Override public Property set(Object values) {
            indexToString.clear();
            stringToIndex.clear();
            ordinalEncoding = Arrays.stream((String[]) values).mapToInt(this::encode).toArray();
            return this;
        }
        @Override public JSONArray getJSON() {
            final JSONArray jsonArray = new JSONArray();
            Arrays.stream(ordinalEncoding, 0, numRows).mapToObj(indexToString::get).forEach(jsonArray::add);
            return jsonArray;
        }
        public String getAsString(final int rowIndex) {
            return indexToString.get(ordinalEncoding[rowIndex]);
        }
        @Override public float getAsFloat(final int rowIndex, final int hyperIndex) {
            throw new IllegalArgumentException("String properties cannot be converted to Float, they must be one-hot-encoded");
        }
        @Override public int getAllocatedRows() { return ordinalEncoding.length; }
        @Override public void setAllocatedRowsUnguarded(final int numRows) {
            ordinalEncoding = Arrays.copyOf(ordinalEncoding, numRows);
        }
        @Override protected void assignNextValue(final Object value) {
            ordinalEncoding[numRows] = encode((String)value);
        }
        @Override protected double[] getValuesAsOrderedDoubles() {
            throw new IllegalArgumentException("String properties cannot be converted to double, they must be one-hot-encoded");
        }
        @Override public boolean isNumeric() { return false; }
        @Override public List<String> getAllLabels() {
            return indexToString.stream()
                    .distinct()
                    .sorted()
                    .collect(Collectors.toList());
        }
        @Override public void oneHot(final List<String> allLabels, final PropertiesTable propertiesTable) {
            final int numProperties = allLabels.size();
            final List<Property> properties = IntStream.range(0, numProperties)
                    .mapToObj(i ->
                            (BooleanArrProperty)propertiesTable.getOrCreateProperty(
                                    name + "=" + allLabels.get(i),
                                    PropertyClass.BooleanArrProperty, numRows
                            )
                    )
                    .collect(Collectors.toList());
            properties.forEach(p -> p.numRows = 0);
            for(int row = 0; row < numRows; ++row) {
                final String value = indexToString.get(ordinalEncoding[row]);
                for(int propertyIndex = 0; propertyIndex < numProperties; ++propertyIndex) {
                    properties.get(propertyIndex).append(value.equals(allLabels.get(propertyIndex)));
                }
            }
        }
    }

    static public class StringMatProperty extends Property {
        private final List<String> indexToString = new ArrayList<>();
        private final Map<String, Integer> stringToIndex = new HashMap<>();
        private int[][] ordinalEncoding;

        StringMatProperty(final String name, final int numValues) {
            super(name);
            ordinalEncoding = new int[numValues][];
        }
        StringMatProperty(final String name) { this(name, DEFAULT_INITIAL_NUM_ALLOCATED_ROWS); }

        private int encode(final String value) {
            final Integer index = stringToIndex.getOrDefault(value, null);
            if(index == null) {
                final int newIndex = indexToString.size();
                indexToString.add(value);
                stringToIndex.put(value, newIndex);
                return newIndex;
            } else {
                return index;
            }
        }

        private int[] encodeRow(final String[] row) {
            return Arrays.stream(row).mapToInt(this::encode).toArray();
        }

        @Override public Property set(Object values) {
            indexToString.clear();
            stringToIndex.clear();
            ordinalEncoding = Arrays.stream((String[][]) values).map(this::encodeRow).toArray(int[][]::new);
            return this;
        }
        @Override public JSONArray getJSON() {
            final JSONArray jsonArray = new JSONArray();
            Arrays.stream(ordinalEncoding, 0, numRows)
                .map(row -> Arrays.stream(row).mapToObj(indexToString::get).toArray(String[]::new))
                .forEach(jsonArray::add);
            return jsonArray;
        }
        @Override public float getAsFloat(final int rowIndex, final int hyperIndex) {
            throw new IllegalArgumentException("String properties cannot be converted to Float, they must be one-hot-encoded");
        }
        @Override protected int[] getArrayOfDistinctNumColumns() {
            return IntStream.range(0, numRows).map(i -> ordinalEncoding[i].length).distinct().toArray();
        }
        @Override public int getAllocatedRows() { return ordinalEncoding.length; }
        @Override public void setAllocatedRowsUnguarded(final int numRows) {
            ordinalEncoding = Arrays.copyOf(ordinalEncoding, numRows);
        }
        @Override protected void assignNextValue(final Object value) {
            ordinalEncoding[numRows] = encodeRow((String[])value);
        }
        @Override protected double[] getValuesAsOrderedDoubles() {
            throw new IllegalArgumentException("String properties cannot be converted to double, they must be one-hot-encoded");
        }
        @Override public boolean isNumeric() { return false; }
        @Override public List<String> getAllLabels() {
            return indexToString.stream()
                    .distinct()
                    .sorted()
                    .collect(Collectors.toList());
        }
        @Override public void oneHot(final List<String> allLabels, final PropertiesTable propertiesTable) {
            final int numColumns = numRows > 0 ? getNumColumns() : 0;
            final int numProperties = allLabels.size();
            final List<BooleanMatProperty> properties = IntStream.range(0, numProperties)
                    .mapToObj(i ->
                            (BooleanMatProperty)propertiesTable.getOrCreateProperty(
                                    name + "=" + allLabels.get(i),
                                    PropertyClass.BooleanMatProperty, 0
                            ).set(new boolean[numRows][numColumns])
                    )
                    .collect(Collectors.toList());
            properties.forEach(p -> p.numRows = 0);
            for(int row = 0; row < numRows; ++row) {
                final int[] rowIndices = ordinalEncoding[row];
                for(int propertyIndex = 0; propertyIndex < numProperties; ++propertyIndex) {
                    final String propertyLabel = allLabels.get(propertyIndex);;
                    final BooleanMatProperty booleanMatProperty = properties.get(propertyIndex);
                    for(int col = 0; col < numColumns; ++col) {
                        final String value = indexToString.get(rowIndices[col]);
                        booleanMatProperty.assignColumnValue(col, value.equals(propertyLabel));
                    }
                    ++booleanMatProperty.numRows;
                }
            }
        }
    }

    @SuppressWarnings("unchecked")
    static public class StringSetArrProperty extends Property {
        /*
        Use ordinal encoding for each unique *combination* of strings during load-in: it achieves substantial
        compression over storing the Set<String>, and is much easier to manage than attempting to one-hot encode on the
        fly while additional categories are being added.
         */
        private final Map<Set<String>, Integer> setsToIndex = new HashMap<>();
        private final List<Set<String>> indexToSets = new ArrayList<>();
        private int[] ordinalEncoding;

        StringSetArrProperty(final String name, final int numValues) {
            super(name);
            ordinalEncoding = new int[numValues];
        }
        StringSetArrProperty(final String name) { this(name, DEFAULT_INITIAL_NUM_ALLOCATED_ROWS); }

        private int encode(final Set<String> value) {
            final Integer index = setsToIndex.getOrDefault(value, null);
            if(index == null) {
                final int newIndex = indexToSets.size();
                indexToSets.add(value);
                setsToIndex.put(value, newIndex);
                return newIndex;
            } else {
                return index;
            }
        }

        @Override protected void assignNextValue(final Object value) {
            ordinalEncoding[numRows] = encode((Set<String>)value);
        }

        @Override public Property set(Object values) {
            setsToIndex.clear();
            indexToSets.clear();
            ordinalEncoding = Arrays.stream((Set<String>[]) values).mapToInt(this::encode).toArray();
            return this;
        }

        @Override public JSONArray getJSON() {
            final JSONArray jsonArray = new JSONArray();
            Arrays.stream(ordinalEncoding, 0, numRows).mapToObj(indexToSets::get).forEach(jsonArray::add);
            return jsonArray;
        }
        @Override public float getAsFloat(final int rowIndex, final int hyperIndex) {
            throw new IllegalArgumentException("String properties cannot be converted to Float, they must be one-hot-encoded");
        }
        @Override public int getAllocatedRows() { return ordinalEncoding.length; }
        @Override public void setAllocatedRowsUnguarded(final int numRows) {
            ordinalEncoding = Arrays.copyOf(ordinalEncoding, numRows);
        }
        @Override protected double[] getValuesAsOrderedDoubles() {
            throw new IllegalArgumentException("String properties cannot be converted to double, they must be one-hot-encoded");
        }
        @Override public boolean isNumeric() { return false; }
        @Override public List<String> getAllLabels() {
            return indexToSets.stream()
                    .flatMap(Collection::stream)
                    .distinct()
                    .sorted()
                    .collect(Collectors.toList());
        }
        @Override public void oneHot(final List<String> allLabels, final PropertiesTable propertiesTable) {
            final int numProperties = allLabels.size();
            final List<Property> properties = IntStream.range(0, numProperties)
                    .mapToObj(i ->
                            (BooleanArrProperty)propertiesTable.getOrCreateProperty(
                                    name + "=" + allLabels.get(i),
                                    PropertyClass.BooleanArrProperty, numRows
                            )
                    )
                    .collect(Collectors.toList());
            properties.forEach(p -> p.numRows = 0);
            for(int row = 0; row < numRows; ++row) {
                final Set<String> value = indexToSets.get(ordinalEncoding[row]);
                for(int propertyIndex = 0; propertyIndex < numProperties; ++propertyIndex) {
                    properties.get(propertyIndex).append(
                        value.contains(allLabels.get(propertyIndex))
                    );
                }
            }
        }
    }

    @SuppressWarnings("unchecked")
    static public class StringSetMatProperty extends Property {
        /*
        Use ordinal encoding for each unique *combination* of strings during load-in: it achieves substantial
        compression over storing the Set<String>, and is much easier to manage than attempting to one-hot encode on the
        fly while additional categories are being added.
         */
        private final Map<Set<String>, Integer> setsToIndex = new HashMap<>();
        private final List<Set<String>> indexToSets = new ArrayList<>();
        private int[][] ordinalEncoding;

        StringSetMatProperty(final String name, final int numValues) {
            super(name);
            ordinalEncoding = new int[numValues][];
        }
        StringSetMatProperty(final String name) { this(name, DEFAULT_INITIAL_NUM_ALLOCATED_ROWS); }

        private int encode(final Set<String> value) {
            final Integer index = setsToIndex.getOrDefault(value, null);
            if(index == null) {
                final int newIndex = indexToSets.size();
                indexToSets.add(value);
                setsToIndex.put(value, newIndex);
                return newIndex;
            } else {
                return index;
            }
        }

        private int[] encodeRow(final Set<String>[] row) {
            return Arrays.stream(row).mapToInt(this::encode).toArray();
        }

        @Override protected void assignNextValue(final Object value) {
            ordinalEncoding[numRows] = encodeRow((Set<String>[])value);
        }

        @Override public Property set(Object values) {
            setsToIndex.clear();
            indexToSets.clear();
            ordinalEncoding = Arrays.stream((Set<String>[][]) values).map(this::encodeRow).toArray(int[][]::new);
            return this;
        }
        @Override public JSONArray getJSON() {
            final JSONArray jsonArray = new JSONArray();
            Arrays.stream(ordinalEncoding, 0, numRows)
                .map(row -> Arrays.stream(row).mapToObj(indexToSets::get).toArray())
                .forEach(jsonArray::add);
            return jsonArray;
        }
        @Override public float getAsFloat(final int rowIndex, final int hyperIndex) {
            throw new IllegalArgumentException("String properties cannot be converted to Float, they must be one-hot-encoded");
        }
        @Override protected int[] getArrayOfDistinctNumColumns() {
            return IntStream.range(0, numRows).map(i -> ordinalEncoding[i].length).distinct().toArray();
        }
        @Override public int getAllocatedRows() { return ordinalEncoding.length; }
        @Override public void setAllocatedRowsUnguarded(final int numRows) {
            ordinalEncoding = Arrays.copyOf(ordinalEncoding, numRows);
        }

        @Override protected double[] getValuesAsOrderedDoubles() {
            throw new IllegalArgumentException("String properties cannot be converted to double, they must be one-hot-encoded");
        }
        @Override public boolean isNumeric() { return false; }
        @Override public List<String> getAllLabels() {
            return indexToSets.stream()
                    .flatMap(Collection::stream)
                    .distinct()
                    .sorted()
                    .collect(Collectors.toList());
        }
        @Override public void oneHot(final List<String> allLabels, final PropertiesTable propertiesTable) {
            final int numProperties = allLabels.size();
            final int numColumns = numRows > 0 ? getNumColumns() : 0;
                final List<BooleanMatProperty> properties = IntStream.range(0, numProperties)
                        .mapToObj(i ->
                                (BooleanMatProperty)propertiesTable.getOrCreateProperty(
                                        name + "=" + allLabels.get(i),
                                        PropertyClass.BooleanMatProperty, numRows
                                ).set(new boolean[numRows][numColumns])
                        )
                        .collect(Collectors.toList());
                properties.forEach(p -> p.numRows = 0);
                for(int row = 0; row < numRows; ++row) {
                    final int[] rowIndices = ordinalEncoding[row];
                    for(int propertyIndex = 0; propertyIndex < numProperties; ++propertyIndex) {
                        final BooleanMatProperty booleanMatProperty = properties.get(propertyIndex);
                        final String propertyLabel = allLabels.get(propertyIndex);;
                        for(int col = 0; col < numColumns; ++col) {
                            final Set<String> value = indexToSets.get(rowIndices[col]);
                            booleanMatProperty.assignColumnValue(col, value.contains(propertyLabel));
                        }
                        ++booleanMatProperty.numRows;
                    }
                }
        }
    }
}