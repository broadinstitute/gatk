package org.broadinstitute.hellbender.tools.walkers.featuremapping;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.exceptions.GATKException;

import java.util.LinkedHashMap;
import java.util.Map;

public abstract class SnvqrFeature {

    private String name;
    private Map<Float, Object> dict;      // value -> category
    private Map<Object, Float> rDict;      // category -> value

    public SnvqrFeature(final String name) {
        this.name = name;
    }

    public float getValue(final MappedFeature mappedFeature, final byte nonCalledBase, SAMFileHeader header) {
        Object v = getObjectValue(mappedFeature, nonCalledBase, header);
        if ( rDict == null ) {
            if ( v instanceof Number ) {
                return ((Number)v).floatValue();
            } else if ( v instanceof Boolean ) {
                return ((Boolean)v).booleanValue() ? 1.0f : 0.0f;
            } else {
                return Float.parseFloat(v.toString());
            }
        } else {
            if ( !rDict.containsKey(v) ) {
                throw new GATKException("snvqrFeature " + name + " can not translate given value: " + v);
            } else {
                return rDict.get(v);
            }
        }
    }

    public Object getObjectValue(final MappedFeature mappedFeature) {
        if ( rDict == null ) {
            return 0.0f;
        } else {
            return rDict.keySet().iterator().next();
        }
    }

    public Object getObjectValue(final MappedFeature mappedFeature, final byte nonCalledBase) {
        return getObjectValue(mappedFeature);
    }

    public Object getObjectValue(final MappedFeature mappedFeature, final byte nonCalledBase, final SAMFileHeader header) {
        return getObjectValue(mappedFeature, nonCalledBase);
    }

    public void setCategoryDictionary(final Map<String, Object> dict) {
        this.rDict = new LinkedHashMap<>();
        this.dict = new LinkedHashMap<>();
        for ( Map.Entry<String, Object> entry : dict.entrySet() ) {
            this.dict.put(Float.parseFloat(entry.getKey()), entry.getValue());
            this.rDict.put(entry.getValue(), Float.parseFloat(entry.getKey()));

        }
    }

    public Object getCategoryValue(float value) {
        if ( dict == null ) {
            return Float.toString(value);
        } else {
            return dict.get(value);
        }
    }

    void setName(String name) {
        this.name = name;
    }

    public boolean isNonCalledBasedIndependent() {
        // if independent of the mapped feature, then it is definitely independent of the non-called-based
        return isMappedFeatureIndependent();
    }

    public boolean isMappedFeatureIndependent() {
        return false;
    }

    public String getName() {
        return this.name;
    }
}
