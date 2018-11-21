package org.broadinstitute.hellbender.tools.walkers.varianteval.util;

import org.broadinstitute.hellbender.exceptions.GATKException;

import java.lang.annotation.Annotation;
import java.lang.reflect.Field;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;


/**
 * @author aaron
 *         <p/>
 *         Class AnalysisModuleScanner
 *         <p/>
 *         Given an analysis, find the annotated fields and methods.  Given this module and
 *         the object, a Mashalling object can serialize or deserialize a analysis module.
 */
public class AnalysisModuleScanner {
    final private static Map<String, Annotation[]> annotationCache = new HashMap<String, Annotation[]>();

    // what we extracted from the class
    private Map<Field, DataPoint> datums = new LinkedHashMap<Field, DataPoint>();   // the data we've discovered
    private Analysis analysis;  // the analysis annotation
    
    private Field moltenField = null;
    private Molten moltenAnnotation = null;
    
    // private storage of the class type
    private final Class<?> cls;

    /**
     * create a report scanner from the passed in class
     * @param cls the target class, annotated with the @Analysis annotation
     */
    public AnalysisModuleScanner(Class<?> cls) {
        this.cls = cls;
        scan(); // scan the passed in class
    }

    /**
     * create a report scanner from the passed in class
     * @param obj the target object, annotated with the @Analysis annotation
     */
    public AnalysisModuleScanner(Object obj) {
        this.cls = obj.getClass();
        scan(); // scan the passed in class
    }

    /** scan the class and find all appropriate fields and tables */
    public void scan() {
        if (cls == null || !cls.isAnnotationPresent(Analysis.class))
            throw new GATKException("The class passed in cannot be null, " + "" +
                                     "and must contain the @Analysis annotation, class " + cls + " was the input");

        // get the annotation off of the class
        analysis = cls.getAnnotation(Analysis.class);
        scanFields();
    }

    /**
     * scan the fields of the class, extracting parameters and table annotations and their associated fields
     */
    private void scanFields() {
        // get the fields from the class, and extract
        for ( Class<?> superCls = cls; superCls != null; superCls=superCls.getSuperclass() ) {
            for (Field f : superCls.getDeclaredFields()) {
                for (Annotation annotation : getAnnotations(f)) {
                    if (annotation.annotationType().equals(DataPoint.class))
                        datums.put(f,(DataPoint) annotation);
                    if ( annotation.annotationType().equals(Molten.class)) {
                        if ( hasMoltenField() )
                            throw new GATKException("Analysis " + analysis.name() + " has multiple @Molten fields, which is forbidden");
                        moltenField = f;
                        moltenAnnotation = (Molten)annotation;
                    }
                }
            }
        }
        
        if ( hasMoltenField() ) {
            if ( datums.size() > 0 )
                throw new GATKException("Analysis " + analysis.name() + " has an @Molten field as well as @DataPoint fields, which is forbidden");
        }
    }

    public Field getMoltenField() {
        return moltenField;
    }

    public Molten getMoltenAnnotation() {
        return moltenAnnotation;
    }

    public boolean hasMoltenField() {
        return getMoltenField() != null;
    }

    private Annotation[] getAnnotations(final Field field) {
        final String fieldName = field.toString();
        Annotation[] annotations = annotationCache.get(fieldName);
        if ( annotations == null ) {
            annotations = field.getAnnotations();
            annotationCache.put(fieldName, annotations);
        }
        return annotations;
    }

    /**
     *
     * @return a map of the datum annotations found
     */
    public Map<Field, DataPoint> getData() {
        return datums;
    }

    /**
     *
     * @return the analysis annotation found
     */
    public Analysis getAnalysis() {
        return analysis;
    }

    public Class<?> getModuleClass() {
        return cls;
    }
}
