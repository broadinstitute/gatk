package org.broadinstitute.hellbender.tools.funcotator.dataSources.XSV;

import htsjdk.tribble.Feature;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.funcotator.DataSourceFuncotationFactory;
import org.broadinstitute.hellbender.tools.funcotator.Funcotation;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation;

import java.util.LinkedHashSet;
import java.util.List;

/**
 *  Factory for creating {@link XSVFuncotation}s by handling `Separated Value` files with arbitrary delimiters
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

    //==================================================================================================================
    // Private Members:

    //==================================================================================================================
    // Constructors:

    //==================================================================================================================
    // Override Methods:

    @Override
    public String getName() {
        return null;
    }

    @Override
    public LinkedHashSet<String> getSupportedFuncotationFields() {
        return null;
    }

    @Override
    public List<Funcotation> createFuncotations(final VariantContext variant, final ReferenceContext referenceContext, final List<Feature> featureList) {
        return null;
    }

    @Override
    public List<Funcotation> createFuncotations(final VariantContext variant, final ReferenceContext referenceContext, final List<Feature> featureList, final List<GencodeFuncotation> gencodeFuncotations) {
        return null;
    }


    //==================================================================================================================
    // Static Methods:

    //==================================================================================================================
    // Instance Methods:

    //==================================================================================================================
    // Helper Data Types:

}
