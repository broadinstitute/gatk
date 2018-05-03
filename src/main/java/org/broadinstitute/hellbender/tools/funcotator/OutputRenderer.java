package org.broadinstitute.hellbender.tools.funcotator;

import htsjdk.variant.variantcontext.VariantContext;

import java.util.LinkedHashMap;
import java.util.List;
import java.util.stream.Collectors;

/**
 * An abstract class to allow for writing output for the Funcotator.
 * Used to output Funcotations to a location of the user's choice.
 * For example, writing out to a VCF file (e.g. {@link org.broadinstitute.hellbender.tools.funcotator.vcfOutput.VcfOutputRenderer}).
 * It is also possible to write to a stream, or other connection (e.g. a UDP port)
 * at the user's discretion.
 * Created by jonn on 8/30/17.
 */
public abstract class OutputRenderer implements AutoCloseable {

    //==================================================================================================================
    /**
     * {@link LinkedHashMap} of manually specified annotations to add to each output in addition to annotations provided
     * to {@link OutputRenderer#write(VariantContext, FuncotationMap)}.
     */
    protected LinkedHashMap<String, String> manualAnnotations;

    /**
     * {@link String} representation of {@link OutputRenderer#manualAnnotations} serialized to the output format of this {@link OutputRenderer}.
     */
    protected String manualAnnotationSerializedString;

    /**
     * {@link List} of the {@link DataSourceFuncotationFactory} objects that are being used in this run of {@link Funcotator}.
     */
    protected List<DataSourceFuncotationFactory> dataSourceFactories;

    //==================================================================================================================

    /**
     * @return A {@link String} containing information about the data sources that are used to create the {@link Funcotation}s by this {@link OutputRenderer}.
     */
    public String getDataSourceInfoString() {

        return dataSourceFactories.stream()
                .map(DataSourceFuncotationFactory::getInfoString)
                .collect(Collectors.joining(" | "));
    }

    /**
     * Close the {@link OutputRenderer}.
     */
    public abstract void close();

    /**
     * Write the given {@code variant} and {@code txToFuncotationMap} to the output file.
     * @param variant {@link VariantContext} to write to the file.
     * @param txToFuncotationMap {@link FuncotationMap} to add to the given {@code variant} on output.
     */
    public abstract void write(final VariantContext variant, final FuncotationMap txToFuncotationMap);
}
