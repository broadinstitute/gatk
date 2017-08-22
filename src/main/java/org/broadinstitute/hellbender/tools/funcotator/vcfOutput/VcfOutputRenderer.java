package org.broadinstitute.hellbender.tools.funcotator.vcfOutput;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.broadinstitute.hellbender.tools.funcotator.DataSourceFuncotationFactory;
import org.broadinstitute.hellbender.tools.funcotator.Funcotation;
import org.broadinstitute.hellbender.tools.funcotator.Funcotator;
import org.broadinstitute.hellbender.tools.funcotator.OutputRenderer;

import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * A Funcotator output renderer for writing to VCF files.
 * Created by jonn on 8/30/17.
 */
public class VcfOutputRenderer extends OutputRenderer {

    /**
     * The name of the field inside the VCF INFO field in which to put {@link Funcotator} results.
     */
    public static final String FUNCOTATOR_VCF_FIELD_NAME = "FUNCOTATION";

    //==================================================================================================================

    private final VariantContextWriter vcfWriter;
    private final VCFHeader existingHeader;
    private final List<DataSourceFuncotationFactory> dataSourceFactories;

    //==================================================================================================================

    public VcfOutputRenderer(final VCFHeader existingHeader,
                             final VariantContextWriter vcfWriter,
                             final List<DataSourceFuncotationFactory> dataSources) {
        this.vcfWriter = vcfWriter;
        this.existingHeader = existingHeader;
        this.dataSourceFactories = dataSources;
    }

    //==================================================================================================================

    @Override
    public void open() {
        final VCFHeader newHeader = createVCFHeader(existingHeader, dataSourceFactories);
        vcfWriter.writeHeader(newHeader);
    }

    @Override
    public void close() {
        vcfWriter.close();
    }

    @Override
    public void write(final VariantContext variant, final List<Funcotation> funcotations) {

        // Create a new variant context builder:
        final VariantContextBuilder variantContextBuilder = new VariantContextBuilder(variant);

        final StringBuilder funcotatorAnnotationStringBuilder = new StringBuilder();

        // Get the old VCF Annotation field and append the new information to it:
        final Object existingAnnotation = variant.getAttribute(FUNCOTATOR_VCF_FIELD_NAME, null);
        if ( existingAnnotation != null) {
            funcotatorAnnotationStringBuilder.append( existingAnnotation.toString() );
            funcotatorAnnotationStringBuilder.append( ',' );
        }

        funcotatorAnnotationStringBuilder.append(
                funcotations.stream().map(Funcotation::serializeToVcfString).collect(Collectors.joining(","))
        );

        // Add our new annotation and render the VariantContext:
        variantContextBuilder.attribute(FUNCOTATOR_VCF_FIELD_NAME, funcotatorAnnotationStringBuilder.toString());

        vcfWriter.add( variantContextBuilder.make() );
    }

    //==================================================================================================================

    /**
     * Create a header for a VCF file.
     * Uses {@link VcfOutputRenderer#dataSourceFactories} to get a list of fields to report producing (preserving their order).
     * @param existingHeader An existing {@link VCFHeader} from which to replicate fields in this new file.
     * @param dataSourceFactories A {@link List} of {@link DataSourceFuncotationFactory} objects from which to pull field names.
     * @return The {@link VCFHeader} object with relevant information for {@link Funcotator}.
     */
    private static VCFHeader createVCFHeader(final VCFHeader existingHeader,
                                             final List<DataSourceFuncotationFactory> dataSourceFactories) {
        final Set<VCFHeaderLine> headerLines = new HashSet<>();

        // Add all lines of our existing VCF header:
        headerLines.addAll( existingHeader.getMetaDataInInputOrder() );

        // Add in the line about Funcotations:
        headerLines.add(new VCFInfoHeaderLine(FUNCOTATOR_VCF_FIELD_NAME, VCFHeaderLineCount.A,
                VCFHeaderLineType.String, "Functional annotation from the Funcotator tool.  Funcotation fields are: "
                + getDataSourceFieldNamesForHeader(dataSourceFactories))
        );

        return new VCFHeader(headerLines);
    }

    /**
     * Creates a {@link String} containing the field names from our {@link VcfOutputRenderer#dataSourceFactories} suitable for putting in the VCF header.
     * @param dataSourceFactories A {@link List} of {@link DataSourceFuncotationFactory} objects from which to pull field names.
     * @return A {@link String} containing the field names from our {@link VcfOutputRenderer#dataSourceFactories} suitable for putting in the VCF header.
     */
    private static String getDataSourceFieldNamesForHeader(final List<DataSourceFuncotationFactory> dataSourceFactories) {
        return dataSourceFactories.stream()
                        .map(DataSourceFuncotationFactory::getSupportedFuncotationFields)
                        .flatMap(List::stream)
                        .map(Object::toString).collect(Collectors.joining(";"));
    }
}
