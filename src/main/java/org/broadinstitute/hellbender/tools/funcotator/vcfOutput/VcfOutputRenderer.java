package org.broadinstitute.hellbender.tools.funcotator.vcfOutput;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.broadinstitute.hellbender.tools.funcotator.DataSourceFuncotationFactory;
import org.broadinstitute.hellbender.tools.funcotator.Funcotation;
import org.broadinstitute.hellbender.tools.funcotator.Funcotator;
import org.broadinstitute.hellbender.tools.funcotator.OutputRenderer;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
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

    /**
     * The delimiter in the header for the list of fields that each funcotation has.
     */
    public static final String HEADER_LISTED_FIELD_DELIMITER = "|";

    /**
     * The delimiter for each field within the VCF Funcotation field.
     */
    public static final String FIELD_DELIMITER = "|";

    /**
     * The delimiter for the `Other Transcript` field within the Funcotation annotation in the VCF.
     */
    public static final String OTHER_TRANSCRIPT_DELIMITER = ";";

    //==================================================================================================================

    private final VariantContextWriter vcfWriter;
    private final VCFHeader existingHeader;

    private final LinkedHashSet<VCFHeaderLine> defaultToolVcfHeaderLines;

    //==================================================================================================================

    public VcfOutputRenderer(final VCFHeader existingHeader,
                             final VariantContextWriter vcfWriter,
                             final List<DataSourceFuncotationFactory> dataSources) {
        this(existingHeader, vcfWriter, dataSources, new LinkedHashMap<>(), new LinkedHashMap<>());
    }

    public VcfOutputRenderer(final VCFHeader existingHeader,
                             final VariantContextWriter vcfWriter,
                             final List<DataSourceFuncotationFactory> dataSources,
                             final LinkedHashMap<String, String> unaccountedForDefaultAnnotations,
                             final LinkedHashMap<String, String> unaccountedForOverrideAnnotations) {
        this(existingHeader, vcfWriter, dataSources, unaccountedForDefaultAnnotations, unaccountedForOverrideAnnotations, new LinkedHashSet<>());
    }

    public VcfOutputRenderer(final VCFHeader existingHeader,
                             final VariantContextWriter vcfWriter,
                             final List<DataSourceFuncotationFactory> dataSources,
                             final LinkedHashMap<String, String> unaccountedForDefaultAnnotations,
                             final LinkedHashMap<String, String> unaccountedForOverrideAnnotations,
                             final Set<VCFHeaderLine> defaultToolVcfHeaderLines) {

        this.vcfWriter = vcfWriter;
        this.existingHeader = existingHeader;
        this.dataSourceFactories = dataSources;

        // Merge the annotations into our manualAnnotations:
        manualAnnotations = new LinkedHashMap<>();
        manualAnnotations.putAll(unaccountedForDefaultAnnotations);
        manualAnnotations.putAll(unaccountedForOverrideAnnotations);

        // Get our default tool VCF header lines:
        this.defaultToolVcfHeaderLines = new LinkedHashSet<>(defaultToolVcfHeaderLines);

        // Cache the manual annotation string so we can pass it easily into any Funcotations:
        manualAnnotationSerializedString = (manualAnnotations.size() != 0 ? String.join( FIELD_DELIMITER, manualAnnotations.values() ) + FIELD_DELIMITER : "");
    }

    //==================================================================================================================

    @Override
    public void open() {
        final VCFHeader newHeader = createVCFHeader();
        vcfWriter.writeHeader(newHeader);
    }

    @Override
    public void close() {
        vcfWriter.close();
    }

    @Override
    public void write(final VariantContext variant, final List<Funcotation> funcotations) {

        // Create a new variant context builder:
        final VariantContextBuilder variantContextOutputBuilder = new VariantContextBuilder(variant);

        final StringBuilder funcotatorAnnotationStringBuilder = new StringBuilder();

        // Get the old VCF Annotation field and append the new information to it:
        final Object existingAnnotation = variant.getAttribute(FUNCOTATOR_VCF_FIELD_NAME, null);
        final List<String> existingAlleleAnnotations;
        if ( existingAnnotation != null) {
            existingAlleleAnnotations = Utils.split(existingAnnotation.toString(), ',');
        }
        else {
            existingAlleleAnnotations = Collections.emptyList();
        }

        // Go through each allele and add it to the writer separately:
        final List<Allele> alternateAlleles = variant.getAlternateAlleles();
        for ( int alleleIndex = 0; alleleIndex < alternateAlleles.size() ; ++alleleIndex ) {

            final Allele altAllele = alternateAlleles.get(alleleIndex);

            if ( alleleIndex < existingAlleleAnnotations.size() ) {
                funcotatorAnnotationStringBuilder.append( existingAlleleAnnotations.get(alleleIndex) );
                funcotatorAnnotationStringBuilder.append(FIELD_DELIMITER);
            }

            funcotatorAnnotationStringBuilder.append(
                    funcotations.stream()
                            .filter(f -> f.getAltAllele().equals(altAllele) )
                            .map(f -> f.serializeToVcfString(manualAnnotationSerializedString))
                            .collect(Collectors.joining(FIELD_DELIMITER))
            );
            funcotatorAnnotationStringBuilder.append(",");
        }

        // We have a trailing "," - we need to remove it:
        funcotatorAnnotationStringBuilder.deleteCharAt(funcotatorAnnotationStringBuilder.length()-1);

        // Add our new annotation:
        variantContextOutputBuilder.attribute(FUNCOTATOR_VCF_FIELD_NAME, funcotatorAnnotationStringBuilder.toString());

        // Add the genotypes from the variant:
        variantContextOutputBuilder.genotypes( variant.getGenotypes() );

        // Render and add our VCF line:
        vcfWriter.add( variantContextOutputBuilder.make() );
    }

    //==================================================================================================================

    /**
     * Create a header for a VCF file.
     * Uses {@link VcfOutputRenderer#dataSourceFactories} to get a list of fields to report producing (preserving their order).
     * Includes fields from the manual annotation maps.
     * @return The {@link VCFHeader} object with relevant information for {@link Funcotator}.
     */
    private VCFHeader createVCFHeader() {

        final Set<VCFHeaderLine> headerLines = new LinkedHashSet<>();

        // Add all lines of our existing VCF header:
        headerLines.addAll( existingHeader.getMetaDataInInputOrder() );

        final String dataSourceFields = getDataSourceFieldNamesForHeader(dataSourceFactories);
        final String manualAnnotationFields = String.join( HEADER_LISTED_FIELD_DELIMITER, manualAnnotations.keySet() );

        // Add in the lines about Funcotations:
        headerLines.addAll(defaultToolVcfHeaderLines);
        headerLines.add(new VCFHeaderLine("Funcotator Version", Funcotator.VERSION + " | " + getDataSourceInfoString()));
        headerLines.add(new VCFInfoHeaderLine(FUNCOTATOR_VCF_FIELD_NAME, VCFHeaderLineCount.A,
                VCFHeaderLineType.String, "Functional annotation from the Funcotator tool.  Funcotation fields are: " +
                manualAnnotationFields + HEADER_LISTED_FIELD_DELIMITER + dataSourceFields)
        );

        // Create a new header and preserve the genotype sample names:
        return new VCFHeader(headerLines, existingHeader.getGenotypeSamples());
    }

    /**
     * Creates a {@link String} containing the field names from our {@link VcfOutputRenderer#dataSourceFactories} suitable for putting in the VCF header.
     * @param dataSourceFactories A {@link List} of {@link DataSourceFuncotationFactory} objects from which to pull field names.
     * @return A {@link String} containing the field names from our {@link VcfOutputRenderer#dataSourceFactories} suitable for putting in the VCF header.
     */
    private static String getDataSourceFieldNamesForHeader(final List<DataSourceFuncotationFactory> dataSourceFactories) {
        return dataSourceFactories.stream()
                        .map(DataSourceFuncotationFactory::getSupportedFuncotationFields)
                        .flatMap(LinkedHashSet::stream)
                        .map(Object::toString).collect(Collectors.joining(HEADER_LISTED_FIELD_DELIMITER));
    }
}
