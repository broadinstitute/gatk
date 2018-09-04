package org.broadinstitute.hellbender.tools.funcotator.vcfOutput;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.apache.commons.lang.StringUtils;
import org.broadinstitute.hellbender.tools.funcotator.*;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.DataSourceUtils;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotationBuilder;
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
    public static final String OTHER_TRANSCRIPT_DELIMITER = "/";

    /**
     * The delimiter to use when separating the information regarding a transcript.
     */
    public static final String ALL_TRANSCRIPT_DELIMITER = "#";

    /**
     * The delimiter to use when starting the information regarding a transcript.
     */
    public static final String START_TRANSCRIPT_DELIMITER = "[";

    /**
     * The delimiter to use when ending the information regarding a transcript.
     */
    public static final String END_TRANSCRIPT_DELIMITER = "]";

    /**
     * The delimiter used between the preamble and the actual fields in the description of the
     *  {@link #FUNCOTATOR_VCF_FIELD_NAME} attribute i the header.
     */
    public static final String DESCRIPTION_PREAMBLE_DELIMITER = ": ";

    //==================================================================================================================

    private final VariantContextWriter vcfWriter;
    private final VCFHeader existingHeader;

    private final LinkedHashSet<VCFHeaderLine> defaultToolVcfHeaderLines;

    //==================================================================================================================

    public VcfOutputRenderer(final VariantContextWriter vcfWriter,
                             final VCFHeader existingHeader,
                             final List<DataSourceFuncotationFactory> dataSources) {
        this(vcfWriter, dataSources, existingHeader, new LinkedHashMap<>(), new LinkedHashMap<>());
    }

    public VcfOutputRenderer(final VariantContextWriter vcfWriter,
                             final List<DataSourceFuncotationFactory> dataSources,
                             final VCFHeader existingHeader,
                             final LinkedHashMap<String, String> unaccountedForDefaultAnnotations,
                             final LinkedHashMap<String, String> unaccountedForOverrideAnnotations) {
        this(vcfWriter, dataSources, existingHeader, unaccountedForDefaultAnnotations, unaccountedForOverrideAnnotations, new LinkedHashSet<>());
    }

    public VcfOutputRenderer(final VariantContextWriter vcfWriter,
                             final List<DataSourceFuncotationFactory> dataSources,
                             final VCFHeader existingHeader,
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

        // Open the output file and set up the header:
        final VCFHeader newHeader = createVCFHeader();
        vcfWriter.writeHeader(newHeader);
    }

    //==================================================================================================================

    @Override
    public void close() {
        vcfWriter.close();
    }

    @Override
    public void write(final VariantContext variant, final FuncotationMap txToFuncotationMap) {

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

            for (final String txId : txToFuncotationMap.getTranscriptList()) {
                funcotatorAnnotationStringBuilder.append(START_TRANSCRIPT_DELIMITER);
                final List<Funcotation> funcotations = txToFuncotationMap.get(txId);
                funcotatorAnnotationStringBuilder.append(
                        funcotations.stream()
                                .filter(f -> f.getAltAllele().equals(altAllele))
                                .filter(f -> f.getFieldNames().size() > 0)
                                .filter(f -> !f.getDataSourceName().equals(FuncotatorConstants.DATASOURCE_NAME_FOR_INPUT_VCFS))
                                .map(VcfOutputRenderer::adjustIndelAlleleInformation)
                                .map(f -> retrieveSanitizedFuncotation(f, manualAnnotationSerializedString))
                                .collect(Collectors.joining(FIELD_DELIMITER))
                );
                funcotatorAnnotationStringBuilder.append(END_TRANSCRIPT_DELIMITER + ALL_TRANSCRIPT_DELIMITER);
            }
            // We have a trailing "#" - we need to remove it:
            funcotatorAnnotationStringBuilder.deleteCharAt(funcotatorAnnotationStringBuilder.length()-1);
            funcotatorAnnotationStringBuilder.append(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR);
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
     * Adjusts the given {@link Funcotation} if it is a {@link GencodeFuncotation}.
     * @param funcotation The {@link Funcotation} to adjust.
     */
    private static Funcotation adjustIndelAlleleInformation(final Funcotation funcotation) {
        if ( funcotation instanceof GencodeFuncotation ) {
            return adjustIndelAlleleInformation((GencodeFuncotation)funcotation);
        }
        return funcotation;
    }

    /**
     * Adjusts the given {@link GencodeFuncotation}'s start, end, reference, and alternate alleles if the variant to
     * which the funcotation is associated was an insertion or deletion.
     * Makes adjustments in place.
     * @param gencodeFuncotation The {@link GencodeFuncotation} to adjust.
     */
    private static GencodeFuncotation adjustIndelAlleleInformation(final GencodeFuncotation gencodeFuncotation) {

        final GencodeFuncotation outFuncotation = new GencodeFuncotationBuilder(gencodeFuncotation).build();

        if ( (gencodeFuncotation.getVariantType().equals(GencodeFuncotation.VariantType.DEL)) ||
                (gencodeFuncotation.getVariantType().equals(GencodeFuncotation.VariantType.INS)) ) {

            final int refAlleleLength = gencodeFuncotation.getRefAllele().length();
            final int altAlleleLength = gencodeFuncotation.getTumorSeqAllele2().length();

            // Check to see if it's an insertion:
            if ( refAlleleLength < altAlleleLength ) {
                // We must:
                //    Remove the first N bases from the ALT_allele where N = length(ref_allele)
                //    Replace the ref_allele with "-"
                //    Replace the Tumor_Seq_Allele1 with "-"
                //    Set the End_Position to be Start_Position + 1 (All Insertions should have length 1 to represent the bases between which the insertion occurs).
                outFuncotation.setTumorSeqAllele2(gencodeFuncotation.getTumorSeqAllele2().substring(refAlleleLength));
                outFuncotation.setRefAllele("-");
                outFuncotation.setEnd(gencodeFuncotation.getStart() + 1);
            }
            // Check to see if it's a deletion:
            else if ( refAlleleLength > altAlleleLength ) {
                // We must:
                //    Remove the first N bases from the REF_allele where N = length(alt_allele)
                //    Remove the first N bases from the Tumor_Seq_Allele1
                //    Replace the alt_allele with "-"
                //    Increment the Start_Position by 1 (start position should be inclusive of the first base deleted)
                //    Increment the End_Position by M-1 where M = length(ref_allele) (end position should be inclusive of the last base deleted)
                outFuncotation.setRefAllele(gencodeFuncotation.getRefAllele().substring(altAlleleLength));
                outFuncotation.setTumorSeqAllele2("-");
//                outputMap.put(MafOutputRendererConstants.FieldName_Tumor_Seq_Allele2, MafOutputRendererConstants.EmptyAllele);
                outFuncotation.setStart(gencodeFuncotation.getStart() + 1);
                outFuncotation.setEnd(gencodeFuncotation.getEnd() + refAlleleLength - 1);
            }
        }

        return outFuncotation;
    }

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

        // Construct (only) the field list delimited by HEADER_LISTED_FIELD_DELIMITER
        final String delimitedFields = StringUtils.isEmpty(manualAnnotationFields) ? dataSourceFields :
                manualAnnotationFields + HEADER_LISTED_FIELD_DELIMITER + dataSourceFields;

        // Add in the lines about Funcotations:
        headerLines.addAll(defaultToolVcfHeaderLines);
        headerLines.add(new VCFHeaderLine("Funcotator Version", Funcotator.VERSION + " | " + getDataSourceInfoString()));
        headerLines.add(new VCFInfoHeaderLine(FUNCOTATOR_VCF_FIELD_NAME, VCFHeaderLineCount.A,
                VCFHeaderLineType.String, "Functional annotation from the Funcotator tool.  Funcotation fields are"
                + DESCRIPTION_PREAMBLE_DELIMITER +
                delimitedFields)
        );

        // Create a new header and preserve the genotype sample names:
        return new VCFHeader(headerLines, existingHeader.getGenotypeSamples());
    }

    /**
     * Creates a {@link String} containing the field names from our {@link VcfOutputRenderer#dataSourceFactories} suitable for putting in the VCF header.
     *
     * Gencode annotations are put first and then the rest.
     *
     * @param dataSourceFactories A {@link List} of {@link DataSourceFuncotationFactory} objects from which to pull field names.
     * @return A {@link String} containing the field names from our {@link VcfOutputRenderer#dataSourceFactories} suitable for putting in the VCF header.
     */
    private static String getDataSourceFieldNamesForHeader(final List<DataSourceFuncotationFactory> dataSourceFactories) {
        return dataSourceFactories.stream().sorted(DataSourceUtils::datasourceComparator)
                        .map(DataSourceFuncotationFactory::getSupportedFuncotationFields)
                        .flatMap(LinkedHashSet::stream)
                        .map(Object::toString).collect(Collectors.joining(HEADER_LISTED_FIELD_DELIMITER));
    }

    private static String retrieveSanitizedFuncotation(final Funcotation funcotation, final String manualAnnotationSerializedString) {
       return funcotation.serializeToVcfString(manualAnnotationSerializedString);

    }
}
