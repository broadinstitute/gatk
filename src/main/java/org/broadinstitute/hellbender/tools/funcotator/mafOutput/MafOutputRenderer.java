package org.broadinstitute.hellbender.tools.funcotator.mafOutput;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.ImmutableMap;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.apache.commons.lang.StringUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.funcotator.*;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation;
import org.broadinstitute.hellbender.tools.funcotator.metadata.SamplePairExtractor;
import org.broadinstitute.hellbender.tools.funcotator.metadata.TumorNormalPair;
import org.broadinstitute.hellbender.tools.funcotator.vcfOutput.VcfOutputRenderer;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.stream.Collectors;

/**
 * A Funcotator output renderer for writing to MAF files.
 * Complies with version 2.4 of the MAF standard found here:
 * https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4
 * Created by jonn on 12/5/17.
 */
public class MafOutputRenderer extends OutputRenderer {

    //==================================================================================================================
    // Public Static Members:

    /**
     * Version of the MAF standard that this {@link MafOutputRenderer} writes.
     */
    public static String VERSION = "2.4";

    //==================================================================================================================
    // Private Static Members:

    static final Logger logger = LogManager.getLogger(MafOutputRenderer.class);

    /**
     * Default set of columns to include in this {@link MafOutputRenderer}.
     * Order of the columns is preserved by the {@link LinkedHashMap}, while still being able to access each field via
     * the associated key.
     */
    private final LinkedHashMap<String, String> defaultMap = new LinkedHashMap<>();

    /**
     * Map for: outputFieldName -> dataSourceFieldName1, dataSourceFieldName2 ...
     *
     * This map informs how to fill out outputFields using the input funcotations.
     *
     * Strings in the value list are in order of priority for use in the outputField.
     * That is, if both dataSourceFieldName1 and dataSourceFieldName2 are present as funcotation fields, then only
     * dataSourceFieldName1 will be used as the outputField.
     */
    private final Map<String, List<String>> outputFieldNameMap = new LinkedHashMap<>();

    //==================================================================================================================
    // Private Members:

    /** Flag to see if the header has been written to the output file yet. */
    private boolean hasWrittenHeader = false;

    /**
     * {@link Path} to the resulting MAF output file for this {@link MafOutputRenderer}.
     */
    private final Path outputFilePath;

    /**
     * {@link java.io.PrintWriter} to which to write the output MAF file.
     */
    private PrintWriter printWriter;

    /**
     * Tool header information to go into the header.
     */
    private final LinkedHashSet<String> toolHeaderLines;

    /**
     * Existing header from the input file to preserve in the output.
     */
    private final VCFHeader inputFileHeader;

    /** Override annotation list from the user. */
    private final LinkedHashMap<String, String> overrideAnnotations;

    /** The tumor normal pairs discovered in the input */
    private final List<TumorNormalPair> tnPairs;

    //==================================================================================================================
    // Constructors:

    /**
     * Create a {@link MafOutputRenderer}.
     * @param outputFilePath {@link Path} to output file (must not be null).
     * @param dataSources {@link List} of {@link DataSourceFuncotationFactory} to back our annotations (must not be null).
     * @param inputFileHeader {@link VCFHeader} of input VCF file to preserve.
     * @param unaccountedForDefaultAnnotations {@link LinkedHashMap} of default annotations that must be added.
     * @param unaccountedForOverrideAnnotations {@link LinkedHashMap} of override annotations that must be added.
     * @param toolHeaderLines Lines to add to the header with information about Funcotator.
     */
    public MafOutputRenderer(final Path outputFilePath,
                             final List<DataSourceFuncotationFactory> dataSources,
                             final VCFHeader inputFileHeader,
                             final LinkedHashMap<String, String> unaccountedForDefaultAnnotations,
                             final LinkedHashMap<String, String> unaccountedForOverrideAnnotations,
                             final Set<String> toolHeaderLines) {

        // Set our internal variables from the input:
        this.outputFilePath = outputFilePath;
        this.toolHeaderLines = new LinkedHashSet<>(toolHeaderLines);
        this.inputFileHeader = inputFileHeader;
        dataSourceFactories = dataSources;

        this.tnPairs = SamplePairExtractor.extractPossibleTumorNormalPairs(this.inputFileHeader);
        if (tnPairs.size() == 0) {
            logger.warn("No tumor/normal pairs were seen, cannot populate the some of the MAF fields (e.g. t_alt_count).  Please add '##tumor_sample=<tumor_sample_name>' and (if applicable) '##normal_sample=<normal_sample_name>' to the input VCF header");
        }

        // TODO: Make this check unnecessary and use the contained information to populate the MAF entries correctly.  (https://github.com/broadinstitute/gatk/issues/4912)
        if (this.tnPairs.size() > 1) {
            throw new UserException.BadInput("Input files with more than one tumor normal pair are currently not supported.  Found: " + tnPairs.stream()
                    .map(tn -> tn.toString())
                    .collect(Collectors.joining("; ")) );
        }

        // Merge the default annotations into our manualAnnotations:
        manualAnnotations = new LinkedHashMap<>();
        if ( unaccountedForDefaultAnnotations != null ) {
            manualAnnotations.putAll(unaccountedForDefaultAnnotations);
        }

        // Handle our override annotations a little differently:
        this.overrideAnnotations = unaccountedForOverrideAnnotations;

        // Fill in our default output map:
        initializeDefaultMapWithKeys();

        // Fill in our default map for outputField -> funcotation name:
        initializeOutputFieldNameMap();

        // TODO: Make this FASTER!
        // Update our defaultMap with the dataSourceFactories:
        for ( final DataSourceFuncotationFactory funcotationFactory : dataSourceFactories ) {
            for ( final String field : funcotationFactory.getSupportedFuncotationFields() ) {

                // Check if it's in our output map already:
                if ( !defaultMap.containsKey(field) ) {
                    boolean found = false;

                    // Now check our field assignment/translation map:
                    for ( final Map.Entry<String, List<String>> entry : outputFieldNameMap.entrySet() ) {
                        if ( entry.getValue().contains(field) ) {
                            found = true;
                            break;
                        }
                    }
                    if ( !found ) {
                        defaultMap.put(field, MafOutputRendererConstants.UNKNOWN_VALUE_STRING);
                    }
                }
            }
        }

        // Update the Default map with manual annotations that we have
        // and remove those annotations from the manual map:
        // NOTE: We must check the outputFieldNameMap in case we override an aliased output field:
        final Iterator<Map.Entry<String, String>> it = manualAnnotations.entrySet().iterator();
        while ( it.hasNext() ) {
            final Map.Entry<String, String> manualAnnotation = it.next();
            if ( defaultMap.containsKey(manualAnnotation.getKey()) ) {
                defaultMap.put(manualAnnotation.getKey(), manualAnnotation.getValue());
                it.remove();
            }
            else {
                // Check the aliased output fields.
                // NOTE: A field may be aliased many times, so we should check ALL fields before
                //       removing the iterator and moving to the next value:
                boolean isAliasedName = false;
                for ( final Map.Entry<String, List<String>> entry : outputFieldNameMap.entrySet() ) {
                    if ( entry.getValue().contains(manualAnnotation.getKey())) {
                        isAliasedName = true;
                        defaultMap.put( entry.getKey(), manualAnnotation.getValue() );
                    }
                }
                if (isAliasedName) {
                    it.remove();
                }
            }
        }

        // Set values for unused fields:
        defaultMap.put(MafOutputRendererConstants.FieldName_Score, MafOutputRendererConstants.UNUSED_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_BAM_File, MafOutputRendererConstants.UNUSED_STRING);

        // Cache the manual annotation string so we can pass it easily into any Funcotations:
        manualAnnotationSerializedString = (manualAnnotations.size() != 0 ? MafOutputRendererConstants.FIELD_DELIMITER + String.join( MafOutputRendererConstants.FIELD_DELIMITER, manualAnnotations.values() ) + MafOutputRendererConstants.FIELD_DELIMITER : "");

        // Open the output object:
        try {
            printWriter = new PrintWriter(Files.newOutputStream(outputFilePath));
        }
        catch (final IOException ex) {
            throw new UserException("Error opening output file path: " + outputFilePath.toUri().toString(), ex);
        }
    }

    //==================================================================================================================
    // Override Methods:

    @Override
    public void close() {
        if ( printWriter != null ) {
            printWriter.flush();
            printWriter.close();
        }
    }

    @Override
    public void write(final VariantContext variant, final FuncotationMap txToFuncotationMap) {

        if (txToFuncotationMap.getTranscriptList().size() > 1) {
            logger.warn("MAF typically does not support multiple transcripts per variant, though this should be able to render (grouped by transcript).  No user action needed.");
        }

        // Add the generated funcotations necessary for a MAF output rendering.
        final List<Funcotation> customMafFuncotations = CustomMafFuncotationCreator.createCustomMafCountFields(variant, tnPairs);
        txToFuncotationMap.getTranscriptList().forEach(txId -> txToFuncotationMap.add(txId, customMafFuncotations));

        // Loop through each alt allele in our variant:
        for ( final Allele altAllele : variant.getAlternateAlleles() ) {

            // Ignore spanning deletions.  Those have no meaning in a MAF.
            if (altAllele.equals(Allele.SPAN_DEL)) {
                continue;
            }

            for (final String txId : txToFuncotationMap.getTranscriptList()) {

                final List<Funcotation> funcotations = txToFuncotationMap.get(txId);

                // Create our output maps:
                final LinkedHashMap<String, Object> outputMap = new LinkedHashMap<>(defaultMap);
                final LinkedHashMap<String, Object> extraFieldOutputMap = new LinkedHashMap<>();

                // Get our funcotations for this allele and add them to the output maps:
                for (final Funcotation funcotation : funcotations) {
                    if (funcotation.getAltAllele().equals(altAllele)) {
                        // Add all the fields from the other funcotations into the extra field output:
                        for (final String field : funcotation.getFieldNames()) {
                            setField(extraFieldOutputMap, field, funcotation.getField(field));
                        }
                    }
                }

                // Now add in our annotation overrides so they can be aliased correctly with the outputFieldNameMap:
                extraFieldOutputMap.putAll(overrideAnnotations);

                // Go through all output fields and see if any of the names in the value list are in our extraFieldOutputMap.
                // For any that match, we remove them from our extraFieldOutputMap and add them to the outputMap with the
                // correct key.
                for (final Map.Entry<String, List<String>> entry : outputFieldNameMap.entrySet()) {
                    for (final String fieldName : entry.getValue()) {
                        if (extraFieldOutputMap.containsKey(fieldName)) {
                            outputMap.put(entry.getKey(), extraFieldOutputMap.remove(fieldName));
                            break;
                        }
                    }
                }

                // Merge our output maps together:
                outputMap.putAll(extraFieldOutputMap);

                // Now translate fields to the field names that MAF likes:
                final LinkedHashMap<String, String> mafCompliantOutputMap = replaceFuncotationValuesWithMafCompliantValues(outputMap);

                // Write our header if we have to:
                if (!hasWrittenHeader) {
                    writeHeader(mafCompliantOutputMap);
                }

                // Write the output (with manual annotations at the end):
                for (final Map.Entry<String, String> entry : mafCompliantOutputMap.entrySet()) {
                    writeString(entry.getValue());
                    writeString(MafOutputRendererConstants.FIELD_DELIMITER);
                }
                writeLine(manualAnnotationSerializedString);
            }
        }
    }

    //==================================================================================================================
    // Instance Methods:

    /**
     * @return The {@link #defaultMap} with currently populated fields.
     */
    @VisibleForTesting
    LinkedHashMap<String, String> getDefaultMap() {
        return defaultMap;
    }

    /**
     * Ensures each value in the given {@code outputMap} is compliant with MAF-expected values.
     * @param outputMap The {@link Map} of output field -> values to be checked for MAF compliance.
     * @return A {@link LinkedHashMap} of output field strings -> values to be written to the MAF file.
     */
    @VisibleForTesting
    LinkedHashMap<String, String> replaceFuncotationValuesWithMafCompliantValues( final Map<String, Object> outputMap ) {

        final LinkedHashMap<String, String> finalOutMap = new LinkedHashMap<>(outputMap.size());

        // Massage individual Key/Value pairs:
        for ( final String key : outputMap.keySet() ) {
            finalOutMap.put(key, mafTransform(key, outputMap.get(key).toString()) );
        }

        // Massage the OtherTranscripts field:
        if ( finalOutMap.containsKey(MafOutputRendererConstants.FieldName_Other_Transcripts) ) {
            finalOutMap.put(
                    MafOutputRendererConstants.FieldName_Other_Transcripts, finalOutMap.get(MafOutputRendererConstants.FieldName_Other_Transcripts)
                            .replaceAll(VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER, MafOutputRendererConstants.OTHER_TRANSCRIPT_DELIMITER)
            );
        }

        // Fix the alleles in the case of INDELS:
        adjustIndelAlleleInformationForMafOutput(finalOutMap);

        return finalOutMap;
    }

    /**
     * Checks the given {@code outputMap} has MAF-correct values for information relating to INDEL alleles.
     * NOTE: The output map is modified in place.
     * @param outputMap The {@link Map} of output field -> values to be checked for MAF compliance.
     */
    @VisibleForTesting
    void adjustIndelAlleleInformationForMafOutput(final LinkedHashMap<String, String> outputMap) {
        // Massage the start/end/alleles in the case of INDELs
        // (Because MAF has different conventions from VCF for start/end positions of INDELs)
        if ( outputMap.containsKey(MafOutputRendererConstants.FieldName_Variant_Type) &&
            (outputMap.get(MafOutputRendererConstants.FieldName_Variant_Type).equals(MafOutputRendererConstants.FieldValue_Variant_Type_Insertion) ||
             outputMap.get(MafOutputRendererConstants.FieldName_Variant_Type).equals(MafOutputRendererConstants.FieldValue_Variant_Type_Deletion)) ) {

            final int refAlleleLength = outputMap.get(MafOutputRendererConstants.FieldName_Reference_Allele).length();
            final int altAlleleLength = outputMap.get(MafOutputRendererConstants.FieldName_Tumor_Seq_Allele2).length();

            // TODO: port these changes to GencodeFuncotationFactory (issue: https://github.com/broadinstitute/gatk/issues/4378)
            // Check to see if it's an insertion:
            if ( refAlleleLength < altAlleleLength ) {
                // We must:
                //    Remove the first N bases from the ALT_allele where N = length(ref_allele)
                //    Replace the ref_allele with "-"
                //    Replace the Tumor_Seq_Allele1 with "-"
                //    Set the End_Position to be Start_Position + 1 (All Insertions should have length 1 to represent the bases between which the insertion occurs).
                outputMap.put(MafOutputRendererConstants.FieldName_Tumor_Seq_Allele2, outputMap.get(MafOutputRendererConstants.FieldName_Tumor_Seq_Allele2).substring(refAlleleLength));
                outputMap.put(MafOutputRendererConstants.FieldName_Reference_Allele,  MafOutputRendererConstants.EmptyAllele);
                outputMap.put(MafOutputRendererConstants.FieldName_Tumor_Seq_Allele1, MafOutputRendererConstants.EmptyAllele);
                outputMap.put(MafOutputRendererConstants.FieldName_End_Position, String.valueOf(Integer.valueOf(outputMap.get(MafOutputRendererConstants.FieldName_Start_Position)) + 1));
            }
            // Check to see if it's a deletion:
            else if ( refAlleleLength > altAlleleLength ) {
                // We must:
                //    Remove the first N bases from the REF_allele where N = length(alt_allele)
                //    Remove the first N bases from the Tumor_Seq_Allele1
                //    Replace the alt_allele with "-"
                //    Increment the Start_Position by 1 (start position should be inclusive of the first base deleted)
                //    Increment the End_Position by M-1 where M = length(ref_allele) (end position should be inclusive of the last base deleted)
                outputMap.put(MafOutputRendererConstants.FieldName_Reference_Allele,  outputMap.get(MafOutputRendererConstants.FieldName_Reference_Allele).substring(altAlleleLength));
                outputMap.put(MafOutputRendererConstants.FieldName_Tumor_Seq_Allele1, outputMap.get(MafOutputRendererConstants.FieldName_Tumor_Seq_Allele1).substring(altAlleleLength));
                outputMap.put(MafOutputRendererConstants.FieldName_Tumor_Seq_Allele2, MafOutputRendererConstants.EmptyAllele);
                outputMap.put(MafOutputRendererConstants.FieldName_Start_Position, String.valueOf(Integer.valueOf(outputMap.get(MafOutputRendererConstants.FieldName_Start_Position)) + 1));
                outputMap.put(MafOutputRendererConstants.FieldName_End_Position, String.valueOf(Integer.valueOf(outputMap.get(MafOutputRendererConstants.FieldName_End_Position)) + refAlleleLength - 1));
            }
        }
    }

    /**
     * Transforms a given {@code value} to the equivalent MAF-valid value based on the given {@code key}.
     * @param key The {@code key} off of which to base the transformation.  This key is the final (transformed) key for output (i.e. the column name in the MAF file).
     * @param value The {@code value} to transform into a MAF-valid value.
     * @return The MAF-valid equivalent of the given {@code value}.
     */
    @VisibleForTesting
    String mafTransform(final String key, final String value) {

        switch (key) {
            case MafOutputRendererConstants.FieldName_Variant_Classification:
                if ( MafOutputRendererConstants.VariantClassificationMap.containsKey(value)) {
                    return MafOutputRendererConstants.VariantClassificationMap.get(value);
                }
                break;
                
            case MafOutputRendererConstants.FieldName_Chromosome:
                if ( value.equals(MafOutputRendererConstants.FieldValue_Gencode_Chromosome_Mito) ) {
                    return MafOutputRendererConstants.FieldValue_Chromosome_Mito;
                }
                else if ( value.toLowerCase().startsWith("chr") ) {
                    return value.substring(3);
                }
                break;
            case MafOutputRendererConstants.FieldName_Other_Transcripts:

                // Use apache commons string utils because it's much, much faster to do this replacement:
                return StringUtils.replaceEachRepeatedly(value,
                        new String[]{GencodeFuncotation.VariantClassification.IN_FRAME_DEL.toString(),
                                GencodeFuncotation.VariantClassification.IN_FRAME_DEL.toString(),
                                GencodeFuncotation.VariantClassification.IN_FRAME_INS.toString(),
                                GencodeFuncotation.VariantClassification.FRAME_SHIFT_INS.toString(),
                                GencodeFuncotation.VariantClassification.FRAME_SHIFT_DEL.toString(),
                                GencodeFuncotation.VariantClassification.MISSENSE.toString(),
                                GencodeFuncotation.VariantClassification.NONSENSE.toString(),
                                GencodeFuncotation.VariantClassification.SILENT.toString(),
                                GencodeFuncotation.VariantClassification.SPLICE_SITE.toString(),
                                GencodeFuncotation.VariantClassification.DE_NOVO_START_IN_FRAME.toString(),
                                GencodeFuncotation.VariantClassification.DE_NOVO_START_OUT_FRAME.toString(),
                                GencodeFuncotation.VariantClassification.START_CODON_SNP.toString(),
                                GencodeFuncotation.VariantClassification.START_CODON_INS.toString(),
                                GencodeFuncotation.VariantClassification.START_CODON_DEL.toString(),
                                GencodeFuncotation.VariantClassification.NONSTOP.toString(),
                                GencodeFuncotation.VariantClassification.FIVE_PRIME_UTR.toString(),
                                GencodeFuncotation.VariantClassification.THREE_PRIME_UTR.toString(),
                                GencodeFuncotation.VariantClassification.FIVE_PRIME_FLANK.toString(),
                                GencodeFuncotation.VariantClassification.INTRON.toString(),
                                GencodeFuncotation.VariantClassification.LINCRNA.toString(),
                        },
                        new String[]{
                                MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.IN_FRAME_DEL.toString()),
                                MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.IN_FRAME_DEL.toString()),
                                MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.IN_FRAME_INS.toString()),
                                MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.FRAME_SHIFT_INS.toString()),
                                MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.FRAME_SHIFT_DEL.toString()),
                                MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.MISSENSE.toString()),
                                MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.NONSENSE.toString()),
                                MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.SILENT.toString()),
                                MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.SPLICE_SITE.toString()),
                                MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.DE_NOVO_START_IN_FRAME.toString()),
                                MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.DE_NOVO_START_OUT_FRAME.toString()),
                                MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.START_CODON_SNP.toString()),
                                MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.START_CODON_INS.toString()),
                                MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.START_CODON_DEL.toString()),
                                MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.NONSTOP.toString()),
                                MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.FIVE_PRIME_UTR.toString()),
                                MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.THREE_PRIME_UTR.toString()),
                                MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.FIVE_PRIME_FLANK.toString()),
                                MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.INTRON.toString()),
                                MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.LINCRNA.toString())
                            }
                        );
        }

        return value;
     }

    /**
     * Set the field in the given {@code outputMap} specified by the given {@code key} to the given {@code value}.
     * Checks the given {@code key} to see if the given {@code value} must be transformed to be a valid MAF field.
     * Will only set this field if the given {@code value} is not {@code null}.
     * @param outputMap {@link Map} of {@link String} to {@link Object} to hold output annotations.  Must not be {@code null}.
     * @param key {@link String} key to add to the output map.  Must not be {@code null}.
     * @param value {@link Object} value for the output map.
     */
    private void setField(final Map<String, Object> outputMap, final String key, final Object value) {
        Utils.nonNull(outputMap);
        Utils.nonNull(key);

        if ( value != null ) {
            outputMap.put(key, value);
        }
        else {
            outputMap.put(key, MafOutputRendererConstants.UNKNOWN_VALUE_STRING);
        }
    }

    /**
     * Write the given line to the {@link #printWriter} and append a newline.
     * @param line The {@link String} to write as a line to the {@link #printWriter}.
     */
    private void writeLine(final String line) {
        writeString(line);
        writeString(System.lineSeparator());
    }

    /**
     * Write the given {@link String} to the {@link #printWriter}.
     * @param s The {@link String} to write to the {@link #printWriter}.
     */
    private void writeString(final String s) {
        printWriter.write(s);
    }

    /**
     * Write the header to the output file.
     * @param outputMap A populated output map from which to derive the header columns.
     */
    protected void writeHeader(final LinkedHashMap<String, String> outputMap) {
        // Write out version:
        writeLine(MafOutputRendererConstants.COMMENT_STRING + "version " + VERSION);
        writeLine(MafOutputRendererConstants.COMMENT_STRING + MafOutputRendererConstants.COMMENT_STRING);

        // Write previous header info:
        for ( final VCFHeaderLine line : inputFileHeader.getMetaDataInInputOrder() ) {
            printWriter.write(MafOutputRendererConstants.COMMENT_STRING);
            printWriter.write(MafOutputRendererConstants.COMMENT_STRING);
            printWriter.write(" ");
            writeLine( line.toString() );
        }

        // Write any default tool header lines:
        for ( final String line : toolHeaderLines ) {
            printWriter.write(MafOutputRendererConstants.COMMENT_STRING);
            printWriter.write(MafOutputRendererConstants.COMMENT_STRING);
            printWriter.write(" ");
            writeLine(line);
        }

        // Write tool name and the data sources with versions:
        printWriter.write(MafOutputRendererConstants.COMMENT_STRING);
        printWriter.write(MafOutputRendererConstants.COMMENT_STRING);
        printWriter.write(" ");
        printWriter.write(" Funcotator ");
        printWriter.write(Funcotator.VERSION);
        printWriter.write(" | Date ");
        printWriter.write(new SimpleDateFormat("yyyymmdd'T'hhmmss").format(new Date()));
        printWriter.write(" | ");
        printWriter.write(getDataSourceInfoString());
        writeLine("");

        // Write the column headers for our output set and our manual annotations:
        printWriter.write( outputMap.keySet().stream().collect(Collectors.joining(MafOutputRendererConstants.FIELD_DELIMITER)) );
        writeLine(MafOutputRendererConstants.FIELD_DELIMITER + manualAnnotations.keySet().stream().collect(Collectors.joining(MafOutputRendererConstants.FIELD_DELIMITER)));

        // Make sure we keep track of the fact that we've now written the header:
        hasWrittenHeader = true;
    }

    /** No field ordering is preserved in the output of this method.
     * @return Immutable copy of the reverse output field map used by this MafOutputRenderer.
     * In other words, return a mapping for the input field name to the MAF column this will produce. Never {@code null}
     */
    public ImmutableMap<String, Set<String>> getReverseOutputFieldNameMap() {
        return  ImmutableMap.copyOf(Utils.getReverseValueToListMap(outputFieldNameMap));
    }

    /**
     * Initializes {@link MafOutputRenderer#defaultMap} with the default keys for the columns in a MAF file.
     */
    protected void initializeDefaultMapWithKeys() {
        
        // Baseline required fields:
        defaultMap.put(MafOutputRendererConstants.FieldName_Hugo_Symbol                   ,               MafOutputRendererConstants.UNKNOWN_VALUE_STRING );
        defaultMap.put(MafOutputRendererConstants.FieldName_Entrez_Gene_Id                ,               MafOutputRendererConstants.UNKNOWN_VALUE_STRING );
        defaultMap.put(MafOutputRendererConstants.FieldName_Center                        ,               MafOutputRendererConstants.UNKNOWN_VALUE_STRING );
        defaultMap.put(MafOutputRendererConstants.FieldName_NCBI_Build                    ,               MafOutputRendererConstants.UNKNOWN_VALUE_STRING );
        defaultMap.put(MafOutputRendererConstants.FieldName_Chromosome                    ,               MafOutputRendererConstants.UNKNOWN_VALUE_STRING );
        defaultMap.put(MafOutputRendererConstants.FieldName_Start_Position                ,               MafOutputRendererConstants.UNKNOWN_VALUE_STRING );
        defaultMap.put(MafOutputRendererConstants.FieldName_End_Position                  ,               MafOutputRendererConstants.UNKNOWN_VALUE_STRING );
        defaultMap.put(MafOutputRendererConstants.FieldName_Strand                        ,               MafOutputRendererConstants.FieldValue_Strand );
        defaultMap.put(MafOutputRendererConstants.FieldName_Variant_Classification        ,               MafOutputRendererConstants.UNKNOWN_VALUE_STRING );
        defaultMap.put(MafOutputRendererConstants.FieldName_Variant_Type                  ,               MafOutputRendererConstants.UNKNOWN_VALUE_STRING );
        defaultMap.put(MafOutputRendererConstants.FieldName_Reference_Allele              ,               MafOutputRendererConstants.UNKNOWN_VALUE_STRING );
        defaultMap.put(MafOutputRendererConstants.FieldName_Tumor_Seq_Allele1             ,               MafOutputRendererConstants.UNKNOWN_VALUE_STRING );
        defaultMap.put(MafOutputRendererConstants.FieldName_Tumor_Seq_Allele2             ,               MafOutputRendererConstants.UNKNOWN_VALUE_STRING );
        defaultMap.put(MafOutputRendererConstants.FieldName_dbSNP_RS                      ,               MafOutputRendererConstants.UNKNOWN_VALUE_STRING );
        defaultMap.put(MafOutputRendererConstants.FieldName_dbSNP_Val_Status              ,               MafOutputRendererConstants.UNKNOWN_VALUE_STRING );
        defaultMap.put(MafOutputRendererConstants.FieldName_Tumor_Sample_Barcode          ,               MafOutputRendererConstants.UNKNOWN_VALUE_STRING );
        defaultMap.put(MafOutputRendererConstants.FieldName_Matched_Norm_Sample_Barcode   ,               MafOutputRendererConstants.UNKNOWN_VALUE_STRING );
        defaultMap.put(MafOutputRendererConstants.FieldName_Match_Norm_Seq_Allele1        ,               MafOutputRendererConstants.UNKNOWN_VALUE_STRING );
        defaultMap.put(MafOutputRendererConstants.FieldName_Match_Norm_Seq_Allele2        ,               MafOutputRendererConstants.UNKNOWN_VALUE_STRING );
        defaultMap.put(MafOutputRendererConstants.FieldName_Tumor_Validation_Allele1      ,               MafOutputRendererConstants.UNKNOWN_VALUE_STRING );
        defaultMap.put(MafOutputRendererConstants.FieldName_Tumor_Validation_Allele2      ,               MafOutputRendererConstants.UNKNOWN_VALUE_STRING );
        defaultMap.put(MafOutputRendererConstants.FieldName_Match_Norm_Validation_Allele1 ,               MafOutputRendererConstants.UNKNOWN_VALUE_STRING );
        defaultMap.put(MafOutputRendererConstants.FieldName_Match_Norm_Validation_Allele2 ,               MafOutputRendererConstants.UNKNOWN_VALUE_STRING );
        defaultMap.put(MafOutputRendererConstants.FieldName_Verification_Status           ,               MafOutputRendererConstants.UNKNOWN_VALUE_STRING );
        defaultMap.put(MafOutputRendererConstants.FieldName_Validation_Status             ,               MafOutputRendererConstants.UNKNOWN_VALUE_STRING );
        defaultMap.put(MafOutputRendererConstants.FieldName_Mutation_Status               ,               MafOutputRendererConstants.UNKNOWN_VALUE_STRING );
        defaultMap.put(MafOutputRendererConstants.FieldName_Sequencing_Phase              ,               MafOutputRendererConstants.UNKNOWN_VALUE_STRING );
        defaultMap.put(MafOutputRendererConstants.FieldName_Sequence_Source               ,               MafOutputRendererConstants.UNKNOWN_VALUE_STRING );
        defaultMap.put(MafOutputRendererConstants.FieldName_Validation_Method             ,               MafOutputRendererConstants.UNKNOWN_VALUE_STRING );
        defaultMap.put(MafOutputRendererConstants.FieldName_Score                         ,               MafOutputRendererConstants.UNKNOWN_VALUE_STRING );
        defaultMap.put(MafOutputRendererConstants.FieldName_BAM_File                      ,               MafOutputRendererConstants.UNKNOWN_VALUE_STRING );
        defaultMap.put(MafOutputRendererConstants.FieldName_Sequencer                     ,               MafOutputRendererConstants.UNKNOWN_VALUE_STRING );
        defaultMap.put(MafOutputRendererConstants.FieldName_Tumor_Sample_UUID             ,               MafOutputRendererConstants.UNKNOWN_VALUE_STRING );
        defaultMap.put(MafOutputRendererConstants.FieldName_Matched_Norm_Sample_UUID      ,               MafOutputRendererConstants.UNKNOWN_VALUE_STRING );

        // Required "optional" fields:
        defaultMap.put(MafOutputRendererConstants.FieldName_Genome_Change                          ,      MafOutputRendererConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_Annotation_Transcript                  ,      MafOutputRendererConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_Transcript_Strand                      ,      MafOutputRendererConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_Transcript_Exon                        ,      MafOutputRendererConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_Transcript_Position                    ,      MafOutputRendererConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_cDNA_Change                            ,      MafOutputRendererConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_Codon_Change                           ,      MafOutputRendererConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_Protein_Change                         ,      MafOutputRendererConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_Other_Transcripts                      ,      MafOutputRendererConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_Refseq_mRNA_Id                         ,      MafOutputRendererConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_Refseq_prot_Id                         ,      MafOutputRendererConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_SwissProt_acc_Id                       ,      MafOutputRendererConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_SwissProt_entry_Id                     ,      MafOutputRendererConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_Description                            ,      MafOutputRendererConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_UniProt_AApos                          ,      MafOutputRendererConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_UniProt_Region                         ,      MafOutputRendererConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_UniProt_Site                           ,      MafOutputRendererConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_UniProt_Natural_Variations             ,      MafOutputRendererConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_UniProt_Experimental_Info              ,      MafOutputRendererConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_GO_Biological_Process                  ,      MafOutputRendererConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_GO_Cellular_Component                  ,      MafOutputRendererConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_GO_Molecular_Function                  ,      MafOutputRendererConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_COSMIC_overlapping_mutations           ,      MafOutputRendererConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_COSMIC_fusion_genes                    ,      MafOutputRendererConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_COSMIC_tissue_types_affected           ,      MafOutputRendererConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_COSMIC_total_alterations_in_gene       ,      MafOutputRendererConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_Tumorscape_Amplification_Peaks         ,      MafOutputRendererConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_Tumorscape_Deletion_Peaks              ,      MafOutputRendererConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_TCGAscape_Amplification_Peaks          ,      MafOutputRendererConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_TCGAscape_Deletion_Peaks               ,      MafOutputRendererConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_DrugBank                               ,      MafOutputRendererConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_ref_context                            ,      MafOutputRendererConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_gc_content                             ,      MafOutputRendererConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_CCLE_ONCOMAP_overlapping_mutations     ,      MafOutputRendererConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_CCLE_ONCOMAP_total_mutations_in_gene   ,      MafOutputRendererConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_CGC_Mutation_Type                      ,      MafOutputRendererConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_CGC_Translocation_Partner              ,      MafOutputRendererConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_CGC_Tumor_Types_Somatic                ,      MafOutputRendererConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_CGC_Tumor_Types_Germline               ,      MafOutputRendererConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_CGC_Other_Diseases                     ,      MafOutputRendererConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_DNARepairGenes_Activity_linked_to_OMIM ,      MafOutputRendererConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_FamilialCancerDatabase_Syndromes       ,      MafOutputRendererConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_MUTSIG_Published_Results               ,      MafOutputRendererConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_OREGANNO_ID                            ,      MafOutputRendererConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_OREGANNO_Values                        ,      MafOutputRendererConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_tumor_f                                ,      MafOutputRendererConstants.UNKNOWN_VALUE_STRING);
        defaultMap.put(MafOutputRendererConstants.FieldName_t_alt_count                            ,      MafOutputRendererConstants.UNKNOWN_VALUE_STRING );
        defaultMap.put(MafOutputRendererConstants.FieldName_t_ref_count                            ,      MafOutputRendererConstants.UNKNOWN_VALUE_STRING );
        defaultMap.put(MafOutputRendererConstants.FieldName_n_alt_count                            ,      MafOutputRendererConstants.UNKNOWN_VALUE_STRING );
        defaultMap.put(MafOutputRendererConstants.FieldName_n_ref_count                            ,      MafOutputRendererConstants.UNKNOWN_VALUE_STRING );

    }

    /**
     * Initializes the output map
     */
    private void initializeOutputFieldNameMap() {

        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Hugo_Symbol                            , MafOutputRendererConstants.OutputFieldNameMap_Hugo_Symbol );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Entrez_Gene_Id                         , MafOutputRendererConstants.OutputFieldNameMap_Entrez_Gene_Id );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Center                                 , MafOutputRendererConstants.OutputFieldNameMap_Center );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_NCBI_Build                             , MafOutputRendererConstants.OutputFieldNameMap_NCBI_Build );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Chromosome                             , MafOutputRendererConstants.OutputFieldNameMap_Chromosome );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Start_Position                         , MafOutputRendererConstants.OutputFieldNameMap_Start_Position );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_End_Position                           , MafOutputRendererConstants.OutputFieldNameMap_End_Position );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Strand                                 , MafOutputRendererConstants.OutputFieldNameMap_Strand );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Variant_Classification                 , MafOutputRendererConstants.OutputFieldNameMap_Variant_Classification );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Variant_Type                           , MafOutputRendererConstants.OutputFieldNameMap_Variant_Type );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Reference_Allele                       , MafOutputRendererConstants.OutputFieldNameMap_Reference_Allele );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Tumor_Seq_Allele1                      , MafOutputRendererConstants.OutputFieldNameMap_Tumor_Seq_Allele1 );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Tumor_Seq_Allele2                      , MafOutputRendererConstants.OutputFieldNameMap_Tumor_Seq_Allele2 );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_dbSNP_RS                               , MafOutputRendererConstants.OutputFieldNameMap_dbSNP_RS );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_dbSNP_Val_Status                       , MafOutputRendererConstants.OutputFieldNameMap_dbSNP_Val_Status );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Tumor_Sample_Barcode                   , MafOutputRendererConstants.OutputFieldNameMap_Tumor_Sample_Barcode );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Matched_Norm_Sample_Barcode            , MafOutputRendererConstants.OutputFieldNameMap_Matched_Norm_Sample_Barcode );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Match_Norm_Seq_Allele1                 , MafOutputRendererConstants.OutputFieldNameMap_Match_Norm_Seq_Allele1 );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Match_Norm_Seq_Allele2                 , MafOutputRendererConstants.OutputFieldNameMap_Match_Norm_Seq_Allele2 );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Tumor_Validation_Allele1               , MafOutputRendererConstants.OutputFieldNameMap_Tumor_Validation_Allele1 );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Tumor_Validation_Allele2               , MafOutputRendererConstants.OutputFieldNameMap_Tumor_Validation_Allele2 );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Match_Norm_Validation_Allele1          , MafOutputRendererConstants.OutputFieldNameMap_Match_Norm_Validation_Allele1 );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Match_Norm_Validation_Allele2          , MafOutputRendererConstants.OutputFieldNameMap_Match_Norm_Validation_Allele2 );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Verification_Status                    , MafOutputRendererConstants.OutputFieldNameMap_Verification_Status );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Validation_Status                      , MafOutputRendererConstants.OutputFieldNameMap_Validation_Status );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Mutation_Status                        , MafOutputRendererConstants.OutputFieldNameMap_Mutation_Status );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Sequencing_Phase                       , MafOutputRendererConstants.OutputFieldNameMap_Sequencing_Phase );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Sequence_Source                        , MafOutputRendererConstants.OutputFieldNameMap_Sequence_Source );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Validation_Method                      , MafOutputRendererConstants.OutputFieldNameMap_Validation_Method );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Score                                  , MafOutputRendererConstants.OutputFieldNameMap_Score );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_BAM_File                               , MafOutputRendererConstants.OutputFieldNameMap_BAM_File );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Sequencer                              , MafOutputRendererConstants.OutputFieldNameMap_Sequencer );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Tumor_Sample_UUID                      , MafOutputRendererConstants.OutputFieldNameMap_Tumor_Sample_UUID );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Matched_Norm_Sample_UUID               , MafOutputRendererConstants.OutputFieldNameMap_Matched_Norm_Sample_UUID );

        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Genome_Change                          , MafOutputRendererConstants.OutputFieldNameMap_Genome_Change );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Annotation_Transcript                  , MafOutputRendererConstants.OutputFieldNameMap_Annotation_Transcript );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Transcript_Strand                      , MafOutputRendererConstants.OutputFieldNameMap_Transcript_Strand );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Transcript_Exon                        , MafOutputRendererConstants.OutputFieldNameMap_Transcript_Exon );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Transcript_Position                    , MafOutputRendererConstants.OutputFieldNameMap_Transcript_Position );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_cDNA_Change                            , MafOutputRendererConstants.OutputFieldNameMap_cDNA_Change );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Codon_Change                           , MafOutputRendererConstants.OutputFieldNameMap_Codon_Change );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Protein_Change                         , MafOutputRendererConstants.OutputFieldNameMap_Protein_Change );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Other_Transcripts                      , MafOutputRendererConstants.OutputFieldNameMap_Other_Transcripts );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Refseq_mRNA_Id                         , MafOutputRendererConstants.OutputFieldNameMap_Refseq_mRNA_Id );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Refseq_prot_Id                         , MafOutputRendererConstants.OutputFieldNameMap_Refseq_prot_Id );

        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_SwissProt_acc_Id                       , MafOutputRendererConstants.OutputFieldNameMap_SwissProt_acc_Id );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_SwissProt_entry_Id                     , MafOutputRendererConstants.OutputFieldNameMap_SwissProt_entry_Id );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Description                            , MafOutputRendererConstants.OutputFieldNameMap_Description );

        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_UniProt_AApos                          , MafOutputRendererConstants.OutputFieldNameMap_UniProt_AApos );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_UniProt_Region                         , MafOutputRendererConstants.OutputFieldNameMap_UniProt_Region );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_UniProt_Site                           , MafOutputRendererConstants.OutputFieldNameMap_UniProt_Site );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_UniProt_Natural_Variations             , MafOutputRendererConstants.OutputFieldNameMap_UniProt_Natural_Variations );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_UniProt_Experimental_Info              , MafOutputRendererConstants.OutputFieldNameMap_UniProt_Experimental_Info );

        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_GO_Biological_Process                  , MafOutputRendererConstants.OutputFieldNameMap_GO_Biological_Process );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_GO_Cellular_Component                  , MafOutputRendererConstants.OutputFieldNameMap_GO_Cellular_Component );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_GO_Molecular_Function                  , MafOutputRendererConstants.OutputFieldNameMap_GO_Molecular_Function );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_COSMIC_overlapping_mutations           , MafOutputRendererConstants.OutputFieldNameMap_COSMIC_overlapping_mutations );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_COSMIC_fusion_genes                    , MafOutputRendererConstants.OutputFieldNameMap_COSMIC_fusion_genes );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_COSMIC_tissue_types_affected           , MafOutputRendererConstants.OutputFieldNameMap_COSMIC_tissue_types_affected );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_COSMIC_total_alterations_in_gene       , MafOutputRendererConstants.OutputFieldNameMap_COSMIC_total_alterations_in_gene );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Tumorscape_Amplification_Peaks         , MafOutputRendererConstants.OutputFieldNameMap_Tumorscape_Amplification_Peaks );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_Tumorscape_Deletion_Peaks              , MafOutputRendererConstants.OutputFieldNameMap_Tumorscape_Deletion_Peaks );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_TCGAscape_Amplification_Peaks          , MafOutputRendererConstants.OutputFieldNameMap_TCGAscape_Amplification_Peaks );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_TCGAscape_Deletion_Peaks               , MafOutputRendererConstants.OutputFieldNameMap_TCGAscape_Deletion_Peaks );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_DrugBank                               , MafOutputRendererConstants.OutputFieldNameMap_DrugBank );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_ref_context                            , MafOutputRendererConstants.OutputFieldNameMap_ref_context );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_gc_content                             , MafOutputRendererConstants.OutputFieldNameMap_gc_content );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_CCLE_ONCOMAP_overlapping_mutations     , MafOutputRendererConstants.OutputFieldNameMap_CCLE_ONCOMAP_overlapping_mutations );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_CCLE_ONCOMAP_total_mutations_in_gene   , MafOutputRendererConstants.OutputFieldNameMap_CCLE_ONCOMAP_total_mutations_in_gene );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_CGC_Mutation_Type                      , MafOutputRendererConstants.OutputFieldNameMap_CGC_Mutation_Type );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_CGC_Translocation_Partner              , MafOutputRendererConstants.OutputFieldNameMap_CGC_Translocation_Partner );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_CGC_Tumor_Types_Somatic                , MafOutputRendererConstants.OutputFieldNameMap_CGC_Tumor_Types_Somatic );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_CGC_Tumor_Types_Germline               , MafOutputRendererConstants.OutputFieldNameMap_CGC_Tumor_Types_Germline );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_CGC_Other_Diseases                     , MafOutputRendererConstants.OutputFieldNameMap_CGC_Other_Diseases );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_DNARepairGenes_Activity_linked_to_OMIM , MafOutputRendererConstants.OutputFieldNameMap_DNARepairGenes_Activity_linked_to_OMIM );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_FamilialCancerDatabase_Syndromes       , MafOutputRendererConstants.OutputFieldNameMap_FamilialCancerDatabase_Syndromes );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_MUTSIG_Published_Results               , MafOutputRendererConstants.OutputFieldNameMap_MUTSIG_Published_Results );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_OREGANNO_ID                            , MafOutputRendererConstants.OutputFieldNameMap_OREGANNO_ID );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_OREGANNO_Values                        , MafOutputRendererConstants.OutputFieldNameMap_OREGANNO_Values );

        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_tumor_f                                , MafOutputRendererConstants.OutputFieldNameMap_tumor_f );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_t_alt_count                            , MafOutputRendererConstants.OutputFieldNameMap_t_alt_count );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_t_ref_count                            , MafOutputRendererConstants.OutputFieldNameMap_t_ref_count );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_n_alt_count                            , MafOutputRendererConstants.OutputFieldNameMap_n_alt_count );
        outputFieldNameMap.put( MafOutputRendererConstants.FieldName_n_ref_count                            , MafOutputRendererConstants.OutputFieldNameMap_n_ref_count );
    }

    //==================================================================================================================
    // Helper Data Types:

    //------------------------------------------------------------------------------------------------------------------
    // Columns:
    
// Required:
// 1	Hugo_Symbol
// 2	Entrez_Gene_Id
// 3	Center
// 4	NCBI_Build
// 5	Chromosome
// 6	Start_Position
// 7	End_Position
// 8	Strand
// 9	Variant_Classification
//10	Variant_Type
//11	Reference_Allele
//12	Tumor_Seq_Allele1
//13	Tumor_Seq_Allele2
//14	dbSNP_RS
//15	dbSNP_Val_Status
//16	Tumor_Sample_Barcode
//17	Matched_Norm_Sample_Barcode
//18	Match_Norm_Seq_Allele1
//19	Match_Norm_Seq_Allele2
//20	Tumor_Validation_Allele1
//21	Tumor_Validation_Allele2
//22	Match_Norm_Validation_Allele1
//23	Match_Norm_Validation_Allele2
//24	Verification_Status
//25	Validation_Status
//26	Mutation_Status
//27	Sequencing_Phase
//28	Sequence_Source
//29	Validation_Method
//30	Score
//31	BAM_File
//32	Sequencer
//33    Tumor_Sample_UUID
//34    Matched_Norm_Sample_UUID

// "Optional" Required:
//35    Genome_Change
//36    Annotation_Transcript
//37    Transcript_Strand
//38    Transcript_Exon
//39    Transcript_Position
//40    cDNA_Change
//41    Codon_Change
//42    Protein_Change
//43    Other_Transcripts
//44    Refseq_mRNA_Id
//45    Refseq_prot_Id
//46    SwissProt_acc_Id
//47    SwissProt_entry_Id
//48    Description
//49    UniProt_AApos
//50    UniProt_Region
//51    UniProt_Site
//52    UniProt_Natural_Variations
//53    UniProt_Experimental_Info
//54    GO_Biological_Process
//55    GO_Cellular_Component
//56    GO_Molecular_Function
//57    COSMIC_overlap
//58    ping_mutations
//59    COSMIC_fusion_genes
//60    COSMIC_tissue_types_affected
//61    COSMIC_total_alterations_in_gene
//62    Tumorscape_Amplification_Peaks
//63    Tumorscape_Deletion_Peaks
//64    TCGAscape_Amplification_Peaks
//65    TCGAscape_Deletion_Peaks
//66    DrugBank
//67    ref_context
//68    gc_content
//69    CCLE_ONCOMAP_overlapping_mutations
//70    CCLE_ONCOMAP_total_mutations_in_gene
//71    CGC_Mutation_Type
//72    CGC_Translocation_Partner
//73    CGC_Tumor_Types_Somatic
//74    CGC_Tumor_Types_Germline
//75    CGC_Other_Diseases
//76    DNARepairGenes_Role
//77    FamilialCancerDatabase_Syndromes
//78    MUTSIG_Published_Results
//79    OREGANNO_ID
//80    OREGANNO_Values
//81    tumor_f

}
