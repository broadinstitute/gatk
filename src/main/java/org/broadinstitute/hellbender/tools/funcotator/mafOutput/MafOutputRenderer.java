package org.broadinstitute.hellbender.tools.funcotator.mafOutput;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.funcotator.DataSourceFuncotationFactory;
import org.broadinstitute.hellbender.tools.funcotator.Funcotation;
import org.broadinstitute.hellbender.tools.funcotator.Funcotator;
import org.broadinstitute.hellbender.tools.funcotator.OutputRenderer;
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

    /**
     * The string representing a comment in a MAF file.
     */
    protected static final String COMMENT_STRING = "#";

    /**
     * Value to insert into columns with unspecified annotations.
     */
    protected static final String UNKNOWN_VALUE_STRING = "__UNKNOWN__";

    /**
     * Value to insert into unused annotation columns.
     */
    protected static final String UNUSED_STRING = "NA";

    /**
     * Default set of columns to include in this {@link MafOutputRenderer}.
     * Order of the columns is preserved by the {@link LinkedHashMap}, while still being able to access each field via
     * the associated key.
     */
    private static final LinkedHashMap<String, String> defaultMap = new LinkedHashMap<>();

    /**
     * Map for: outputFieldName -> dataSourceFieldName1, dataSourceFieldName2 ...
     *
     * This map informs how to fill out outputFields using the input funcotations.
     *
     * Strings in the value list are in order of priority for use in the outputField.
     * That is, if both dataSourceFieldName1 and dataSourceFieldName2 are present as funcotation fields, then only
     * dataSourceFieldName1 will be used as the outputField.
     */
    private static final Map<String, List<String>> outputFieldNameMap = new LinkedHashMap<>();

    /**
     * Delimiter for fields in the output MAF file.
     */
    private static final String FIELD_DELIMITER = "\t";

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

    //==================================================================================================================
    // Constructors:

    /**
     * Create a {@link MafOutputRenderer}.
     * @param outputFilePath {@link Path} to output file (must not be null).
     * @param dataSources {@link List} of {@link DataSourceFuncotationFactory} to back our annotations (must not be null).
     * @param unaccountedForDefaultAnnotations {@link LinkedHashMap} of default annotations that must be added.
     * @param unaccountedForOverrideAnnotations {@link LinkedHashMap} of override annotations that must be added.
     */
    public MafOutputRenderer(final Path outputFilePath,
                             final List<DataSourceFuncotationFactory> dataSources,
                             final LinkedHashMap<String, String> unaccountedForDefaultAnnotations,
                             final LinkedHashMap<String, String> unaccountedForOverrideAnnotations) {
        this(outputFilePath, dataSources, unaccountedForDefaultAnnotations, unaccountedForOverrideAnnotations, new VCFHeader(Collections.emptySet()), new LinkedHashSet<>());
    }

    /**
     * Create a {@link MafOutputRenderer}.
     * @param outputFilePath {@link Path} to output file (must not be null).
     * @param dataSources {@link List} of {@link DataSourceFuncotationFactory} to back our annotations (must not be null).
     * @param unaccountedForDefaultAnnotations {@link LinkedHashMap} of default annotations that must be added.
     * @param unaccountedForOverrideAnnotations {@link LinkedHashMap} of override annotations that must be added.
     * @param inputFileHeader {@link VCFHeader} of input VCF file to preserve.
     * @param toolHeaderLines Lines to add to the header with information about Funcotator.
     */
    public MafOutputRenderer(final Path outputFilePath,
                             final List<DataSourceFuncotationFactory> dataSources,
                             final LinkedHashMap<String, String> unaccountedForDefaultAnnotations,
                             final LinkedHashMap<String, String> unaccountedForOverrideAnnotations,
                             final VCFHeader inputFileHeader,
                             final Set<String> toolHeaderLines) {

        // Set our internal variables from the input:
        this.outputFilePath = outputFilePath;
        this.toolHeaderLines = new LinkedHashSet<>(toolHeaderLines);
        this.inputFileHeader = inputFileHeader;
        dataSourceFactories = dataSources;

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
                        defaultMap.put(field, UNKNOWN_VALUE_STRING);
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
        defaultMap.put("Score", UNUSED_STRING);
        defaultMap.put("BAM_File", UNUSED_STRING);

        // Cache the manual annotation string so we can pass it easily into any Funcotations:
        manualAnnotationSerializedString = (manualAnnotations.size() != 0 ? FIELD_DELIMITER + String.join( FIELD_DELIMITER, manualAnnotations.values() ) + FIELD_DELIMITER : "");
    }

    //==================================================================================================================
    // Override Methods:

    @Override
    public void open() {
        try {
            printWriter = new PrintWriter(Files.newOutputStream(outputFilePath));
        }
        catch (final IOException ex) {
            throw new UserException("Error opening output file path: " + outputFilePath.toUri().toString(), ex);
        }
    }

    @Override
    public void close() {
        printWriter.flush();
        printWriter.close();
    }

    @Override
    public void write(final VariantContext variant, final List<Funcotation> funcotations) {

        // Loop through each alt allele in our variant:
        for ( final Allele altAllele : variant.getAlternateAlleles() ) {

            // Create our output maps:
            final LinkedHashMap<String, Object> outputMap = new LinkedHashMap<>(defaultMap);
            final LinkedHashMap<String, Object> extraFieldOutputMap = new LinkedHashMap<>();

            // Get our funcotations for this allele and add them to the output maps:
            for ( final Funcotation funcotation : funcotations ) {
                if ( funcotation.getAltAllele().equals(altAllele) ) {
                    // Add all the fields from the other funcotations into the extra field output:
                    for ( final String field : funcotation.getFieldNames() ) {
                        setField(extraFieldOutputMap, field, funcotation.getField(field));
                    }
                }
            }

            // Now add in our annotation overrides so they can be aliased correctly with the outputFieldNameMap:
            extraFieldOutputMap.putAll(overrideAnnotations);

            // Go through all output fields and see if any of the names in the value list are in our extraFieldOutputMap.
            // For any that match, we remove them from our extraFieldOutputMap and add them to the outputMap with the
            // correct key.
            for ( final Map.Entry<String, List<String>> entry : outputFieldNameMap.entrySet() ) {
                for ( final String fieldName : entry.getValue() ) {
                    if ( extraFieldOutputMap.containsKey(fieldName) ) {
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
            if ( ! hasWrittenHeader ) {
                writeHeader(mafCompliantOutputMap);
            }

            // Write the output (with manual annotations at the end):
            writeLine(
                    mafCompliantOutputMap.values().stream().collect(Collectors.joining(FIELD_DELIMITER))
                    + manualAnnotationSerializedString
            );
        }
    }


    //==================================================================================================================
    // Static Methods:

    //==================================================================================================================
    // Instance Methods:

    /**
     * Ensures each value in the given {@code outputMap} is compliant with MAF-expected values.
     * @param outputMap The {@link Map} of output field -> values to be checked for MAF compliance.
     * @return A {@link LinkedHashMap} of output field strings -> values to be written to the MAF file.
     */
    private LinkedHashMap<String, String> replaceFuncotationValuesWithMafCompliantValues(final Map<String, Object> outputMap ) {

        final LinkedHashMap<String, String> finalOutMap = new LinkedHashMap<>(outputMap.size());

        // Massage individual Key/Value pairs:
        for ( final String key : outputMap.keySet() ) {
            finalOutMap.put(key, mafTransform(key, outputMap.get(key).toString()) );
        }

        // Massage the OtherTranscripts field:
        finalOutMap.put("Other_Transcripts", finalOutMap.get("Other_Transcripts").replaceAll(VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER, "|"));

        // Massage the start/end/alleles in the case of INDELs
        // (Because MAF has different conventions from VCF for start/end positions of INDELs)
        if ( finalOutMap.get("Variant_Type").equals("INS") || finalOutMap.get("Variant_Type").equals("DEL") ) {

            final int refAlleleLength = finalOutMap.get("Reference_Allele").length();
            final int altAlleleLength = finalOutMap.get("Tumor_Seq_Allele2").length();

            // TODO: port these changes to GencodeFuncotationFactory (issue: https://github.com/broadinstitute/gatk/issues/4378)
            // Check to see if it's an insertion:
            if ( refAlleleLength < altAlleleLength ) {
                // We must:
                //    Remove the first N bases from the ALT_allele where N = length(ref_allele)
                //    Replace the ref_allele with "-"
                //    Replace the Tumor_Seq_Allele1 with "-"
                //    Set the End_Position to be Start_Position + 1 (All Insertions should have length 1 to represent the bases between which the insertion occurs).
                finalOutMap.put("Tumor_Seq_Allele2", finalOutMap.get("Tumor_Seq_Allele2").substring(refAlleleLength));
                finalOutMap.put("Reference_Allele", "-");
                finalOutMap.put("Tumor_Seq_Allele1", "-");
                finalOutMap.put("End_Position", String.valueOf(Integer.valueOf(finalOutMap.get("Start_Position")) + 1));
            }
            // Check to see if it's a deletion:
            else if ( refAlleleLength > altAlleleLength ) {
                // We must:
                //    Remove the first N bases from the REF_allele where N = length(alt_allele)
                //    Remove the first N bases from the Tumor_Seq_Allele1
                //    Replace the alt_allele with "-"
                //    Increment the Start_Position by 1 (start position should be inclusive of the first base deleted)
                //    Increment the End_Position by M-1 where M = length(ref_allele) (end position should be inclusive of the last base deleted)
                finalOutMap.put("Reference_Allele", finalOutMap.get("Reference_Allele").substring(altAlleleLength));
                finalOutMap.put("Tumor_Seq_Allele1", finalOutMap.get("Tumor_Seq_Allele1").substring(altAlleleLength));
                finalOutMap.put("Tumor_Seq_Allele2", "-");
                finalOutMap.put("Start_Position", String.valueOf(Integer.valueOf(finalOutMap.get("Start_Position")) + 1));
                finalOutMap.put("End_Position", String.valueOf(Integer.valueOf(finalOutMap.get("End_Position")) + refAlleleLength - 1));
            }
        }

        return finalOutMap;
    }

    /**
     * Transforms a given {@code value} to the equivalent MAF-valid value based on the given {@code key}.
     * @param key The {@code key} off of which to base the transformation.  This key is the final (transformed) key for output (i.e. the column name in the MAF file).
     * @param value The {@code value} to transform into a MAF-valid value.
     * @return The MAF-valid equivalent of the given {@code value}.
     */
    private String mafTransform(final String key, final String value) {
        if ( key.equals("Variant_Classification") ) {
            switch(value) {
                case "IN_FRAME_DEL":             return "In_Frame_Del";
                case "IN_FRAME_INS":             return "In_Frame_Ins";
                case "FRAME_SHIFT_INS":          return "Frame_Shift_Ins";
                case "FRAME_SHIFT_DEL":          return "Frame_Shift_Del";
                case "MISSENSE":                 return "Missense_Mutation";
                case "NONSENSE":                 return "Nonsense_Mutation";
                case "SILENT":                   return "Silent";
                case "SPLICE_SITE":              return "Splice_Site";
                case "DE_NOVO_START_IN_FRAME":
                case "DE_NOVO_START_OUT_FRAME":
                case "START_CODON_SNP":
                case "START_CODON_INS":
                case "START_CODON_DEL":          return "Translation_Start_Site";
                case "NONSTOP":                  return "Nonstop_Mutation";
                case "FIVE_PRIME_UTR":           return "5'UTR";
                case "THREE_PRIME_UTR":          return "3'UTR";
                case "FIVE_PRIME_FLANK":         return "5'Flank";
                case "INTRON":                   return "Intron";
                case "LINCRNA":                  return "RNA";
            }
        }
        else if ( key.equals("Chromosome") ) {
            if ( value.equals("chrM") ) {
                return "MT";
            }
            else if (value.startsWith("chr")) {
                return value.substring(3);
            }
        }
        else if ( key.equals("Other_Transcripts") ) {
            return value.replaceAll("IN_FRAME_DEL", "In_Frame_Del")
                .replaceAll( "IN_FRAME_DEL", "In_Frame_Del" )
                .replaceAll( "IN_FRAME_INS", "In_Frame_Ins" )
                .replaceAll( "FRAME_SHIFT_INS", "Frame_Shift_Ins" )
                .replaceAll( "FRAME_SHIFT_DEL", "Frame_Shift_Del" )
                .replaceAll( "MISSENSE", "Missense_Mutation" )
                .replaceAll( "NONSENSE", "Nonsense_Mutation" )
                .replaceAll( "SILENT", "Silent" )
                .replaceAll( "SPLICE_SITE", "Splice_Site" )
                .replaceAll( "DE_NOVO_START_IN_FRAME", "Translation_Start_Site" )
                .replaceAll( "DE_NOVO_START_OUT_FRAME", "Translation_Start_Site" )
                .replaceAll( "START_CODON_SNP", "Translation_Start_Site" )
                .replaceAll( "START_CODON_INS", "Translation_Start_Site" )
                .replaceAll( "START_CODON_DEL", "Translation_Start_Site" )
                .replaceAll( "NONSTOP", "Nonstop_Mutation" )
                .replaceAll( "FIVE_PRIME_UTR", "5'UTR" )
                .replaceAll( "THREE_PRIME_UTR", "3'UTR" )
                .replaceAll( "FIVE_PRIME_FLANK", "5'Flank" )
                .replaceAll( "INTRON", "Intron" )
                .replaceAll( "LINCRNA", "RNA" );
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
            outputMap.put(key, UNKNOWN_VALUE_STRING);
        }
    }

    /**
     * Write the given line to the {@link #printWriter}.
     * @param line The {@link String} to write as a line to the {@link #printWriter}.
     */
    private void writeLine(final String line) {
        printWriter.write(line + System.lineSeparator());
    }

    /**
     * Write the header to the output file.
     * @param outputMap A populated output map from which to derive the header columns.
     */
    protected void writeHeader(final LinkedHashMap<String, String> outputMap) {
        // Write out version:
        writeLine(COMMENT_STRING + "version " + VERSION);
        writeLine(COMMENT_STRING + COMMENT_STRING);

        // Write previous header info:
        for ( final VCFHeaderLine line : inputFileHeader.getMetaDataInInputOrder() ) {
            printWriter.write(COMMENT_STRING);
            printWriter.write(COMMENT_STRING);
            printWriter.write(" ");
            writeLine( line.toString() );
        }

        // Write any default tool header lines:
        for ( final String line : toolHeaderLines ) {
            printWriter.write(COMMENT_STRING);
            printWriter.write(COMMENT_STRING);
            printWriter.write(" ");
            writeLine(line);
        }

        // Write tool name and the data sources with versions:
        printWriter.write(COMMENT_STRING);
        printWriter.write(COMMENT_STRING);
        printWriter.write(" ");
        printWriter.write(" Funcotator ");
        printWriter.write(Funcotator.VERSION);
        printWriter.write(" | Date ");
        printWriter.write(new SimpleDateFormat("yyyymmdd'T'hhmmss").format(new Date()));
        printWriter.write(getDataSourceInfoString());
        writeLine("");

        // Write the column headers for our output set and our manual annotations:
        printWriter.write( outputMap.keySet().stream().collect(Collectors.joining(FIELD_DELIMITER)) );
        writeLine(FIELD_DELIMITER + manualAnnotations.keySet().stream().collect(Collectors.joining(FIELD_DELIMITER)));

        // Make sure we keep track of the fact that we've now written the header:
        hasWrittenHeader = true;
    }

    /**
     * Initializes {@link MafOutputRenderer#defaultMap} with the default keys for the columns in a MAF file.
     */
    protected void initializeDefaultMapWithKeys() {
        
        // Baseline required fields:
        defaultMap.put("Hugo_Symbol",                                 UNKNOWN_VALUE_STRING );
        defaultMap.put("Entrez_Gene_Id",                              UNKNOWN_VALUE_STRING );
        defaultMap.put("Center",                                      UNKNOWN_VALUE_STRING );
        defaultMap.put("NCBI_Build",                                  UNKNOWN_VALUE_STRING );
        defaultMap.put("Chromosome",                                  UNKNOWN_VALUE_STRING );
        defaultMap.put("Start_Position",                              UNKNOWN_VALUE_STRING );
        defaultMap.put("End_Position",                                UNKNOWN_VALUE_STRING );
        defaultMap.put("Strand",                                      "+" );
        defaultMap.put("Variant_Classification",                      UNKNOWN_VALUE_STRING );
        defaultMap.put("Variant_Type",                                UNKNOWN_VALUE_STRING );
        defaultMap.put("Reference_Allele",                            UNKNOWN_VALUE_STRING );
        defaultMap.put("Tumor_Seq_Allele1",                           UNKNOWN_VALUE_STRING );
        defaultMap.put("Tumor_Seq_Allele2",                           UNKNOWN_VALUE_STRING );
        defaultMap.put("dbSNP_RS",                                    UNKNOWN_VALUE_STRING );
        defaultMap.put("dbSNP_Val_Status",                            UNKNOWN_VALUE_STRING );
        defaultMap.put("Tumor_Sample_Barcode",                        UNKNOWN_VALUE_STRING );
        defaultMap.put("Matched_Norm_Sample_Barcode",                 UNKNOWN_VALUE_STRING );
        defaultMap.put("Match_Norm_Seq_Allele1",                      UNKNOWN_VALUE_STRING );
        defaultMap.put("Match_Norm_Seq_Allele2",                      UNKNOWN_VALUE_STRING );
        defaultMap.put("Tumor_Validation_Allele1",                    UNKNOWN_VALUE_STRING );
        defaultMap.put("Tumor_Validation_Allele2",                    UNKNOWN_VALUE_STRING );
        defaultMap.put("Match_Norm_Validation_Allele1",               UNKNOWN_VALUE_STRING );
        defaultMap.put("Match_Norm_Validation_Allele2",               UNKNOWN_VALUE_STRING );
        defaultMap.put("Verification_Status",                         UNKNOWN_VALUE_STRING );
        defaultMap.put("Validation_Status",                           UNKNOWN_VALUE_STRING );
        defaultMap.put("Mutation_Status",                             UNKNOWN_VALUE_STRING );
        defaultMap.put("Sequencing_Phase",                            UNKNOWN_VALUE_STRING );
        defaultMap.put("Sequence_Source",                             UNKNOWN_VALUE_STRING );
        defaultMap.put("Validation_Method",                           UNKNOWN_VALUE_STRING );
        defaultMap.put("Score",                                       UNKNOWN_VALUE_STRING );
        defaultMap.put("BAM_File",                                    UNKNOWN_VALUE_STRING );
        defaultMap.put("Sequencer",                                   UNKNOWN_VALUE_STRING );
        defaultMap.put("Tumor_Sample_UUID",                           UNKNOWN_VALUE_STRING );
        defaultMap.put("Matched_Norm_Sample_UUID",                    UNKNOWN_VALUE_STRING );

        // Required "optional" fields:
        defaultMap.put("Genome_Change",                               UNKNOWN_VALUE_STRING);
        defaultMap.put("Annotation_Transcript",                       UNKNOWN_VALUE_STRING);
        defaultMap.put("Transcript_Strand",                           UNKNOWN_VALUE_STRING);
        defaultMap.put("Transcript_Exon",                             UNKNOWN_VALUE_STRING);
        defaultMap.put("Transcript_Position",                         UNKNOWN_VALUE_STRING);
        defaultMap.put("cDNA_Change",                                 UNKNOWN_VALUE_STRING);
        defaultMap.put("Codon_Change",                                UNKNOWN_VALUE_STRING);
        defaultMap.put("Protein_Change",                              UNKNOWN_VALUE_STRING);
        defaultMap.put("Other_Transcripts",                           UNKNOWN_VALUE_STRING);
        defaultMap.put("Refseq_mRNA_Id",                              UNKNOWN_VALUE_STRING);
        defaultMap.put("Refseq_prot_Id",                              UNKNOWN_VALUE_STRING);
        defaultMap.put("SwissProt_acc_Id",                            UNKNOWN_VALUE_STRING);
        defaultMap.put("SwissProt_entry_Id",                          UNKNOWN_VALUE_STRING);
        defaultMap.put("Description",                                 UNKNOWN_VALUE_STRING);
        defaultMap.put("UniProt_AApos",                               UNKNOWN_VALUE_STRING);
        defaultMap.put("UniProt_Region",                              UNKNOWN_VALUE_STRING);
        defaultMap.put("UniProt_Site",                                UNKNOWN_VALUE_STRING);
        defaultMap.put("UniProt_Natural_Variations",                  UNKNOWN_VALUE_STRING);
        defaultMap.put("UniProt_Experimental_Info",                   UNKNOWN_VALUE_STRING);
        defaultMap.put("GO_Biological_Process",                       UNKNOWN_VALUE_STRING);
        defaultMap.put("GO_Cellular_Component",                       UNKNOWN_VALUE_STRING);
        defaultMap.put("GO_Molecular_Function",                       UNKNOWN_VALUE_STRING);
        defaultMap.put("COSMIC_overlapping_mutations",                UNKNOWN_VALUE_STRING);
        defaultMap.put("COSMIC_fusion_genes",                         UNKNOWN_VALUE_STRING);
        defaultMap.put("COSMIC_tissue_types_affected",                UNKNOWN_VALUE_STRING);
        defaultMap.put("COSMIC_total_alterations_in_gene",            UNKNOWN_VALUE_STRING);
        defaultMap.put("Tumorscape_Amplification_Peaks",              UNKNOWN_VALUE_STRING);
        defaultMap.put("Tumorscape_Deletion_Peaks",                   UNKNOWN_VALUE_STRING);
        defaultMap.put("TCGAscape_Amplification_Peaks",               UNKNOWN_VALUE_STRING);
        defaultMap.put("TCGAscape_Deletion_Peaks",                    UNKNOWN_VALUE_STRING);
        defaultMap.put("DrugBank",                                    UNKNOWN_VALUE_STRING);
        defaultMap.put("ref_context",                                 UNKNOWN_VALUE_STRING);
        defaultMap.put("gc_content",                                  UNKNOWN_VALUE_STRING);
        defaultMap.put("CCLE_ONCOMAP_overlapping_mutations",          UNKNOWN_VALUE_STRING);
        defaultMap.put("CCLE_ONCOMAP_total_mutations_in_gene",        UNKNOWN_VALUE_STRING);
        defaultMap.put("CGC_Mutation_Type",                           UNKNOWN_VALUE_STRING);
        defaultMap.put("CGC_Translocation_Partner",                   UNKNOWN_VALUE_STRING);
        defaultMap.put("CGC_Tumor_Types_Somatic",                     UNKNOWN_VALUE_STRING);
        defaultMap.put("CGC_Tumor_Types_Germline",                    UNKNOWN_VALUE_STRING);
        defaultMap.put("CGC_Other_Diseases",                          UNKNOWN_VALUE_STRING);
        defaultMap.put("DNARepairGenes_Activity_linked_to_OMIM",      UNKNOWN_VALUE_STRING);
        defaultMap.put("FamilialCancerDatabase_Syndromes",            UNKNOWN_VALUE_STRING);
        defaultMap.put("MUTSIG_Published_Results",                    UNKNOWN_VALUE_STRING);
        defaultMap.put("OREGANNO_ID",                                 UNKNOWN_VALUE_STRING);
        defaultMap.put("OREGANNO_Values",                             UNKNOWN_VALUE_STRING);
    }

    /**
     * Initializes the output map
     */
    private void initializeOutputFieldNameMap() {

        outputFieldNameMap.put( "Hugo_Symbol",                             Arrays.asList("Hugo_Symbol", "Gencode_19_hugoSymbol", "gene", "Gene") );
        outputFieldNameMap.put( "Entrez_Gene_Id",                          Arrays.asList("Entrez_Gene_Id", "HGNC_Entrez_Gene_ID", "HGNC_Entrez Gene ID", "HGNC_Entrez_Gene_ID(supplied_by_NCBI)", "HGNC_Entrez Gene ID(supplied by NCBI)", "entrez_id", "gene_id") );
        outputFieldNameMap.put( "Center",                                  Arrays.asList("Center", "center") );
        outputFieldNameMap.put( "NCBI_Build",                              Arrays.asList("NCBI_Build", "Gencode_19_ncbiBuild", "ncbi_build") );
        outputFieldNameMap.put( "Chromosome",                              Arrays.asList("Chromosome", "Gencode_19_chromosome", "chr", "contig", "chromosome", "chrom", "Chrom") );
        outputFieldNameMap.put( "Start_Position",                          Arrays.asList("Start_position", "Gencode_19_start", "start", "Start", "start_pos", "pos") );
        outputFieldNameMap.put( "End_Position",                            Arrays.asList("End_position", "Gencode_19_end", "end", "End", "end_pos") );
        outputFieldNameMap.put( "Strand",                                  Collections.singletonList("Strand") );
        outputFieldNameMap.put( "Variant_Classification",                  Arrays.asList("Variant_Classification", "Gencode_19_variantClassification", "variant_classification") );
        outputFieldNameMap.put( "Variant_Type",                            Arrays.asList("Variant_Type", "Gencode_19_variantType", "variant_type") );
        outputFieldNameMap.put( "Reference_Allele",                        Arrays.asList("Reference_Allele", "Gencode_19_refAllele", "ref", "ref_allele", "reference_allele") );
        outputFieldNameMap.put( "Tumor_Seq_Allele1",                       Arrays.asList("Tumor_Seq_Allele1", "Gencode_19_tumorSeqAllele1", "ref", "ref_allele", "reference_allele") );
        outputFieldNameMap.put( "Tumor_Seq_Allele2",                       Arrays.asList("Tumor_Seq_Allele2", "Gencode_19_tumorSeqAllele2", "alt", "alt_allele", "alt2", "alt_allele2", "alternate_allele2", "observed_allele2", "alternate_allele", "observed_allele", "alt1", "alt_allele1", "alternate_allele1", "observed_allele1") );
        outputFieldNameMap.put( "dbSNP_RS",                                Arrays.asList("dbSNP_RS", "dbsnp_rs") );
        outputFieldNameMap.put( "dbSNP_Val_Status",                        Arrays.asList("dbSNP_Val_Status", "dbsnp_val_status") );
        outputFieldNameMap.put( "Tumor_Sample_Barcode",                    Arrays.asList("Tumor_Sample_Barcode", "tumor_barcode", "tumor_id", "case_barcode", "case_id", "tumor_name") );
        outputFieldNameMap.put( "Matched_Norm_Sample_Barcode",             Arrays.asList("Matched_Norm_Sample_Barcode", "normal_barcode", "normal_id", "control_barcode", "control_id", "normal_name", "sample_name") );
        outputFieldNameMap.put( "Match_Norm_Seq_Allele1",                  Arrays.asList("Match_Norm_Seq_Allele1", "Match_Norm_Seq_Allele1") );
        outputFieldNameMap.put( "Match_Norm_Seq_Allele2",                  Arrays.asList("Match_Norm_Seq_Allele2", "Match_Norm_Seq_Allele2") );
        outputFieldNameMap.put( "Tumor_Validation_Allele1",                Arrays.asList("Tumor_Validation_Allele1", "Tumor_Validation_Allele1") );
        outputFieldNameMap.put( "Tumor_Validation_Allele2",                Arrays.asList("Tumor_Validation_Allele2", "Tumor_Validation_Allele2") );
        outputFieldNameMap.put( "Match_Norm_Validation_Allele1",           Arrays.asList("Match_Norm_Validation_Allele1", "Match_Norm_Validation_Allele1") );
        outputFieldNameMap.put( "Match_Norm_Validation_Allele2",           Arrays.asList("Match_Norm_Validation_Allele2", "Match_Norm_Validation_Allele2") );
        outputFieldNameMap.put( "Verification_Status",                     Arrays.asList("Verification_Status", "Verification_Status") );
        outputFieldNameMap.put( "Validation_Status",                       Arrays.asList("Validation_Status", "validation_status") );
        outputFieldNameMap.put( "Mutation_Status",                         Arrays.asList("Mutation_Status", "status") );
        outputFieldNameMap.put( "Sequencing_Phase",                        Arrays.asList("Sequencing_Phase", "phase") );
        outputFieldNameMap.put( "Sequence_Source",                         Arrays.asList("Sequence_Source", "source") );
        outputFieldNameMap.put( "Validation_Method",                       Arrays.asList("Validation_Method", "Validation_Method") );
        outputFieldNameMap.put( "Score",                                   Arrays.asList("Score", "Score") );
        outputFieldNameMap.put( "BAM_file",                                Arrays.asList("BAM_file", "bam", "bam_file") );
        outputFieldNameMap.put( "Sequencer",                               Arrays.asList("Sequencer", "sequencer", "platform") );
        outputFieldNameMap.put( "Tumor_Sample_UUID",                       Arrays.asList("Tumor_Sample_UUID", "tumor_uuid", "case_uuid", "tumor_barcode", "tumor_id", "case_barcode", "case_id", "tumor_name", "Tumor_Sample_Barcode") );
        outputFieldNameMap.put( "Matched_Norm_Sample_UUID",                Arrays.asList("Matched_Norm_Sample_UUID", "normal_uuid", "control_uuid", "normal_barcode", "normal_id", "control_barcode", "control_id", "normal_name", "sample_name", "Matched_Norm_Sample_Barcode") );

        outputFieldNameMap.put( "Genome_Change",                           Arrays.asList("Genome_Change", "Gencode_19_genomeChange", "genome_change") );
        outputFieldNameMap.put( "Annotation_Transcript",                   Arrays.asList("Annotation_Transcript", "Gencode_19_annotationTranscript", "annotation_transcript", "transcript_id") );
        outputFieldNameMap.put( "Transcript_Strand",                       Arrays.asList("Transcript_Strand", "Gencode_19_transcriptStrand", "transcript_strand") );
        outputFieldNameMap.put( "Transcript_Exon",                         Arrays.asList("Transcript_Exon", "Gencode_19_transcriptExon", "transcript_exon") );
        outputFieldNameMap.put( "Transcript_Position",                     Arrays.asList("Transcript_Position", "Gencode_19_transcriptPos", "transcript_position") );
        outputFieldNameMap.put( "cDNA_Change",                             Arrays.asList("cDNA_Change", "Gencode_19_cDnaChange", "transcript_change") );
        outputFieldNameMap.put( "Codon_Change",                            Arrays.asList("Codon_Change", "Gencode_19_codonChange", "codon_change") );
        outputFieldNameMap.put( "Protein_Change",                          Arrays.asList("Protein_Change", "Gencode_19_proteinChange", "protein_change") );
        outputFieldNameMap.put( "Other_Transcripts",                       Arrays.asList("Other_Transcripts", "Gencode_19_otherTranscripts", "other_transcripts") );
        outputFieldNameMap.put( "Refseq_mRNA_Id",                          Arrays.asList("Refseq_mRNA_Id", "Gencode_XRefSeq_mRNA_id", "gencode_xref_refseq_mRNA_id", "ENSEMBL_RefSeq_mRNA_accession", "RefSeq_mRNA_Id", "HGNC_RefSeq IDs") );
        outputFieldNameMap.put( "Refseq_prot_Id",                          Arrays.asList("Refseq_prot_Id", "Gencode_XRefSeq_prot_acc", "gencode_xref_refseq_prot_acc", "ENSEMBL_RefSeq_protein_accession", "RefSeq_prot_Id") );

        outputFieldNameMap.put( "SwissProt_acc_Id",                        Arrays.asList("SwissProt_acc_Id", "Simple_Uniprot_uniprot_accession", "uniprot_accession", "UniProt_uniprot_accession") );
        outputFieldNameMap.put( "SwissProt_entry_Id",                      Arrays.asList("SwissProt_entry_Id", "Simple_Uniprot_uniprot_entry_name", "uniprot_entry_name", "UniProt_uniprot_entry_name") );
        outputFieldNameMap.put( "Description",                             Arrays.asList("Description", "RefSeq_Description", "HGNC_Approved_Name", "HGNC_Approved Name") );

        outputFieldNameMap.put( "UniProt_AApos",                           Arrays.asList("UniProt_AApos", "UniProt_AAxform_aapos", "uniprot_AA_pos") );
        outputFieldNameMap.put( "UniProt_Region",                          Arrays.asList("UniProt_Region", "UniProt_AA_region") );
        outputFieldNameMap.put( "UniProt_Site",                            Arrays.asList("UniProt_Site", "UniProt_AA_site") );
        outputFieldNameMap.put( "UniProt_Natural_Variations",              Arrays.asList("UniProt_Natural_Variations", "UniProt_AA_natural_variation") );
        outputFieldNameMap.put( "UniProt_Experimental_Info",               Arrays.asList("UniProt_Experimental_Info", "UniProt_AA_experimental_info") );

        outputFieldNameMap.put( "GO_Biological_Process",                   Arrays.asList("GO_Biological_Process", "Simple_Uniprot_GO_Biological_Process", "UniProt_GO_Biological_Process") );
        outputFieldNameMap.put( "GO_Cellular_Component",                   Arrays.asList("GO_Cellular_Component", "Simple_Uniprot_GO_Cellular_Component", "UniProt_GO_Cellular_Component") );
        outputFieldNameMap.put( "GO_Molecular_Function",                   Arrays.asList("GO_Molecular_Function", "Simple_Uniprot_GO_Molecular_Function", "UniProt_GO_Molecular_Function") );
        outputFieldNameMap.put( "COSMIC_overlapping_mutations",            Arrays.asList("Cosmic_overlapping_mutations", "COSMIC_overlapping_mutations", "COSMIC_overlapping_mutation_AAs") );
        outputFieldNameMap.put( "COSMIC_fusion_genes",                     Arrays.asList("COSMIC_fusion_genes", "CosmicFusion_fusion_genes", "COSMIC_FusionGenes_fusion_genes") );
        outputFieldNameMap.put( "COSMIC_tissue_types_affected",            Arrays.asList("CosmicTissue_tissue_types_affected", "COSMIC_tissue_types_affected", "COSMIC_Tissue_tissue_types_affected") );
        outputFieldNameMap.put( "COSMIC_total_alterations_in_gene",        Arrays.asList("CosmicTissue_total_alterations_in_gene", "COSMIC_total_alterations_in_gene", "COSMIC_Tissue_total_alterations_in_gene") );
        outputFieldNameMap.put( "Tumorscape_Amplification_Peaks",          Arrays.asList("Tumorscape_Amplification_Peaks", "TUMORScape_Amplification_Peaks") );
        outputFieldNameMap.put( "Tumorscape_Deletion_Peaks",               Arrays.asList("Tumorscape_Deletion_Peaks", "TUMORScape_Deletion_Peaks") );
        outputFieldNameMap.put( "TCGAscape_Amplification_Peaks",           Arrays.asList("TCGAscape_Amplification_Peaks", "TCGAScape_Amplification_Peaks") );
        outputFieldNameMap.put( "TCGAscape_Deletion_Peaks",                Arrays.asList("TCGAscape_Deletion_Peaks", "TCGAScape_Deletion_Peaks") );
        outputFieldNameMap.put( "DrugBank",                                Arrays.asList("DrugBank", "Simple_Uniprot_DrugBank", "UniProt_DrugBank") );
        outputFieldNameMap.put( "ref_context",                             Arrays.asList("ref_context", "Gencode_19_referenceContext", "ref_context") );
        outputFieldNameMap.put( "gc_content",                              Arrays.asList("gc_content", "Gencode_19_gcContent", "gc_content") );
        outputFieldNameMap.put( "CCLE_ONCOMAP_overlapping_mutations",      Arrays.asList("CCLE_ONCOMAP_overlapping_mutations", "CCLE_By_GP_overlapping_mutations") );
        outputFieldNameMap.put( "CCLE_ONCOMAP_total_mutations_in_gene",    Arrays.asList("CCLE_ONCOMAP_total_mutations_in_gene", "CCLE_By_Gene_total_mutations_in_gene") );
        outputFieldNameMap.put( "CGC_Mutation_Type",                       Arrays.asList("CGC_Mutation_Type", "CGC_Mutation Type") );
        outputFieldNameMap.put( "CGC_Translocation_Partner",               Arrays.asList("CGC_Translocation_Partner", "CGC_Translocation Partner") );
        outputFieldNameMap.put( "CGC_Tumor_Types_Somatic",                 Arrays.asList("CGC_Tumor_Types_Somatic", "CGC_Tumour Types  (Somatic Mutations)") );
        outputFieldNameMap.put( "CGC_Tumor_Types_Germline",                Arrays.asList("CGC_Tumor_Types_Germline", "CGC_Tumour Types (Germline Mutations)") );
        outputFieldNameMap.put( "CGC_Other_Diseases",                      Arrays.asList("CGC_Other_Diseases", "CGC_Other Syndrome/Disease") );
        outputFieldNameMap.put( "DNARepairGenes_Activity_linked_to_OMIM",  Collections.singletonList("DNARepairGenes_Activity_linked_to_OMIM") );
        outputFieldNameMap.put( "FamilialCancerDatabase_Syndromes",        Arrays.asList("FamilialCancerDatabase_Syndromes", "Familial_Cancer_Genes_Syndrome") );
        outputFieldNameMap.put( "MUTSIG_Published_Results",                Arrays.asList("MUTSIG_Published_Results", "MutSig Published Results_Published_Results") );
        outputFieldNameMap.put( "OREGANNO_ID",                             Arrays.asList("OREGANNO_ID", "Oreganno_ID", "ORegAnno_ID") );
        outputFieldNameMap.put( "OREGANNO_Values",                         Arrays.asList("OREGANNO_Values", "Oreganno_Values", "ORegAnno_Values") );

        outputFieldNameMap.put( "tumor_f",                                 Arrays.asList("tumor_f", "sample_allelic_fraction", "AF") );
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

}
