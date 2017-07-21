package org.broadinstitute.hellbender.utils.codecs.GENCODE;

import htsjdk.samtools.util.Locus;
import htsjdk.tribble.Feature;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.UserException;
import scala.tools.nsc.Global;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

/**
 * A GTF Feature represents one row of a GTF File.
 * The specification of a GTF file is defined here:
 * http://mblab.wustl.edu/GTF22.html
 *
 * Currently only supports version 26 or greater!
 *
 * Created by jonn on 7/21/17.
 */
public abstract class GencodeGtfFeature implements Feature, Comparable<GencodeGtfFeature> {

    protected final Logger logger = LogManager.getLogger(this.getClass());

    // Metadata fields:

    /**
     * The relative order of this Feature.
     * Normally it is a line number indicating the position in the original data file of this feature.
     */
    private long                    featureOrderNumber              = -1;

    // Required base GTF Fields:
    private String                  chromosomeName;
    private AnnotationSource        annotationSource;
    private FeatureType             featureType;
    private int                     genomicStartLocation;
    private int                     genomicEndLocation;
    private GenomicStrand           genomicStrand;
    private GenomicPhase            genomicPhase;

    // "Required" GENCODE GTF Fields:
    private String                  geneId                          = null;
    private String                  transcriptId                    = null;
    private GeneTranscriptStatus    geneStatus                      = null;
    private GeneTranscriptType      geneType                        = null;
    private String                  geneName                        = null;
    private GeneTranscriptType      transcriptType                  = null;
    private GeneTranscriptStatus    transcriptStatus                = null;
    private String                  transcriptName                  = null;
    private int                     exonNumber                      = -1;
    private String                  exonId                          = null;
    private LocusLevel              locusLevel                      = null;

    // Optional GENCODE GTF Fields:
    private ArrayList<OptionalField<?>> optionalFields              = new ArrayList<>();

    // Optional General GTF Fields:
    private String anonymousOptionalFields                          = null;

    // ================================================================================================

    /**
     * Populate this GencodeGtfFeature with the given data.
     */
    protected GencodeGtfFeature(String[] gtfFields) {
        chromosomeName          = gtfFields[0];
        annotationSource        = AnnotationSource.valueOf( gtfFields[1] );
        featureType             = GencodeGtfFeature.FeatureType.valueOf( gtfFields[2].toLowerCase() );
        genomicStartLocation    = Integer.valueOf( gtfFields[3] );
        genomicEndLocation      = Integer.valueOf( gtfFields[4] );
        genomicStrand           = GenomicStrand.getEnum( gtfFields[6] );
        genomicPhase            = GenomicPhase.getEnum( gtfFields[7] );

        // Get the extra fields from the last column:
        String[] extraFields    = gtfFields[8].split(";");

        StringBuilder anonymousOptionalFieldBuilder = new StringBuilder();

        // Now there are "optional" fields to go through (some actually required, some actually optional),
        // But we need to match up the field names to the fields themselves:
        for ( String extraField : extraFields ) {

            String[] fieldParts = extraField.trim().split(" ");

            String fieldName = fieldParts[0].trim();

            // The value of the field may be between two quotes.
            // We remove them here.
            String fieldValue = fieldParts[1].trim().replaceAll("\"", "");

            OptionalField<?> optionalField = null;

            switch (fieldName) {
                // Find the right field to set:
                case "gene_id":
                    geneId = fieldValue;
                    break;
                case "transcript_id":
                    transcriptId = fieldValue;
                    break;
                case "gene_type":
                    geneType = GeneTranscriptType.getEnum(fieldValue);
                    break;
                case "gene_status":
                    geneStatus = GeneTranscriptStatus.valueOf(fieldValue);
                    break;
                case "gene_name":
                    geneName = fieldValue;
                    break;
                case "transcript_type":
                    transcriptType = GeneTranscriptType.getEnum(fieldValue);
                    break;
                case "transcript_status":
                    transcriptStatus = GeneTranscriptStatus.valueOf(fieldValue);
                    break;
                case "transcript_name":
                    transcriptName = fieldValue;
                    break;
                case "exon_number":
                    exonNumber = Integer.valueOf(fieldValue);
                    break;
                case "exon_id":
                    exonId = fieldValue;
                    break;
                case "level":
                    locusLevel = LocusLevel.getEnum(fieldValue);
                    break;
                case "tag":
                    optionalField = new OptionalField<>(fieldName, FeatureTag.getEnum(fieldValue));
                    break;
                case "ccdsid":
                    optionalField = new OptionalField<>(fieldName, fieldValue);
                    break;
                case "havana_gene":
                    optionalField = new OptionalField<>(fieldName, fieldValue);
                    break;
                case "havana_transcript":
                    optionalField = new OptionalField<>(fieldName, fieldValue);
                    break;
                case "protein_id":
                    optionalField = new OptionalField<>(fieldName, fieldValue);
                    break;
                case "ont":
                    optionalField = new OptionalField<>(fieldName, fieldValue);
                    break;
                case "transcript_support_level":
                    optionalField = new OptionalField<>(fieldName, TranscriptSupportLevel.getEnum(fieldValue));
                    break;
                case "remap_status":
                    optionalField = new OptionalField<>(fieldName, RemapStatus.valueOf(fieldValue));
                    break;
                case "remap_original_id":
                    optionalField = new OptionalField<>(fieldName, fieldValue);
                    break;
                case "remap_original_location":
                    optionalField = new OptionalField<>(fieldName, Long.valueOf(fieldValue));
                    break;
                case "remap_num_mappings":
                    optionalField = new OptionalField<>(fieldName, Long.valueOf(fieldValue));
                    break;
                case "remap_target_status":
                    optionalField = new OptionalField<>(fieldName, RemapTargetStatus.getEnum(fieldValue));
                    break;
                case "remap_substituted_missing_target":
                    optionalField = new OptionalField<>(fieldName, fieldValue);
                    break;
                default:
                    anonymousOptionalFieldBuilder.append(extraField + ";");
            }

            // If the optional field was good, we add it:
            if ( optionalField != null ) {
                optionalFields.add(optionalField);
            }
        }

        // Save our anonymous optional fields:
        if ( anonymousOptionalFieldBuilder.length() != 0 ) {
            anonymousOptionalFields = anonymousOptionalFieldBuilder.toString();
        }
    }

    /**
     * Populate this GencodeGtfFeature with the given data.
     */
    protected GencodeGtfFeature(long featureOrderNumber,
                            String chromosomeName,
                            AnnotationSource annotationSource,
                            FeatureType featureType,
                            int genomicStartLocation,
                            int genomicEndLocation,
                            GenomicStrand genomicStrand,
                            GenomicPhase genomicPhase,
                            String geneId,
                            String transcriptId,
                            GeneTranscriptType geneType,
                            GeneTranscriptStatus geneStatus,
                            String geneName,
                            GeneTranscriptType transcriptType,
                            GeneTranscriptStatus transcriptStatus,
                            String transcriptName,
                            int exonNumber,
                            String exonId,
                            LocusLevel locusLevel,
                            ArrayList<OptionalField<?>> optionalFields,
                            String anonymousOptionalFields) {

        this.featureOrderNumber = featureOrderNumber;
        this.chromosomeName = chromosomeName;
        this.annotationSource = annotationSource;
        this.featureType = featureType;
        this.genomicStartLocation = genomicStartLocation;
        this.genomicEndLocation = genomicEndLocation;
        this.genomicStrand = genomicStrand;
        this.genomicPhase = genomicPhase;
        this.geneId = geneId;
        this.transcriptId = transcriptId;
        this.geneType = geneType;
        this.geneStatus = geneStatus;
        this.geneName = geneName;
        this.transcriptType = transcriptType;
        this.transcriptStatus = transcriptStatus;
        this.transcriptName = transcriptName;
        this.exonNumber = exonNumber;
        this.exonId = exonId;
        this.locusLevel = locusLevel;

        if ( optionalFields != null ) {
            this.optionalFields = optionalFields;
        }

        this.anonymousOptionalFields = anonymousOptionalFields;
    }

    // ================================================================================================

    public static GencodeGtfFeature create(long featureOrderNumber,
                                           String chromosomeName,
                                           AnnotationSource annotationSource,
                                           FeatureType featureType,
                                           int genomicStartLocation,
                                           int genomicEndLocation,
                                           GenomicStrand genomicStrand,
                                           GenomicPhase genomicPhase,
                                           String geneId,
                                           String transcriptId,
                                           GeneTranscriptType geneType,
                                           GeneTranscriptStatus geneStatus,
                                           String geneName,
                                           GeneTranscriptType transcriptType,
                                           GeneTranscriptStatus transcriptStatus,
                                           String transcriptName,
                                           int exonNumber,
                                           String exonId,
                                           LocusLevel locusLevel,
                                           ArrayList<OptionalField<?>> optionalFields,
                                           String anonymousOptionalFields) {

        GencodeGtfFeature feature;

        // Figure out which kind of feature to make:
        switch (featureType) {
            case gene:
                feature = GencodeGtfGeneFeature.create(featureOrderNumber, chromosomeName, annotationSource, featureType, genomicStartLocation, genomicEndLocation, genomicStrand, genomicPhase, geneId, transcriptId, geneType, geneStatus, geneName, transcriptType, transcriptStatus, transcriptName, exonNumber, exonId, locusLevel, optionalFields, anonymousOptionalFields);
                break;
            case transcript:
                feature = GencodeGtfTranscriptFeature.create(featureOrderNumber, chromosomeName, annotationSource, featureType, genomicStartLocation, genomicEndLocation, genomicStrand, genomicPhase, geneId, transcriptId, geneType, geneStatus, geneName, transcriptType, transcriptStatus, transcriptName, exonNumber, exonId, locusLevel, optionalFields, anonymousOptionalFields);
                break;
            case exon:
                feature = GencodeGtfExonFeature.create(featureOrderNumber, chromosomeName, annotationSource, featureType, genomicStartLocation, genomicEndLocation, genomicStrand, genomicPhase, geneId, transcriptId, geneType, geneStatus, geneName, transcriptType, transcriptStatus, transcriptName, exonNumber, exonId, locusLevel, optionalFields, anonymousOptionalFields);
                break;
            case cds:
                feature = GencodeGtfCDSFeature.create(featureOrderNumber, chromosomeName, annotationSource, featureType, genomicStartLocation, genomicEndLocation, genomicStrand, genomicPhase, geneId, transcriptId, geneType, geneStatus, geneName, transcriptType, transcriptStatus, transcriptName, exonNumber, exonId, locusLevel, optionalFields, anonymousOptionalFields);
                break;
            case utr:
                feature = GencodeGtfUTRFeature.create(featureOrderNumber, chromosomeName, annotationSource, featureType, genomicStartLocation, genomicEndLocation, genomicStrand, genomicPhase, geneId, transcriptId, geneType, geneStatus, geneName, transcriptType, transcriptStatus, transcriptName, exonNumber, exonId, locusLevel, optionalFields, anonymousOptionalFields);
                break;
            case start_codon:
                feature = GencodeGtfStartCodonFeature.create(featureOrderNumber, chromosomeName, annotationSource, featureType, genomicStartLocation, genomicEndLocation, genomicStrand, genomicPhase, geneId, transcriptId, geneType, geneStatus, geneName, transcriptType, transcriptStatus, transcriptName, exonNumber, exonId, locusLevel, optionalFields, anonymousOptionalFields);
                break;
            case stop_codon:
                feature = GencodeGtfStopCodonFeature.create(featureOrderNumber, chromosomeName, annotationSource, featureType, genomicStartLocation, genomicEndLocation, genomicStrand, genomicPhase, geneId, transcriptId, geneType, geneStatus, geneName, transcriptType, transcriptStatus, transcriptName, exonNumber, exonId, locusLevel, optionalFields, anonymousOptionalFields);
                break;
            case selenocysteine:
                feature = GencodeGtfSelenocysteineFeature.create(featureOrderNumber, chromosomeName, annotationSource, featureType, genomicStartLocation, genomicEndLocation, genomicStrand, genomicPhase, geneId, transcriptId, geneType, geneStatus, geneName, transcriptType, transcriptStatus, transcriptName, exonNumber, exonId, locusLevel, optionalFields, anonymousOptionalFields);
                break;
            default:
                throw new UserException.MalformedFile("Unknown type of GencodeGtfFeature: " + featureType);

        }

        return feature;

    }

    /**
     * Create a {@link GencodeGtfFeature} based on a line from a Gencode GTF File.
     * @param gtfLine A line from a Gencode GTF File to convert into a {@link GencodeGtfFeature} object.
     * @return The {@link GencodeGtfFeature} representing the information in {@code gtfLine}
     */
    public static GencodeGtfFeature create(String gtfLine) {
        return create(gtfLine.split("\t"));
    }

    /**
     * Create a {@link GencodeGtfFeature} based on a line from a Gencode GTF File.
     * @param gtfFields A line from a Gencode GTF File split on the {@code \t} character.
     * @return The {@link GencodeGtfFeature} representing the information in {@code gtfLine}
     */
    public static GencodeGtfFeature create(String[] gtfFields) {

        // Ensure that the input data are superficially well-formed:
        if ( gtfFields.length != GencodeGtfCodec.NUM_COLUMNS ) {
            throw new UserException.MalformedFile("Invalid number of fields in the given GENCODE line " +
                    " - Given: " + gtfFields.length + " Expected: " + GencodeGtfCodec.NUM_COLUMNS);
        }

        GencodeGtfFeature feature = null;

        String featureType = gtfFields[2].toLowerCase();

        // Figure out which kind of feature to make:
        if ( featureType.equals("gene") ) {
            feature = GencodeGtfGeneFeature.create(gtfFields);
        }
        else if ( featureType.equals("transcript") ) {
            feature = GencodeGtfTranscriptFeature.create(gtfFields);
        }
        else if ( featureType.equals("exon") ) {
            feature = GencodeGtfExonFeature.create(gtfFields);
        }
        else if ( featureType.equals("cds") ) {
            feature = GencodeGtfCDSFeature.create(gtfFields);
        }
        else if ( featureType.equals("utr") ) {
            feature = GencodeGtfUTRFeature.create(gtfFields);
        }
        else if ( featureType.equals("start_codon") ) {
            feature = GencodeGtfStartCodonFeature.create(gtfFields);
        }
        else if ( featureType.equals("stop_codon") ) {
            feature = GencodeGtfStopCodonFeature.create(gtfFields);
        }
        else if ( featureType.equals("selenocysteine") ) {
            feature = GencodeGtfSelenocysteineFeature.create(gtfFields);
        }
        else {
            throw new UserException.MalformedFile("Unknown type of GencodeGtfFeature: " + featureType);
        }

        return feature;
    }

    // ================================================================================================

    @Override
    public String getContig() {
        return chromosomeName;
    }

    @Override
    public int getStart() {
        return genomicStartLocation;
    }

    @Override
    public int getEnd() {
        return genomicEndLocation;
    }

    // ================================================================================================

    @Override
    public boolean equals(Object other) {

        boolean isEqual = other instanceof GencodeGtfFeature;

        if ( isEqual ) {

            GencodeGtfFeature otherFeature = (GencodeGtfFeature) other;

            // It goes like this:
            //      If the field is a primitive, just check it.
            //      If the field is an object:
            //          check if both this and the other field is null
            //          if they are not, make sure this one is not null and do a comparison against the other.
            // All compacted together in short-circuited comparisons.
            isEqual =
                    (featureOrderNumber == otherFeature.featureOrderNumber) &&
                    (((chromosomeName == null) && (otherFeature.chromosomeName == null)) || ((chromosomeName != null) && chromosomeName.equals(otherFeature.chromosomeName))) &&
                    (((annotationSource == null) && (otherFeature.annotationSource == null)) || ((annotationSource != null) && annotationSource.equals(otherFeature.annotationSource))) &&
                    (((featureType == null) && (otherFeature.featureType == null)) || ((featureType != null) && featureType.equals(otherFeature.featureType))) &&
                    (genomicStartLocation == otherFeature.genomicStartLocation) &&
                    (genomicEndLocation == otherFeature.genomicEndLocation) &&
                    (((genomicStrand == null) && (otherFeature.genomicStrand == null)) || ((genomicStrand != null) && genomicStrand.equals(otherFeature.genomicStrand))) &&
                    (((genomicPhase == null) && (otherFeature.genomicPhase == null)) || ((genomicPhase != null) && genomicPhase.equals(otherFeature.genomicPhase))) &&
                    (((geneId == null) && (otherFeature.geneId == null)) || ((geneId != null) && geneId.equals(otherFeature.geneId))) &&
                    (((transcriptId == null) && (otherFeature.transcriptId == null)) || ((transcriptId != null) && transcriptId.equals(otherFeature.transcriptId))) &&
                    (((geneType == null) && (otherFeature.geneType == null)) || ((geneType != null) && geneType.equals(otherFeature.geneType))) &&
                    (geneStatus == otherFeature.geneStatus) &&
                    (((geneName == null) && (otherFeature.geneName == null)) || ((geneName != null) && geneName.equals(otherFeature.geneName))) &&
                    (((transcriptType == null) && (otherFeature.transcriptType == null)) || ((transcriptType != null) && transcriptType.equals(otherFeature.transcriptType))) &&
                    (transcriptStatus == otherFeature.transcriptStatus) &&
                    (((transcriptName == null) && (otherFeature.transcriptName == null)) || ((transcriptName != null) && transcriptName.equals(otherFeature.transcriptName))) &&
                    (exonNumber == otherFeature.exonNumber) &&
                    (((exonId == null) && (otherFeature.exonId == null)) || ((exonId != null) && exonId.equals(otherFeature.exonId))) &&
                    (((locusLevel == null) && (otherFeature.locusLevel == null)) || ((locusLevel != null) && locusLevel.equals(otherFeature.locusLevel))) &&
                    (((anonymousOptionalFields == null) && (otherFeature.anonymousOptionalFields == null)) || ((anonymousOptionalFields != null) && anonymousOptionalFields.equals(otherFeature.anonymousOptionalFields)) &&
                    optionalFields.equals(otherFeature.optionalFields));
        }

        return isEqual;
    }

    /**
     * Get all the features from this {@link GencodeGtfFeature} itself.
     * This is useful to get any subfeatures included in this {@link GencodeGtfFeature}.
     * @return A {@link List} of the features represented in this {@link GencodeGtfFeature}.
     */
    protected List<GencodeGtfFeature> getAllFeatures() {
        ArrayList<GencodeGtfFeature> list = new ArrayList<>();
        list.add(this);
        return list;
    }

    /**
     * Serializes this {@link GencodeGtfFeature} to a string.
     * @return a {@link String} representing this {@link GencodeGtfFeature}
     */
    private String serializeToString() {

        StringBuilder stringBuilder = new StringBuilder();

        stringBuilder.append( chromosomeName );
        stringBuilder.append( '\t' );
        stringBuilder.append( annotationSource );
        stringBuilder.append( '\t' );
        stringBuilder.append( featureType );
        stringBuilder.append( '\t' );
        stringBuilder.append( genomicStartLocation );
        stringBuilder.append( '\t' );
        stringBuilder.append( genomicEndLocation );
        stringBuilder.append( "\t.\t" );
        stringBuilder.append( genomicStrand );
        stringBuilder.append( '\t' );
        stringBuilder.append( genomicPhase );
        stringBuilder.append( '\t' );

        if ( geneId != null ) {
            stringBuilder.append("gene_id \"");
            stringBuilder.append(geneId);
            stringBuilder.append( "\"; " );
        }
        if ( transcriptId != null) {
            stringBuilder.append("transcript_id \"");
            stringBuilder.append(transcriptId);
            stringBuilder.append( "\"; " );
        }
        if ( geneType != null ) {
            stringBuilder.append("gene_type \"");
            stringBuilder.append(geneType);
            stringBuilder.append( "\"; " );
        }
        if ( geneStatus != null ) {
            stringBuilder.append("gene_status \"");
            stringBuilder.append(geneStatus);
            stringBuilder.append( "\"; " );
        }
        if ( geneName != null ) {
            stringBuilder.append("gene_name \"");
            stringBuilder.append(geneName);
            stringBuilder.append( "\"; " );
        }
        if ( transcriptType != null ) {
            stringBuilder.append("transcript_type \"");
            stringBuilder.append(transcriptType);
            stringBuilder.append( "\"; " );
        }
        if ( transcriptStatus != null ) {
            stringBuilder.append("transcript_status \"");
            stringBuilder.append(transcriptStatus);
            stringBuilder.append( "\"; " );
        }
        if ( transcriptName != null ) {
            stringBuilder.append("transcript_name \"");
            stringBuilder.append(transcriptName);
            stringBuilder.append( "\"; " );
        }
        if ( exonNumber != -1 ) {
            stringBuilder.append("exon_number ");
            stringBuilder.append(exonNumber);
            stringBuilder.append( "; " );
        }
        if ( exonId != null) {
            stringBuilder.append("exon_id \"");
            stringBuilder.append(exonId);
            stringBuilder.append( "\"; ");
        }
        if (locusLevel != null) {
            stringBuilder.append("level ");
            stringBuilder.append(locusLevel);
            stringBuilder.append("; ");
        }

        // = = = = = = = = = = = = = = = = = = = = = = =

        // Output our optional fields:
        stringBuilder.append(
                optionalFields.stream().map(Object::toString).collect(Collectors.joining(" "))
        );

        if ( anonymousOptionalFields != null ) {
            stringBuilder.append(anonymousOptionalFields);
        }

        return stringBuilder.toString().trim();
    }

    @Override
    public String toString() {
        StringBuilder stringBuilder = new StringBuilder();

        List<GencodeGtfFeature> features = getAllFeatures();
        Collections.sort( features );

        for ( GencodeGtfFeature feature : features ) {
            stringBuilder.append( feature.serializeToString() );
            stringBuilder.append("\n");
        }

        return stringBuilder.toString().trim();
    }

    @Override
    public int hashCode() {
        return this.serializeToString().hashCode();
    }

    // ================================================================================================

    public long getFeatureOrderNumber() { return featureOrderNumber; }

    public String getChromosomeName() {
        return chromosomeName;
    }

    public AnnotationSource getAnnotationSource() {
        return annotationSource;
    }

    public FeatureType getFeatureType() {
        return featureType;
    }

    public int getGenomicStartLocation() {
        return genomicStartLocation;
    }

    public int getGenomicEndLocation() {
        return genomicEndLocation;
    }

    public GenomicStrand getGenomicStrand() {
        return genomicStrand;
    }

    public GenomicPhase getGenomicPhase() {
        return genomicPhase;
    }

    public String getGeneId() {
        return geneId;
    }

    public String getTranscriptId() {
        return transcriptId;
    }

    public GeneTranscriptType getGeneType() {
        return geneType;
    }

    public String getGeneName() {
        return geneName;
    }

    public GeneTranscriptType getTranscriptType() {
        return transcriptType;
    }

    public String getTranscriptName() {
        return transcriptName;
    }

    public int getExonNumber() {
        return exonNumber;
    }

    public String getExonId() {
        return exonId;
    }

    public LocusLevel getLocusLevel() {
        return locusLevel;
    }

    public ArrayList<OptionalField<?>> getOptionalFields() {
        return optionalFields;
    }

    public String getAnonymousOptionalFields() {
        return anonymousOptionalFields;
    }

    public OptionalField<?> getOptionalField(String key) {
        for (OptionalField<?> optionalField : optionalFields) {
            if ( optionalField.getName().equals(key) ) {
                return optionalField;
            }
        }
        return null;
    }

    /**
     * Comparable interface implementation for {@link GencodeGtfFeature}.
     *
     * Order is determined by {@link #featureOrderNumber}
     *
     * @param other {@link GencodeGtfFeature} to which to compare
     * @return -1 if this < other; 0 if this == other; 1 if this > other
     */
    @Override
    public int compareTo(GencodeGtfFeature other) {
        return (int)(featureOrderNumber - other.featureOrderNumber);
    }

    /**
     * Checks if {@code other} is contained within this {@link GencodeGtfFeature}.
     * Comparison is made using {@link #genomicStartLocation} and {@link #genomicEndLocation} (both ends inclusive).
     * @param other {@link GencodeGtfFeature} of which to check the bounds.
     * @return true if {@code other} is contained within the bounds of this {@link GencodeGtfFeature}, false otherwise.
     */
    public boolean contains(GencodeGtfFeature other) {
        return (other.getStart() >= getStart()) && (other.getEnd() <= getEnd());
    }

    public void setFeatureOrderNumber(long featureOrderNumber) {
        this.featureOrderNumber = featureOrderNumber;
    }

    // ================================================================================================

    static public class OptionalField<T> {

        private String name;
        private T value;

        public OptionalField(String name, T value) {
            this.name = name;
            this.value = value;
        }

        public String getName() {
            return name;
        }

        public void setName(String name) {
            this.name = name;
        }

        public T getValue() {
            return value;
        }

        public void setValue(T value) {
            this.value = value;
        }

        @Override
        public String toString() {

            StringBuilder sb = new StringBuilder();

            sb.append(name);
            sb.append(" ");

            // We need to do some formatting for the numbers / non-numbers in the field:
            String valueString = value.toString();
            if ( valueString.matches("\\d\\d*") ) {
                sb.append(valueString);
                sb.append(";");
            }
            else {
                sb.append("\"");
                sb.append(valueString);
                sb.append("\";");
            }

            return sb.toString();
        }

        @Override
        public int hashCode() {
            return this.toString().hashCode();
        }

        @Override
        public boolean equals(Object other) {

            if ( !(other instanceof OptionalField) ) {
                return false;
            }

            OptionalField<?> otherOptionalField = (OptionalField<?>) other;

            return (name.equals(otherOptionalField.name)) &&
                    (value.equals(otherOptionalField.value));
        }
    }

    // ================================================================================================

    public enum AnnotationSource {
        ENSEMBL,
        HAVANA
    }

    public enum FeatureType {
        gene("gene"),
        transcript("transcript"),
        selenocysteine("Selenocysteine"),
        exon("exon"),
        cds("CDS"),
        start_codon("start_codon"),
        stop_codon("stop_codon"),
        utr("UTR");

        private String serialized;

        FeatureType(String serializedValue) { serialized = serializedValue; }

        @Override
        public String toString() { return serialized; }

        public static FeatureType getEnum(String s) {
            for( FeatureType val : values() ) {
                if(val.serialized.equalsIgnoreCase(s)) {
                    return val;
                }
            }
            throw new IllegalArgumentException();
        }
    }

    public enum GenomicStrand {
        FORWARD("+"),
        BACKWARD("-");

        private String serialized;

        GenomicStrand(String serializedValue) {
            serialized = serializedValue;
        }

        @Override
        public String toString() {
            return serialized;
        }

        public static GenomicStrand getEnum(String s) {
            for( GenomicStrand val : values() ) {
                if(val.serialized.equalsIgnoreCase(s)) {
                    return val;
                }
            }
            throw new IllegalArgumentException();
        }
    }

    public enum GenomicPhase {
        ZERO("0"),
        ONE ("1"),
        TWO ("2"),
        DOT (".");

        private String serialized;

        GenomicPhase(String serializedValue) {
            serialized = serializedValue;
        }

        @Override
        public String toString() {
            return serialized;
        }

        public static GenomicPhase getEnum(String s) {
            for( GenomicPhase val : values() ) {
                if(val.serialized.equalsIgnoreCase(s)) {
                    return val;
                }
            }
            throw new IllegalArgumentException();
        }
    }

    public enum GeneTranscriptType {
        // Immunoglobulin (Ig) variable chain and T-cell receptor (TcR) genes imported or annotated according to the IMGT (http://www.imgt.org/)
        IG_C_gene,
        IG_D_gene,
        IG_J_gene,
        IG_LV_gene,
        IG_V_gene,
        TR_C_gene,
        TR_J_gene,
        TR_V_gene,
        TR_D_gene,

        // Inactivated immunoglobulin gene.
        IG_pseudogene,
        IG_C_pseudogene,
        IG_J_pseudogene,
        IG_V_pseudogene,
        TR_V_pseudogene,
        TR_J_pseudogene,

        // Non-coding RNA predicted using sequences from Rfam (http://rfam.xfam.org/) and miRBase (http://www.mirbase.org/)
        Mt_rRNA,
        Mt_tRNA,
        miRNA,
        misc_RNA,
        rRNA,
        scRNA,
        snRNA,
        snoRNA,
        ribozyme,
        sRNA,
        scaRNA,

        // Non-coding RNA predicted to be pseudogene by the Ensembl pipeline
        Mt_tRNA_pseudogene,
        tRNA_pseudogene,
        snoRNA_pseudogene,
        snRNA_pseudogene,
        scRNA_pseudogene,
        rRNA_pseudogene,
        misc_RNA_pseudogene,
        miRNA_pseudogene,

        // To be Experimentally Confirmed. This is used for non-spliced EST clusters that have polyA features. This category has been specifically created for the ENCODE project to highlight regions that could indicate the presence of protein coding genes that require experimental validation, either by 5' RACE or RT-PCR to extend the transcripts, or by confirming expression of the putatively-encoded peptide with specific antibodies.
        TEC,

        // If the coding sequence (following the appropriate reference) of a transcript finishes >50bp from a downstream splice site then it is tagged as NMD. If the variant does not cover the full reference coding sequence then it is annotated as NMD if NMD is unavoidable i.e. no matter what the exon structure of the missing portion is the transcript will be subject to NMD.
        nonsense_mediated_decay,

        // Transcript that has polyA features (including signal) without a prior stop codon in the CDS, i.e. a non-genomic polyA tail attached directly to the CDS without 3' UTR. These transcripts are subject to degradation.
        non_stop_decay,

        // Alternatively spliced transcript believed to contain intronic sequence relative to other, coding, variants.
        retained_intron,

        // Contains an open reading frame (ORF).
        protein_coding,

        // Doesn't contain an ORF.
        processed_transcript,

        // Transcript which is known from the literature to not be protein coding.
        non_coding,

        // Transcript believed to be protein coding, but with more than one possible open reading frame.
        ambiguous_orf,

        // Long non-coding transcript in introns of a coding gene that does not overlap any exons.
        sense_intronic,

        // Long non-coding transcript that contains a coding gene in its intron on the same strand.
        sense_overlapping,

        // Has transcripts that overlap the genomic span (i.e. exon or introns) of a protein-coding locus on the opposite strand.
        antisense,

        known_ncrna,

        // Have homology to proteins but generally suffer from a disrupted coding sequence and an active homologous gene can be found at another locus. Sometimes these entries have an intact coding sequence or an open but truncated ORF, in which case there is other evidence used (for example genomic polyA stretches at the 3' end) to classify them as a pseudogene. Can be further classified as one of the following.
        pseudogene,

        // Pseudogene that lack introns and is thought to arise from reverse transcription of mRNA followed by reinsertion of DNA into the genome.
        processed_pseudogene,

        // Pseudogene owing to a SNP/DIP but in other individuals/haplotypes/strains the gene is translated.
        polymorphic_pseudogene,

        // Pseudogene owing to a reverse transcribed and re-inserted sequence.
        retrotransposed,

        // Pseudogene where protein homology or genomic structure indicates a pseudogene, but the presence of locus-specific transcripts indicates expression.
        transcribed_processed_pseudogene,
        transcribed_unprocessed_pseudogene,
        transcribed_unitary_pseudogene,

        // Pseudogene that has mass spec data suggesting that it is also translated.
        translated_processed_pseudogene,
        translated_unprocessed_pseudogene,

        // A species specific unprocessed pseudogene without a parent gene, as it has an active orthologue in another species.
        unitary_pseudogene,

        // Pseudogene that can contain introns since produced by gene duplication.
        unprocessed_pseudogene,

        // Used to tag mistakes in the public databases (Ensembl/SwissProt/Trembl)
        artifact,

        // Long, intervening noncoding (linc) RNA that can be found in evolutionarily conserved, intergenic regions.
        lincRNA,

        // Unspliced lncRNA that is several kb in size.
        macro_lncRNA,

        // Transcript where ditag and/or published experimental data strongly supports the existence of short non-coding transcripts transcribed from the 3'UTR.
        three_prime_overlapping_ncRNA,

        // Otherwise viable coding region omitted from this alternatively spliced transcript because the splice variation affects a region coding for a protein domain.
        disrupted_domain,

        // Short non coding RNA gene that forms part of the vault ribonucleoprotein complex.
        vaultRNA,

        // A non-coding locus that originates from within the promoter region of a protein-coding gene, with transcription proceeding in the opposite direction on the other strand.
        bidirectional_promoter_lncRNA;

        public static GeneTranscriptType getEnum(String s) {

            if (s.startsWith("3")) {
                s = "three_" + s.substring(1);
            }

            // Looks like sometimes RNA is spelled `rna`.
            // Here's a fix for that:
            s = s.replace( "rna", "RNA" );

            return GeneTranscriptType.valueOf(s);
        }

        @Override
        public String toString() {
            String s = super.toString();

            if ( s.startsWith("three_") ) {
                s = "3" + s.substring(7);
            }
            return s;
        }

    }

    public enum GeneTranscriptStatus {
        KNOWN,
        NOVEL,
        PUTATIVE
    }

    public enum LocusLevel {
        // Verified locus
        ONE("1"),
        // Manually annotated locus
        TWO("2"),
        // Automatically annotated locus
        THREE("3");

        private String serialized;

        LocusLevel(String serializedValue) {
            serialized = serializedValue;
        }

        @Override
        public String toString() {
            return serialized;
        }

        public static LocusLevel getEnum(String s) {
            for( LocusLevel val : values() ) {
                if(val.serialized.equalsIgnoreCase(s)) {
                    return val;
                }
            }
            throw new IllegalArgumentException();
        }
    }

    public enum FeatureTag {
        // 3' end extended based on RNA-seq data.
        three_nested_supported_extension,

        // 3' end extended based on RNA-seq data.
        three_standard_supported_extension,

        // annotated based on RNA-seq data.
        fourfivefour_RNA_Seq_supported,

        // 5' end extended based on RNA-seq data.
        five_nested_supported_extension,

        // 5' end extended based on RNA-seq data.
        five_standard_supported_extension,

        // shares an identical CDS but has alternative 5' UTR with respect to a reference variant.
        alternative_3_UTR,

        // shares an identical CDS but has alternative 3' UTR with respect to a reference variant.
        alternative_5_UTR,

        // (This flag corresponds to the older flag "appris_principal") Where the transcript expected to code for the main
        appris_principal_1,

        // (This flag corresponds to the older flag "appris_candidate_ccds") Where the APPRIS core modules are unable to choose a
        appris_principal_2,

        // Where the APPRIS core modules are unable to choose a clear principal variant and there more than one of the variants
        appris_principal_3,

        // (This flag corresponds to the Ensembl 78 flag "appris_candidate_longest_ccds") Where the APPRIS core modules are unable
        appris_principal_4,

        // (This flag corresponds to the Ensembl 78 flag "appris_candidate_longest_seq") Where the APPRIS core modules are unable
        appris_principal_5,

        // Candidate transcript(s) models that are conserved in at least three tested non-primate species.
        appris_alternative_1,

        // Candidate transcript(s) models that appear to be conserved in fewer than three tested non-primate species.
        appris_alternative_2,

        // ranscript expected to code for the main functional isoform based on a range of protein features (APPRIS pipeline).
        appris_principal,

        // where there is no single 'appris_principal' variant the main functional isoform will be translated from one of the
        appris_candidate,

        // he "appris_candidate" transcript that has an unique CCDS.
        appris_candidate_ccds,

        // where there is no 'appris_principal' variant, the candidate with highest APPRIS score is selected as the primary
        appris_candidate_highest_score,

        // where there is no 'appris_principal' variant, the longest of the 'appris_candidate' variants is selected as the primary
        appris_candidate_longest,

        // he "appris_candidate" transcripts where there are several CCDS, in this case APPRIS labels the longest CCDS.
        appris_candidate_longest_ccds,

        // where there is no "appris_candidate_ccds" or "appris_candidate_longest_ccds" variant, the longest protein of the
        appris_candidate_longest_seq,

        // identifies a subset of representative transcripts for each gene; prioritises full-length protein coding transcripts
        basic,

        // ranscript contains two confidently annotated CDSs. Support may come from eg proteomic data, cross-species conservation
        bicistronic,

        // ranscript 5' end overlaps ENCODE or Fantom CAGE cluster.
        CAGE_supported_TSS,

        // member of the consensus CDS gene set, confirming coding regions between ENSEMBL, UCSC, NCBI and HAVANA.
        CCDS,

        // he coding region end could not be confirmed.
        cds_end_NF,

        // he coding region start could not be confirmed.
        cds_start_NF,

        // ranscript QC checked using dotplot to identify features eg splice junctions, end of homology.
        dotter_confirmed,

        // an upstream ATG is used where a downstream ATG seems more evolutionary conserved.
        downstream_ATG,

        // ranscript was tested and confirmed experimentally.
        exp_conf,

        // locus consists of non-overlapping transcript fragments either because of genome assembly issues (i.e., gaps or
        fragmented_locus,

        // ranscript model contains all possible in-frame exons supported by homology, experimental evidence or conservation, but
        inferred_exon_combination,

        // ranscript model is not supported by a single piece of transcript evidence. May be supported by multiple fragments of
        inferred_transcript_model,

        // ranscript supported by transcript evidence that, while ampping best-in-genome, shows regions of poor sequence quality.
        low_sequence_quality,

        // he mRNA end could not be confirmed.
        mRNA_end_NF,

        // he mRNA start could not be confirmed.
        mRNA_start_NF,

        // in-frame type of variation where, at the acceptor site, some variants splice after the first AG and others after the
        NAGNAG_splice_site,

        // he locus is a host for small non-coding RNAs.
        ncRNA_host,

        // annotated based on RNA-seq data.
        nested_454_RNA_Seq_supported,

        // he transcript looks like it is subject to NMD but publications, experiments or conservation support the translation of
        NMD_exception,

        // codon if the transcript were longer but cannot currently be annotated as NMD as does not fulfil all criteria - most
        NMD_likely_if_extended,

        // he CDS has a non-ATG start and its validity is supported by publication or conservation.
        non_ATG_start,

        // he transcript has a non-canonical splice site conserved in other species.
        non_canonical_conserved,

        // he transcript has a non-canonical splice site explained by a genomic sequencing error.
        non_canonical_genome_sequence_error,

        // he transcript has a non-canonical splice site explained by other reasons.
        non_canonical_other,

        // he transcript has a non-canonical splice site explained by a SNP.
        non_canonical_polymorphism,

        // he transcript has a non-canonical splice site that needs experimental confirmation.
        non_canonical_TEC,

        // he transcript has a non-canonical splice site explained by a U12 intron (i.e. AT-AC splice site).
        non_canonical_U12,

        // a splice variant for which supporting evidence has not been submitted to databases, i.e. the model is based on
        non_submitted_evidence,

        // a transcript is supported by evidence from same species paralogous loci.
        not_best_in_genome_evidence,

        // evidence from other species was used to build model.
        not_organism_supported,

        // protein-coding locus with no paralogues or orthologs.
        orphan,

        // exon(s) of the locus overlap exon(s) of a readthrough transcript or a transcript belonging to another locus.
        overlapping_locus,

        // a low confidence upstream ATG existing in other coding variant would lead to NMD in this trancript, that uses the high
        overlapping_uORF,

        // annotation in the pseudo-autosomal region, which is duplicated between chromosomes X and Y.
        PAR,

        // member of the pseudogene set predicted by YALE, UCSC and HAVANA.
        pseudo_consens,

        // a transcript that overlaps two or more independent loci but is considered to belong to a third, separate locus.
        readthrough_transcript,

        // locus overlaps a sequence error or an assembly error in the reference genome that affects its annotation (e.g., 1 or
        reference_genome_error,

        // internal intron of CDS portion of transcript is retained.
        retained_intron_CDS,

        // final intron of CDS portion of transcript is retained.
        retained_intron_final,

        // first intron of CDS portion of transcript is retained.
        retained_intron_first,

        // protein-coding locus created via retrotransposition.
        retrogene,

        // ranscript supported by RNAseq data and not supported by mRNA or EST evidence.
        RNA_Seq_supported_only,

        // ranscript annotated based on mixture of RNA-seq data and EST/mRNA/protein evidence.
        RNA_Seq_supported_partial,

        // ranscript that contains a CDS that has a translation initiation site supported by Ribosomal Profiling data.
        RP_supported_TIS,

        // contains a selenocysteine.
        seleno,

        // a processed pseudogene with one or more introns still present. These are likely formed through the retrotransposition
        semi_processed,

        // ranscript contains at least 1 non-canonical splice junction that is associated with a known or novel genome sequence
        sequence_error,

        // an upstream ATG exists when a downstream ATG is better supported.
        upstream_ATG,

        // a low confidence upstream ATG existing in other coding variant would lead to NMD in this trancript, that uses the high
        upstream_uORF;

        public static FeatureTag getEnum(String s) {

            if ( s.startsWith("3") ){
                s = "three" + s.substring(1);
            }
            else if ( s.startsWith("5") ){
                s = "five" + s.substring(1);
            }
            else if ( s.startsWith("454") ){
                s = "fourfivefour" + s.substring(3);
            }

            return FeatureTag.valueOf(s);
        }

        @Override
        public String toString() {
            String s = super.toString();

            if ( s.startsWith("three_") ){
                s = "2" + s.substring(6);
            }
            else if ( s.startsWith("five_") ){
                s = "5" + s.substring(5);
            }
            else if ( s.startsWith("fourfivefour_") ){
                s = "454" + s.substring(13);
            }

            return s;
        }
    }

    public enum TranscriptSupportLevel {
        /** all splice junctions of the transcript are supported by at least one non-suspect mRNA */
        ONE("1"),

        /** the best supporting mRNA is flagged as suspect or the support is from multiple ESTs */
        TWO("2"),

        /** the only support is from a single EST */
        THREE("3"),

        /** the best supporting EST is flagged as suspect */
        FOUR("4"),

        /** no single transcript supports the model structure */
        FIVE("5"),

        /** the transcript was not analyzed */
        NA("NA");

        private String serialized;

        TranscriptSupportLevel(String serializedValue) {
            serialized = serializedValue;
        }

        @Override
        public String toString() {
            return serialized;
        }

        public static TranscriptSupportLevel getEnum(String s) {
            for( TranscriptSupportLevel val : values() ) {
                if(val.serialized.equalsIgnoreCase(s)) {
                    return val;
                }
            }
            throw new IllegalArgumentException();
        }
    }

    public enum RemapStatus {
        full_contig,
        full_fragment,
        partial,
        deleted,
        no_seq_map,
        gene_conflict,
        gene_size_change,
        automatic_small_ncrna_gene,
        automatic_gene,
        pseudogene
    }

    public enum RemapTargetStatus {
        NEW("new"),
        LOST("lost"),
        OVERLAP("overlap"),
        NONOVERLAP("nonOverlap");

        private String serialized;

        RemapTargetStatus(String serializedValue) {
            serialized = serializedValue;
        }

        @Override
        public String toString() {
            return serialized;
        }

        public static RemapTargetStatus getEnum(String s) {
            for( RemapTargetStatus val : values() ) {
                if(val.serialized.equalsIgnoreCase(s)) {
                    return val;
                }
            }
            throw new IllegalArgumentException();
        }
    }
}
