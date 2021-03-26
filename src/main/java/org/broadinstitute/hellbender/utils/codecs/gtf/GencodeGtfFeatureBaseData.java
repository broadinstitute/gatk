package org.broadinstitute.hellbender.utils.codecs.gtf;

import htsjdk.tribble.annotation.Strand;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.ArrayList;
import java.util.List;
import java.util.Objects;

/**
 * Struct-like container class for the fields in a {@link GencodeGtfFeature}
 *
 * This is designed to be a basic dummy class to make feature instantiation easier.
 *
 * Created by jonn on 8/11/17.
 */
final public class GencodeGtfFeatureBaseData {

    /**
     * The relative order of this Feature.
     * Normally it is a line number indicating the position in the original data file of this feature.
     */
    public int                                     featureOrderNumber      = GencodeGtfFeature.NO_FEATURE_ORDER;

    /**
     * The source file type from which the data in this {@link GencodeGtfFeatureBaseData} originated.
     * This should always correspond to a value returned by {@link AbstractGtfCodec#getGtfFileType()}.
     */
    public String gtfSourceFileType;

    /**
     * Location of this feature on the genome.
     */
    public SimpleInterval                          genomicPosition;

    /**
     * Keyword identifying the source of this feature.
     *
     * There is no restriction on which keywords can be used here, but some
     * known annotation sources are defined in their GTF parent objects:
     *      {@link GencodeGtfFeature#ANNOTATION_SOURCE_ENSEMBL}
     *      {@link GencodeGtfFeature#ANNOTATION_SOURCE_HAVANA}
     *      {@link GencodeGtfFeature#ANNOTATION_SOURCE_ENA}
     */
    public String annotationSource;

    /**
     * Type of this feature.
     */
    public GencodeGtfFeature.FeatureType           featureType;

    //TODO: Make this a Strand, not a GenomicStrand.
    /**
     * Which strand this feature is on.
     */
    public Strand                                  genomicStrand;

    /**
     * Frame/phase of this feature.
     */
    public GencodeGtfFeature.GenomicPhase          genomicPhase;

    // "Required" GENCODE GTF Fields:
    public String                                   geneId                  = null;
    public String                                   transcriptId            = null;
    public GencodeGtfFeature.GeneTranscriptStatus   geneStatus              = null;
    /**
     * There are no formal definitions for what can be a valid geneType and the
     * number of possible values seem to increase as new versions of Gencode are released.
     */
    public String                                   geneType                = null;
    public String                                   geneName                = null;
    /**
     * There are no formal definitions for what can be a valid transcriptType and the
     * number of possible values seem to increase as new versions of Gencode are released.
     */
    public String                                   transcriptType          = null;
    public GencodeGtfFeature.GeneTranscriptStatus   transcriptStatus        = null;
    public String                                   transcriptName          = null;
    public int                                      exonNumber              = GencodeGtfFeature.NO_EXON_NUMBER;
    public String                                   exonId                  = null;
    public GencodeGtfFeature.LocusLevel             locusLevel              = null;

    /**
     * Optional GENCODE GTF Fields.
     * For details, see the following:
     *     https://www.gencodegenes.org/data_format.html
     *     https://www.gencodegenes.org/gencode_tags.html
     */
    public List<GencodeGtfFeature.OptionalField<?>> optionalFields          = new ArrayList<>();

    /**
     * Additional optional GTF fields.
     */
    public String                                   anonymousOptionalFields = null;

    public GencodeGtfFeatureBaseData() {}

    public GencodeGtfFeatureBaseData(
            final String gtfSourceFileType,
            final int featureOrderNumber,
            final String chromosomeName,
            final String annotationSource,
            final GencodeGtfFeature.FeatureType featureType,
            final int genomicStartLocation,
            final int genomicEndLocation,
            final Strand genomicStrand,
            final GencodeGtfFeature.GenomicPhase genomicPhase,
            final String geneId,
            final String transcriptId,
            final String geneType,
            final GencodeGtfFeature.GeneTranscriptStatus geneStatus,
            final String geneName,
            final String transcriptType,
            final GencodeGtfFeature.GeneTranscriptStatus transcriptStatus,
            final String transcriptName,
            final int exonNumber,
            final String exonId,
            final GencodeGtfFeature.LocusLevel locusLevel,
            final List<GencodeGtfFeature.OptionalField<?>> optionalFields,
            final String anonymousOptionalFields
    ) {
        this.gtfSourceFileType = gtfSourceFileType;
        this.featureOrderNumber = featureOrderNumber;

        this.genomicPosition = new SimpleInterval(chromosomeName, genomicStartLocation, genomicEndLocation);

        this.annotationSource = annotationSource;
        this.featureType = featureType;
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

    @Override
    public boolean equals(final Object that) {

        if (that == null) {
            return false;
        }
        else if ( this == that ) {
            return true;
        }

        boolean isEqual = that instanceof GencodeGtfFeatureBaseData;

        if ( isEqual ) {

            final GencodeGtfFeatureBaseData thatBaseData = (GencodeGtfFeatureBaseData) that;

            isEqual =
                    Objects.equals(gtfSourceFileType,       thatBaseData.gtfSourceFileType)         &&
                    Objects.equals(featureOrderNumber,      thatBaseData.featureOrderNumber)        &&
                    Objects.equals(genomicPosition,         thatBaseData.genomicPosition)           &&
                    Objects.equals(annotationSource,        thatBaseData.annotationSource)          &&
                    Objects.equals(featureType,             thatBaseData.featureType)               &&
                    Objects.equals(genomicStrand,           thatBaseData.genomicStrand)             &&
                    Objects.equals(genomicPhase,            thatBaseData.genomicPhase)              &&
                    Objects.equals(geneId,                  thatBaseData.geneId)                    &&
                    Objects.equals(transcriptId,            thatBaseData.transcriptId)              &&
                    Objects.equals(geneType,                thatBaseData.geneType)                  &&
                    Objects.equals(geneStatus,              thatBaseData.geneStatus)                &&
                    Objects.equals(geneName,                thatBaseData.geneName)                  &&
                    Objects.equals(transcriptType,          thatBaseData.transcriptType)            &&
                    Objects.equals(transcriptStatus,        thatBaseData.transcriptStatus)          &&
                    Objects.equals(transcriptName,          thatBaseData.transcriptName)            &&
                    Objects.equals(exonNumber,              thatBaseData.exonNumber)                &&
                    Objects.equals(exonId,                  thatBaseData.exonId)                    &&
                    Objects.equals(locusLevel,              thatBaseData.locusLevel)                &&
                    Objects.equals(anonymousOptionalFields, thatBaseData.anonymousOptionalFields)   &&
                    Objects.equals(optionalFields,          thatBaseData.optionalFields);
        }

        return isEqual;
    }

    @Override
    public int hashCode() {
        int result = featureOrderNumber;
        result = 31 * result + (gtfSourceFileType != null ? gtfSourceFileType.hashCode() : 0);
        result = 31 * result + (genomicPosition != null ? genomicPosition.hashCode() : 0);
        result = 31 * result + (annotationSource != null ? annotationSource.hashCode() : 0);
        result = 31 * result + (featureType != null ? featureType.hashCode() : 0);
        result = 31 * result + (genomicStrand != null ? genomicStrand.hashCode() : 0);
        result = 31 * result + (genomicPhase != null ? genomicPhase.hashCode() : 0);
        result = 31 * result + (geneId != null ? geneId.hashCode() : 0);
        result = 31 * result + (transcriptId != null ? transcriptId.hashCode() : 0);
        result = 31 * result + (geneStatus != null ? geneStatus.hashCode() : 0);
        result = 31 * result + (geneType != null ? geneType.hashCode() : 0);
        result = 31 * result + (geneName != null ? geneName.hashCode() : 0);
        result = 31 * result + (transcriptType != null ? transcriptType.hashCode() : 0);
        result = 31 * result + (transcriptStatus != null ? transcriptStatus.hashCode() : 0);
        result = 31 * result + (transcriptName != null ? transcriptName.hashCode() : 0);
        result = 31 * result + exonNumber;
        result = 31 * result + (exonId != null ? exonId.hashCode() : 0);
        result = 31 * result + (locusLevel != null ? locusLevel.hashCode() : 0);
        result = 31 * result + (optionalFields != null ? optionalFields.hashCode() : 0);
        result = 31 * result + (anonymousOptionalFields != null ? anonymousOptionalFields.hashCode() : 0);
        return result;
    }
}
