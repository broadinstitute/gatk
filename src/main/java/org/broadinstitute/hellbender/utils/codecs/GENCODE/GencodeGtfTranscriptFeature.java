package org.broadinstitute.hellbender.utils.codecs.GENCODE;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

/**
 * A Gencode GTF Feature representing a transcript.
 *
 * A GTF Feature represents one row of a GTF File.
 * The specification of a GTF file is defined here:
 * http://mblab.wustl.edu/GTF22.html
 *
 * Created by jonn on 7/25/17.
 */
final public class GencodeGtfTranscriptFeature extends GencodeGtfFeature {

    private GencodeGtfTranscriptFeature(String[] gtfFields) {
        super(gtfFields);
    }

    public static GencodeGtfFeature create(String[] gtfFields) {
        return new GencodeGtfTranscriptFeature(gtfFields);
    }

    private GencodeGtfTranscriptFeature(long featureOrderNumber,
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

        super(featureOrderNumber, chromosomeName, annotationSource, featureType, genomicStartLocation, genomicEndLocation, genomicStrand, genomicPhase, geneId, transcriptId, geneType, geneStatus, geneName, transcriptType, transcriptStatus, transcriptName, exonNumber, exonId, locusLevel, optionalFields, anonymousOptionalFields);
    }

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

        return new GencodeGtfTranscriptFeature(featureOrderNumber, chromosomeName, annotationSource, featureType, genomicStartLocation, genomicEndLocation, genomicStrand, genomicPhase, geneId, transcriptId, geneType, geneStatus, geneName, transcriptType, transcriptStatus, transcriptName, exonNumber, exonId, locusLevel, optionalFields, anonymousOptionalFields);
    }

    // ================================================================================================

    private ArrayList<GencodeGtfExonFeature>            exons = new ArrayList<>();
    private ArrayList<GencodeGtfSelenocysteineFeature>  selenocysteines = new ArrayList<>();
    private ArrayList<GencodeGtfUTRFeature>             utrs = new ArrayList<>();

    // ================================================================================================

    public ArrayList<GencodeGtfExonFeature> getExons() {
        return exons;
    }

    public void addExon( GencodeGtfExonFeature exon ) {
        exons.add(exon);
    }

    public ArrayList<GencodeGtfSelenocysteineFeature> getSelenocysteines() {
        return selenocysteines;
    }

    public void addSelenocysteine( GencodeGtfSelenocysteineFeature selenocysteine ) {
        selenocysteines.add(selenocysteine);
    }

    public ArrayList<GencodeGtfUTRFeature> getUtrs() {
        return utrs;
    }

    public void addUtr( GencodeGtfUTRFeature utr ) { utrs.add(utr); }

    @Override
    protected List<GencodeGtfFeature> getAllFeatures() {
        ArrayList<GencodeGtfFeature> list = new ArrayList<>();
        list.add(this);

        for ( GencodeGtfExonFeature feature : exons ) {
            list.addAll(feature.getAllFeatures());
        }

        for ( GencodeGtfSelenocysteineFeature feature : selenocysteines ) {
            list.addAll(feature.getAllFeatures());
        }

        for ( GencodeGtfUTRFeature feature : utrs ) {
            list.addAll(feature.getAllFeatures());
        }

        return list;
    }

    @Override
    public boolean equals(Object other) {
        if ( (!(other instanceof GencodeGtfTranscriptFeature)) ) {
            return false;
        }

        GencodeGtfTranscriptFeature otherTranscript = (GencodeGtfTranscriptFeature) other;

        if ( (!super.equals(otherTranscript)) ||
             (!exons.equals(otherTranscript.exons)) ||
             (!selenocysteines.equals(otherTranscript.selenocysteines)) ||
             (!utrs.equals(otherTranscript.utrs))) {
            return false;
        }

        return true;
    }
}
