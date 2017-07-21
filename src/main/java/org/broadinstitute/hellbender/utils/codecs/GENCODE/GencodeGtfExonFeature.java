package org.broadinstitute.hellbender.utils.codecs.GENCODE;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

/**
 * A Gencode GTF Feature representing an exon.
 *
 * A GTF Feature represents one row of a GTF File.
 * The specification of a GTF file is defined here:
 * http://mblab.wustl.edu/GTF22.html
 *
 * Created by jonn on 7/25/17.
 */
final public class GencodeGtfExonFeature extends GencodeGtfFeature {

    private GencodeGtfExonFeature(String[] gtfFields) {
        super(gtfFields);
    }

    public static GencodeGtfFeature create(String[] gtfFields) {
        return new GencodeGtfExonFeature(gtfFields);
    }

    private GencodeGtfExonFeature(long featureOrderNumber,
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

        return new GencodeGtfExonFeature(featureOrderNumber, chromosomeName, annotationSource, featureType, genomicStartLocation, genomicEndLocation, genomicStrand, genomicPhase, geneId, transcriptId, geneType, geneStatus, geneName, transcriptType, transcriptStatus, transcriptName, exonNumber, exonId, locusLevel, optionalFields, anonymousOptionalFields);
    }

    // ============================================================================================================


    GencodeGtfCDSFeature                        cds = null;
    GencodeGtfStartCodonFeature                 startCodon = null;
    GencodeGtfStopCodonFeature                  stopCodon = null;

    // ============================================================================================================

    public GencodeGtfCDSFeature getCds() {
        return cds;
    }

    public GencodeGtfStartCodonFeature getStartCodon() {
        return startCodon;
    }

    public GencodeGtfStopCodonFeature getStopCodon() {
        return stopCodon;
    }

    // ============================================================================================================

    public void setCds(GencodeGtfCDSFeature cds) {
        this.cds = cds;
    }

    public void setStartCodon(GencodeGtfStartCodonFeature startCodon) {
        this.startCodon = startCodon;
    }

    public void setStopCodon(GencodeGtfStopCodonFeature stopCodon) {
        this.stopCodon = stopCodon;
    }


    @Override
    protected List<GencodeGtfFeature> getAllFeatures() {
        ArrayList<GencodeGtfFeature> list = new ArrayList<>();
        list.add(this);

        if ( cds != null ) { list.add(cds) ; }
        if ( startCodon != null ) { list.add(startCodon) ; }
        if ( stopCodon != null ) { list.add(stopCodon) ; }

        return list;
    }

    @Override
    public boolean equals(Object other) {
        if (!(other instanceof GencodeGtfExonFeature)) {
            return false;
        }

        GencodeGtfExonFeature otherExon = (GencodeGtfExonFeature) other;

        if ( !super.equals(otherExon) ) {
            return false;
        }

        if (!(((cds == null) && (otherExon.cds == null)) || ((cds != null) && cds.equals(otherExon.cds))) ||
            !(((startCodon == null) && (otherExon.startCodon == null)) || ((startCodon != null) && startCodon.equals(otherExon.startCodon))) ||
            !(((stopCodon == null) && (otherExon.stopCodon == null)) || ((stopCodon != null) && stopCodon.equals(otherExon.stopCodon))) ) {
            return false;
        }

        return true;
    }
}
