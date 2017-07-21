package org.broadinstitute.hellbender.utils.codecs.GENCODE;

import java.util.ArrayList;
import java.util.List;

/**
 * A Gencode GTF Feature representing a gene.
 *
 * A GTF Feature represents one row of a GTF File.
 * The specification of a GTF file is defined here:
 * http://mblab.wustl.edu/GTF22.html
 *
 * Created by jonn on 7/25/17.
 */
final public class GencodeGtfGeneFeature extends GencodeGtfFeature {

    private GencodeGtfGeneFeature(String[] gtfFields) {
        super(gtfFields);
    }

    public static GencodeGtfFeature create(String[] gtfFields) {
        return new GencodeGtfGeneFeature(gtfFields);
    }

    private GencodeGtfGeneFeature(long lineNumber,
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

        super(lineNumber, chromosomeName, annotationSource, featureType, genomicStartLocation, genomicEndLocation, genomicStrand, genomicPhase, geneId, transcriptId, geneType, geneStatus, geneName, transcriptType, transcriptStatus, transcriptName, exonNumber, exonId, locusLevel, optionalFields, anonymousOptionalFields);
    }

    public static GencodeGtfFeature create(long lineNumber,
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

        return new GencodeGtfGeneFeature(lineNumber, chromosomeName, annotationSource, featureType, genomicStartLocation, genomicEndLocation, genomicStrand, genomicPhase, geneId, transcriptId, geneType, geneStatus, geneName, transcriptType, transcriptStatus, transcriptName, exonNumber, exonId, locusLevel, optionalFields, anonymousOptionalFields);
    }

    // ================================================================================================

    private ArrayList<GencodeGtfTranscriptFeature> transcripts = new ArrayList<>();

    // ================================================================================================

    public void addTranscript(GencodeGtfTranscriptFeature transcript) { transcripts.add(transcript); }

    public ArrayList<GencodeGtfTranscriptFeature> getTranscripts() {
        return transcripts;
    }

    @Override
    protected List<GencodeGtfFeature> getAllFeatures() {
        ArrayList<GencodeGtfFeature> list = new ArrayList<>();
        list.add(this);

        for ( GencodeGtfTranscriptFeature transcript : transcripts ) {
            list.addAll(transcript.getAllFeatures());
        }

        return list;
    }

    @Override
    public boolean equals(Object other) {
        if ( (!(other instanceof GencodeGtfGeneFeature)) ) {
            return false;
        }

        GencodeGtfGeneFeature otherGene = (GencodeGtfGeneFeature) other;

        if ( !super.equals(otherGene) ) {
            return false;
        }

        for ( int i = 0 ; i < transcripts.size(); ++i ) {
            if ( !transcripts.get(i).equals(otherGene.transcripts.get(i))) {
                return false;
            }
        }


        return true;
    }

}
