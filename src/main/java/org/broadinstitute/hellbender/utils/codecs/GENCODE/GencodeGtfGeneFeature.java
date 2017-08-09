package org.broadinstitute.hellbender.utils.codecs.GENCODE;

import java.util.ArrayList;
import java.util.List;
import java.util.Objects;

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

    private GencodeGtfGeneFeature(final GencodeGtfFeatureBaseData baseData) {
        super(baseData);
    }

    public static GencodeGtfFeature create(final GencodeGtfFeatureBaseData baseData) {
        return new GencodeGtfGeneFeature(baseData);
    }

    // ================================================================================================

    private List<GencodeGtfTranscriptFeature> transcripts = new ArrayList<>();

    // ================================================================================================

    public void addTranscript(GencodeGtfTranscriptFeature transcript) { transcripts.add(transcript); }

    public List<GencodeGtfTranscriptFeature> getTranscripts() {
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

        if (other == null) {
            return false;
        }
        else if ( this == other ) {
            return true;
        }

        if ( (!(other instanceof GencodeGtfGeneFeature)) ) {
            return false;
        }

        GencodeGtfGeneFeature otherGene = (GencodeGtfGeneFeature) other;

        if ( (!super.equals(otherGene)) ) {
            return false;
        }

        if ( (!Objects.equals(transcripts, otherGene.transcripts)) ) {
            return false;
        }

        return true;
    }

    @Override
    public int hashCode() {
        int result = super.hashCode();
        result = 31 * result + (transcripts != null ? transcripts.hashCode() : 0);
        return result;
    }
}
