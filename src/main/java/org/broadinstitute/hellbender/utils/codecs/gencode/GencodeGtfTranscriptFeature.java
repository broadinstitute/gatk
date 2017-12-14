package org.broadinstitute.hellbender.utils.codecs.gencode;

import java.util.ArrayList;
import java.util.List;
import java.util.Objects;

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

    private GencodeGtfTranscriptFeature(final String[] gtfFields) {
        super(gtfFields);
    }

    public static GencodeGtfFeature create(final String[] gtfFields) {
        return new GencodeGtfTranscriptFeature(gtfFields);
    }

    private GencodeGtfTranscriptFeature(final GencodeGtfFeatureBaseData baseData) {
        super(baseData);
    }

    public static GencodeGtfFeature create(final GencodeGtfFeatureBaseData baseData) {
        return new GencodeGtfTranscriptFeature(baseData);
    }

    // ================================================================================================

    private final List<GencodeGtfExonFeature>            exons = new ArrayList<>();
    private final List<GencodeGtfSelenocysteineFeature>  selenocysteines = new ArrayList<>();
    private final List<GencodeGtfUTRFeature>             utrs = new ArrayList<>();

    // ================================================================================================

    public List<GencodeGtfExonFeature> getExons() {
        return exons;
    }

    public void addExon( final GencodeGtfExonFeature exon ) {
        exons.add(exon);
    }

    public List<GencodeGtfSelenocysteineFeature> getSelenocysteines() {
        return selenocysteines;
    }

    public void addSelenocysteine( final GencodeGtfSelenocysteineFeature selenocysteine ) {
        selenocysteines.add(selenocysteine);
    }

    public List<GencodeGtfUTRFeature> getUtrs() {
        return utrs;
    }

    public void addUtr( final GencodeGtfUTRFeature utr ) { utrs.add(utr); }

    @Override
    protected List<GencodeGtfFeature> getAllFeatures() {
        final ArrayList<GencodeGtfFeature> list = new ArrayList<>();
        list.add(this);

        for ( final GencodeGtfExonFeature feature : exons ) {
            list.addAll(feature.getAllFeatures());
        }

        for ( final GencodeGtfSelenocysteineFeature feature : selenocysteines ) {
            list.addAll(feature.getAllFeatures());
        }

        for ( final GencodeGtfUTRFeature feature : utrs ) {
            list.addAll(feature.getAllFeatures());
        }

        return list;
    }

    @Override
    public boolean equals(final Object other) {
        if (other == null) {
            return false;
        }
        else if ( this == other ) {
            return true;
        }

        if ( (!(other instanceof GencodeGtfTranscriptFeature)) ) {
            return false;
        }

        final GencodeGtfTranscriptFeature otherTranscript = (GencodeGtfTranscriptFeature) other;

        if ( !super.equals(otherTranscript) ) {
            return false;
        }

        if ( (!Objects.equals(exons, otherTranscript.exons)) ||
             (!Objects.equals(selenocysteines, otherTranscript.selenocysteines)) ||
             (!Objects.equals(utrs, otherTranscript.utrs)) ) {
            return false;
        }

        return true;
    }

    @Override
    public int hashCode() {
        int result = super.hashCode();
        result = 31 * result + (exons != null ? exons.hashCode() : 0);
        result = 31 * result + (selenocysteines != null ? selenocysteines.hashCode() : 0);
        result = 31 * result + (utrs != null ? utrs.hashCode() : 0);
        return result;
    }
}
