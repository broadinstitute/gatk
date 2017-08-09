package org.broadinstitute.hellbender.utils.codecs.GENCODE;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Objects;
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

    private GencodeGtfTranscriptFeature(final GencodeGtfFeatureBaseData baseData) {
        super(baseData);
    }

    public static GencodeGtfFeature create(final GencodeGtfFeatureBaseData baseData) {
        return new GencodeGtfTranscriptFeature(baseData);
    }

    // ================================================================================================

    private List<GencodeGtfExonFeature>            exons = new ArrayList<>();
    private List<GencodeGtfSelenocysteineFeature>  selenocysteines = new ArrayList<>();
    private List<GencodeGtfUTRFeature>             utrs = new ArrayList<>();

    // ================================================================================================

    public List<GencodeGtfExonFeature> getExons() {
        return exons;
    }

    public void addExon( GencodeGtfExonFeature exon ) {
        exons.add(exon);
    }

    public List<GencodeGtfSelenocysteineFeature> getSelenocysteines() {
        return selenocysteines;
    }

    public void addSelenocysteine( GencodeGtfSelenocysteineFeature selenocysteine ) {
        selenocysteines.add(selenocysteine);
    }

    public List<GencodeGtfUTRFeature> getUtrs() {
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
        if (other == null) {
            return false;
        }
        else if ( this == other ) {
            return true;
        }

        if ( (!(other instanceof GencodeGtfTranscriptFeature)) ) {
            return false;
        }

        GencodeGtfTranscriptFeature otherTranscript = (GencodeGtfTranscriptFeature) other;

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
