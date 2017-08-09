package org.broadinstitute.hellbender.utils.codecs.GENCODE;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Objects;
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

    private GencodeGtfExonFeature(final GencodeGtfFeatureBaseData baseData) {
        super(baseData);
    }

    public static GencodeGtfFeature create(final GencodeGtfFeatureBaseData baseData) {
        return new GencodeGtfExonFeature(baseData);
    }

    // ============================================================================================================


    private GencodeGtfCDSFeature                        cds = null;
    private GencodeGtfStartCodonFeature                 startCodon = null;
    private GencodeGtfStopCodonFeature                  stopCodon = null;

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
        if (this.cds != null) {
            throw new RuntimeException("Attempting to set cds, but it already contains a value!");
        }
        this.cds = cds;
    }

    public void setStartCodon(GencodeGtfStartCodonFeature startCodon) {
        if (this.startCodon != null) {
            throw new RuntimeException("Attempting to set cds, but it already contains a value!");
        }
        this.startCodon = startCodon;
    }

    public void setStopCodon(GencodeGtfStopCodonFeature stopCodon) {
        if (this.stopCodon != null) {
            throw new RuntimeException("Attempting to set cds, but it already contains a value!");
        }
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
        if (other == null) {
            return false;
        }
        else if ( this == other ) {
            return true;
        }

        if (!(other instanceof GencodeGtfExonFeature)) {
            return false;
        }

        GencodeGtfExonFeature otherExon = (GencodeGtfExonFeature) other;

        if ( !super.equals(otherExon) ) {
            return false;
        }

        if ( (!Objects.equals(cds, otherExon.cds)) ||
             (!Objects.equals(startCodon, otherExon.startCodon)) ||
             (!Objects.equals(stopCodon, otherExon.stopCodon)) ) {
            return false;
        }

        return true;
    }

    @Override
    public int hashCode() {
        int result = super.hashCode();
        result = 31 * result + (cds != null ? cds.hashCode() : 0);
        result = 31 * result + (startCodon != null ? startCodon.hashCode() : 0);
        result = 31 * result + (stopCodon != null ? stopCodon.hashCode() : 0);
        return result;
    }
}
