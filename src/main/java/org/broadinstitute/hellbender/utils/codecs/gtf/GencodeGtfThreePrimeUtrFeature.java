package org.broadinstitute.hellbender.utils.codecs.gtf;

/**
 * A Gencode GTF Feature representing a three prime untranslated region.
 *
 * A GTF Feature represents one row of a GTF File.
 * The specification of a GTF file is defined here:
 * http://mblab.wustl.edu/GTF22.html
 *
 * Created by Hailey on 8/9/21.
 */
final public class GencodeGtfThreePrimeUtrFeature extends GencodeGtfFeature {

    private GencodeGtfThreePrimeUtrFeature(final String[] gtfFields, final String gtfFileType) { super(gtfFields, gtfFileType); }

    public static GencodeGtfFeature create(final String[] gtfFields, final String gtfFileType) {
        return new GencodeGtfThreePrimeUtrFeature(gtfFields, gtfFileType);
    }

    private GencodeGtfThreePrimeUtrFeature(final GencodeGtfFeatureBaseData baseData) { super(baseData); }

    public static GencodeGtfFeature create(final GencodeGtfFeatureBaseData baseData) {
        return new GencodeGtfThreePrimeUtrFeature(baseData);
    }
}


