package org.broadinstitute.hellbender.utils.codecs.gtf;

/**
 * A Gencode GTF Feature representing a five prime untranslated region.
 *
 * A GTF Feature represents one row of a GTF File.
 * The specification of a GTF file is defined here:
 * http://mblab.wustl.edu/GTF22.html
 *
 * Created by Hailey on 8/9/21.
 */
final public class GencodeGtfFivePrimeUtrFeature extends GencodeGtfFeature {

    private GencodeGtfFivePrimeUtrFeature(final String[] gtfFields, final String gtfFileType) { super(gtfFields, gtfFileType); }

    public static GencodeGtfFeature create(final String[] gtfFields, final String gtfFileType) {
        return new GencodeGtfFivePrimeUtrFeature(gtfFields, gtfFileType);
    }

    private GencodeGtfFivePrimeUtrFeature(final GencodeGtfFeatureBaseData baseData) { super(baseData); }

    public static GencodeGtfFeature create(final GencodeGtfFeatureBaseData baseData) {
        return new GencodeGtfFivePrimeUtrFeature(baseData);
    }
}


