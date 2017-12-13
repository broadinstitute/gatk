package org.broadinstitute.hellbender.utils.codecs.gencode;

/**
 * A Gencode GTF Feature representing a CDS.
 *
 * A GTF Feature represents one row of a GTF File.
 * The specification of a GTF file is defined here:
 * http://mblab.wustl.edu/GTF22.html
 *
 * Created by jonn on 7/25/17.
 */
final public class GencodeGtfCDSFeature extends GencodeGtfFeature {

    private GencodeGtfCDSFeature(final String[] gtfFields) {
        super(gtfFields);
    }

    public static GencodeGtfFeature create(final String[] gtfFields) {
        return new GencodeGtfCDSFeature(gtfFields);
    }

    private GencodeGtfCDSFeature(final GencodeGtfFeatureBaseData baseData) {
        super(baseData);
    }

    public static GencodeGtfFeature create(final GencodeGtfFeatureBaseData baseData) {
        return new GencodeGtfCDSFeature(baseData);
    }
}
