package org.broadinstitute.hellbender.utils.codecs.gencode;

/**
 * A Gencode GTF Feature representing a stop codon.
 *
 * A GTF Feature represents one row of a GTF File.
 * The specification of a GTF file is defined here:
 * http://mblab.wustl.edu/GTF22.html
 *
 * Created by jonn on 7/25/17.
 */
final public class GencodeGtfStopCodonFeature extends GencodeGtfFeature {

    private GencodeGtfStopCodonFeature(final String[] gtfFields) {
        super(gtfFields);
    }

    public static GencodeGtfFeature create(final String[] gtfFields) {
        return new GencodeGtfStopCodonFeature(gtfFields);
    }

    private GencodeGtfStopCodonFeature(final GencodeGtfFeatureBaseData baseData) {
        super(baseData);
    }

    public static GencodeGtfFeature create(final GencodeGtfFeatureBaseData baseData) {
        return new GencodeGtfStopCodonFeature(baseData);
    }
}
