package org.broadinstitute.hellbender.utils.codecs.gencode;

/**
 * A Gencode GTF Feature representing a start codon.
 *
 * A GTF Feature represents one row of a GTF File.
 * The specification of a GTF file is defined here:
 * http://mblab.wustl.edu/GTF22.html
 *
 * Created by jonn on 7/25/17.
 */
final public class GencodeGtfStartCodonFeature extends GencodeGtfFeature {

    private GencodeGtfStartCodonFeature(final String[] gtfFields) {
        super(gtfFields);
    }

    public static GencodeGtfFeature create(final String[] gtfFields) {
        return new GencodeGtfStartCodonFeature(gtfFields);
    }

    private GencodeGtfStartCodonFeature(final GencodeGtfFeatureBaseData baseData) {
        super(baseData);
    }

    public static GencodeGtfFeature create(final GencodeGtfFeatureBaseData baseData) {
        return new GencodeGtfStartCodonFeature(baseData);
    }
}
