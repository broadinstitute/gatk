package org.broadinstitute.hellbender.tools.gvs.common;

import java.util.Map;

/**
 * Provides GA4GH RefGet accessions (SQ.* digests) for chromosome names.
 *
 * <p>The accession is used when computing VRS allele and location IDs via
 * {@link VrsIdComputer}. It uniquely identifies the reference sequence globally
 * and is the SHA-512t24u digest of the chromosome's full sequence.</p>
 *
 * <p>Currently implemented as static lookup tables for GRCh38/hg38 and GRCh37/hg19.
 * Both "chr"-prefixed (e.g. "chr1") and bare (e.g. "1") contig names are accepted
 * for both builds. If a different implementation is needed in the future (e.g.
 * querying a live RefGet server), replace the body of {@link #getAccession} while
 * keeping this class's interface.</p>
 *
 * <p>Accessions were obtained from the local SeqRepo instance, e.g.:
 * <pre>
 *   curl http://127.0.0.1:5000/seqrepo/1/metadata/refseq:NC_000001.11   # GRCh38 chr1
 *   curl http://127.0.0.1:5000/seqrepo/1/metadata/GRCh37:1              # GRCh37 chr1
 * </pre>
 * and verified against the ga4gh:SQ.* alias in the response.</p>
 */
public class RefGetAccessionProvider {

    public enum GenomeBuild {
        GRCH38, GRCH37;

        /** Accept common string aliases (case-insensitive). */
        public static GenomeBuild fromString(final String s) {
            switch (s.toLowerCase()) {
                case "grch38": case "hg38": case "38": return GRCH38;
                case "grch37": case "hg19": case "37": case "19": return GRCH37;
                default: throw new IllegalArgumentException("Unknown genome build: '" + s + "'");
            }
        }
    }

    /** GRCh38 (hg38) chromosome → GA4GH SQ accession (keyed with "chr" prefix). */
    private static final Map<String, String> GRCH38_ACCESSIONS = Map.ofEntries(
        Map.entry("chr1",  "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO"),
        Map.entry("chr2",  "SQ.pnAqCRBrTsUoBghSD1yp_jXWSmlbdh4g"),
        Map.entry("chr3",  "SQ.Zu7h9AggXxhTaGVsy7h_EZSChSZGcmgX"),
        Map.entry("chr4",  "SQ.HxuclGHh0XCDuF8x6yQrpHUBL7ZntAHc"),
        Map.entry("chr5",  "SQ.aUiQCzCPZ2d0csHbMSbh2NzInhonSXwI"),
        Map.entry("chr6",  "SQ.0iKlIQk2oZLoeOG9P1riRU6hvL5Ux8TV"),
        Map.entry("chr7",  "SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul"),
        Map.entry("chr8",  "SQ.209Z7zJ-mFypBEWLk4rNC6S_OxY5p7bs"),
        Map.entry("chr9",  "SQ.KEO-4XBcm1cxeo_DIQ8_ofqGUkp4iZhI"),
        Map.entry("chr10", "SQ.ss8r_wB0-b9r44TQTMmVTI92884QvBiB"),
        Map.entry("chr11", "SQ.2NkFm8HK88MqeNkCgj78KidCAXgnsfV1"),
        Map.entry("chr12", "SQ.6wlJpONE3oNb4D69ULmEXhqyDZ4vwNfl"),
        Map.entry("chr13", "SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT"),
        Map.entry("chr14", "SQ.eK4D2MosgK_ivBkgi6FVPg5UXs1bYESm"),
        Map.entry("chr15", "SQ.AsXvWL1-2i5U_buw6_niVIxD6zTbAuS6"),
        Map.entry("chr16", "SQ.yC_0RBj3fgBlvgyAuycbzdubtLxq-rE0"),
        Map.entry("chr17", "SQ.dLZ15tNO1Ur0IcGjwc3Sdi_0A6Yf4zm7"),
        Map.entry("chr18", "SQ.vWwFhJ5lQDMhh-czg06YtlWqu0lvFAZV"),
        Map.entry("chr19", "SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl"),
        Map.entry("chr20", "SQ.-A1QmD_MatoqxvgVxBLZTONHz9-c7nQo"),
        Map.entry("chr21", "SQ.5ZUqxCmDDgN4xTRbaSjN8LwgZironmB8"),
        Map.entry("chr22", "SQ.7B7SHsmchAR0dFcDCuSFjJAo7tX87krQ"),
        Map.entry("chrX",  "SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP"),
        Map.entry("chrY",  "SQ.8_liLu1aycC0tPQPFmUaGXJLDs5SbPZ5"),
        Map.entry("chrM",  "SQ.k3grVkjY-hoWcCUojHw6VU6GE3MZ8Sct")
    );

    /**
     * GRCh37 (hg19) chromosome → GA4GH SQ accession (keyed with "chr" prefix).
     * Note: chrM uses the same rCRS sequence as GRCh38 (NC_012920.1 / GRCh37.p13:MT).
     */
    private static final Map<String, String> GRCH37_ACCESSIONS = Map.ofEntries(
        Map.entry("chr1",  "SQ.S_KjnFVz-FE7M0W6yoaUDgYxLPc1jyWU"),
        Map.entry("chr2",  "SQ.9KdcA9ZpY1Cpvxvg8bMSLYDUpsX6GDLO"),
        Map.entry("chr3",  "SQ.VNBualIltAyi2AI_uXcKU7M9XUOuA7MS"),
        Map.entry("chr4",  "SQ.iy7Zfceb5_VGtTQzJ-v5JpPbpeifHD_V"),
        Map.entry("chr5",  "SQ.vbjOdMfHJvTjK_nqvFvpaSKhZillW0SX"),
        Map.entry("chr6",  "SQ.KqaUhJMW3CDjhoVtBetdEKT1n6hM-7Ek"),
        Map.entry("chr7",  "SQ.IW78mgV5Cqf6M24hy52hPjyyo5tCCd86"),
        Map.entry("chr8",  "SQ.tTm7wmhz0G4lpt8wPspcNkAD_qiminj6"),
        Map.entry("chr9",  "SQ.HBckYGQ4wYG9APHLpjoQ9UUe9v7NxExt"),
        Map.entry("chr10", "SQ.-BOZ8Esn8J88qDwNiSEwUr5425UXdiGX"),
        Map.entry("chr11", "SQ.XXi2_O1ly-CCOi3HP5TypAw7LtC6niFG"),
        Map.entry("chr12", "SQ.105bBysLoDFQHhajooTAUyUkNiZ8LJEH"),
        Map.entry("chr13", "SQ.Ewb9qlgTqN6e_XQiRVYpoUfZJHXeiUfH"),
        Map.entry("chr14", "SQ.5Ji6FGEKfejK1U6BMScqrdKJK8GqmIGf"),
        Map.entry("chr15", "SQ.zIMZb3Ft7RdWa5XYq0PxIlezLY2ccCgt"),
        Map.entry("chr16", "SQ.W6wLoIFOn4G7cjopxPxYNk2lcEqhLQFb"),
        Map.entry("chr17", "SQ.AjWXsI7AkTK35XW9pgd3UbjpC3MAevlz"),
        Map.entry("chr18", "SQ.BTj4BDaaHYoPhD3oY2GdwC_l0uqZ92UD"),
        Map.entry("chr19", "SQ.ItRDD47aMoioDCNW_occY5fWKZBKlxCX"),
        Map.entry("chr20", "SQ.iy_UbUrvECxFRX5LPTH_KPojdlT7BKsf"),
        Map.entry("chr21", "SQ.LpTaNW-hwuY_yARP0rtarCnpCQLkgVCg"),
        Map.entry("chr22", "SQ.XOgHwwR3Upfp5sZYk6ZKzvV25a4RBVu8"),
        Map.entry("chrX",  "SQ.v7noePfnNpK8ghYXEqZ9NukMXW7YeNsm"),
        Map.entry("chrY",  "SQ.BT7QyW5iXaX_1PSX-msSGYsqRdMKqkj-"),
        Map.entry("chrM",  "SQ.k3grVkjY-hoWcCUojHw6VU6GE3MZ8Sct")  // same rCRS as GRCh38
    );

    private RefGetAccessionProvider() {}

    /**
     * Normalize any contig name to a "chr"-prefixed key for internal map lookup.
     * Handles: "1"→"chr1", "chr1"→"chr1", "X"→"chrX", "MT"→"chrM", "chrMT"→"chrM".
     */
    private static String normalize(final String contig) {
        String c = contig.startsWith("chr") ? contig : "chr" + contig;
        // GRCh37 uses "MT" for mitochondria; normalize to "chrM"
        if (c.equals("chrMT")) {
            c = "chrM";
        }
        return c;
    }

    /**
     * Returns the GA4GH RefGet accession for the given GRCh38 contig name.
     * Equivalent to {@code getAccession(contig, GenomeBuild.GRCH38)}.
     *
     * @param contig contig name (e.g. "chr20", "20")
     * @return the SQ.* accession string, never null
     * @throws IllegalArgumentException if the contig is not recognized
     */
    public static String getAccession(final String contig) {
        return getAccession(contig, GenomeBuild.GRCH38);
    }

    /**
     * Returns the GA4GH RefGet accession for the given contig and genome build.
     *
     * <p>Both "chr"-prefixed (e.g. "chr1") and bare (e.g. "1") contig names are
     * accepted for both builds. GRCh37's "MT" is normalized to "chrM".</p>
     *
     * @param contig contig name from a VCF or VariantContext
     * @param build  {@link GenomeBuild#GRCH38} or {@link GenomeBuild#GRCH37}
     * @return the SQ.* accession string, never null
     * @throws IllegalArgumentException if the contig is not recognized for the given build
     */
    public static String getAccession(final String contig, final GenomeBuild build) {
        final String key = normalize(contig);
        final Map<String, String> map = build == GenomeBuild.GRCH37 ? GRCH37_ACCESSIONS : GRCH38_ACCESSIONS;
        final String accession = map.get(key);
        if (accession == null) {
            throw new IllegalArgumentException(
                "No " + build + " RefGet accession found for contig: '" + contig + "'. " +
                "Only standard chromosomes (chr1-22, chrX, chrY, chrM) are supported.");
        }
        return accession;
    }

    /**
     * Returns true if the given contig name is known for GRCh38
     * (i.e. {@link #getAccession(String)} will not throw).
     */
    public static boolean isKnownContig(final String contig) {
        return isKnownContig(contig, GenomeBuild.GRCH38);
    }

    /**
     * Returns true if the given contig name is known for the given build
     * (i.e. {@link #getAccession(String, GenomeBuild)} will not throw).
     */
    public static boolean isKnownContig(final String contig, final GenomeBuild build) {
        final String key = normalize(contig);
        final Map<String, String> map = build == GenomeBuild.GRCH37 ? GRCH37_ACCESSIONS : GRCH38_ACCESSIONS;
        return map.containsKey(key);
    }
}

