package org.broadinstitute.hellbender.tools.gvs.common;

import java.util.Map;

/**
 * Provides GA4GH RefGet accessions (SQ.* digests) for chromosome names.
 *
 * <p>The accession is used when computing VRS allele and location IDs via
 * {@link VrsIdComputer}. It uniquely identifies the reference sequence globally
 * and is the SHA-512t24u digest of the chromosome's full sequence.</p>
 *
 * <p>Currently implemented as a static lookup table for GRCh38, the only reference GVS supports. Contig names
 * must match the hg38 FASTA exactly ("chr1" to "chr22", "chrX", "chrY", "chrM"); bare names such as "1" are not
 * recognized, so a GRCh37 contig can never be given a GRCh38 accession. If a different implementation is needed
 * in the future (e.g. querying a live RefGet server), replace the body of {@link #getAccession} while keeping
 * this class's interface.</p>
 *
 * <p>Each accession is the digest of the chromosome's RefSeq sequence: NC_000001.11 for
 * chr1, and so on. These are the accessions VRS tools
 * such as vrs-python use. To re-derive one, download the sequence from NCBI, uppercase it,
 * and base64url-encode the first 24 bytes of its SHA-512:
 * <pre>
 *   curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&amp;id=NC_000021.9&amp;rettype=fasta"
 * </pre>
 * They are not digests of the GVS reference FASTA, which is the GRCh38 analysis set. The two
 * differ wherever the analysis set N-masks sequence, as on chr21.</p>
 */
public class RefGetAccessionProvider {

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


    private RefGetAccessionProvider() {}

    /**
     * Returns the GA4GH RefGet accession for the given GRCh38 contig name.
     *
     * @param contig contig name as in the hg38 FASTA (e.g. "chr20")
     * @return the SQ.* accession string, never null
     * @throws IllegalArgumentException if the contig is not recognized
     */
    public static String getAccession(final String contig) {
        final String accession = GRCH38_ACCESSIONS.get(contig);
        if (accession == null) {
            throw new IllegalArgumentException(
                "No GRCh38 RefGet accession found for contig: '" + contig + "'. " +
                "Only standard chromosomes (chr1-22, chrX, chrY, chrM) are supported.");
        }
        return accession;
    }

    /**
     * Returns true if the given contig name is known for GRCh38
     * (i.e. {@link #getAccession(String)} will not throw).
     */
    public static boolean isKnownContig(final String contig) {
        return GRCH38_ACCESSIONS.containsKey(contig);
    }
}
