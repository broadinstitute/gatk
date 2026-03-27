package org.broadinstitute.hellbender.tools.gvs.common;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

public class RefGetAccessionProviderTest {

    // -------------------------------------------------------------------------
    // GRCh38 known-good values verified against local seqrepo instance
    // -------------------------------------------------------------------------

    @Test
    public void testChr1GRCh38() {
        Assert.assertEquals(RefGetAccessionProvider.getAccession("chr1"),
            "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO");
    }

    @Test
    public void testChr19GRCh38() {
        // This value is also used as the constant in VrsIdComputerTest
        Assert.assertEquals(RefGetAccessionProvider.getAccession("chr19"),
            "SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl");
    }

    @Test
    public void testChr7GRCh38() {
        Assert.assertEquals(RefGetAccessionProvider.getAccession("chr7"),
            "SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul");
    }

    @Test
    public void testChrX_GRCh38() {
        Assert.assertEquals(RefGetAccessionProvider.getAccession("chrX"),
            "SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP");
    }

    @Test
    public void testChrY_GRCh38() {
        Assert.assertEquals(RefGetAccessionProvider.getAccession("chrY"),
            "SQ.8_liLu1aycC0tPQPFmUaGXJLDs5SbPZ5");
    }

    @Test
    public void testChrM_GRCh38() {
        Assert.assertEquals(RefGetAccessionProvider.getAccession("chrM"),
            "SQ.k3grVkjY-hoWcCUojHw6VU6GE3MZ8Sct");
    }

    // -------------------------------------------------------------------------
    // GRCh37 known-good values verified against local seqrepo instance
    // -------------------------------------------------------------------------

    @Test
    public void testChr1GRCh37() {
        Assert.assertEquals(
            RefGetAccessionProvider.getAccession("chr1", RefGetAccessionProvider.GenomeBuild.GRCH37),
            "SQ.S_KjnFVz-FE7M0W6yoaUDgYxLPc1jyWU");
    }

    @Test
    public void testChr19GRCh37() {
        Assert.assertEquals(
            RefGetAccessionProvider.getAccession("chr19", RefGetAccessionProvider.GenomeBuild.GRCH37),
            "SQ.ItRDD47aMoioDCNW_occY5fWKZBKlxCX");
    }

    @Test
    public void testChrM_GRCh37_SameAsGRCh38() {
        // Both builds share the same rCRS mitochondrial sequence (NC_012920.1)
        final String grch38 = RefGetAccessionProvider.getAccession("chrM", RefGetAccessionProvider.GenomeBuild.GRCH38);
        final String grch37 = RefGetAccessionProvider.getAccession("chrM", RefGetAccessionProvider.GenomeBuild.GRCH37);
        Assert.assertEquals(grch37, grch38,
            "chrM should have the same SQ accession in both builds (rCRS)");
        Assert.assertEquals(grch37, "SQ.k3grVkjY-hoWcCUojHw6VU6GE3MZ8Sct");
    }

    @Test
    public void testGRCh37DiffersFromGRCh38() {
        // Autosomes must have different accessions between builds
        final String g38 = RefGetAccessionProvider.getAccession("chr1");
        final String g37 = RefGetAccessionProvider.getAccession("chr1", RefGetAccessionProvider.GenomeBuild.GRCH37);
        Assert.assertNotEquals(g37, g38, "chr1 should differ between GRCh37 and GRCh38");
    }

    // -------------------------------------------------------------------------
    // GenomeBuild.fromString
    // -------------------------------------------------------------------------

    @Test
    public void testGenomeBuildFromString() {
        Assert.assertEquals(RefGetAccessionProvider.GenomeBuild.fromString("GRCh38"),
            RefGetAccessionProvider.GenomeBuild.GRCH38);
        Assert.assertEquals(RefGetAccessionProvider.GenomeBuild.fromString("hg38"),
            RefGetAccessionProvider.GenomeBuild.GRCH38);
        Assert.assertEquals(RefGetAccessionProvider.GenomeBuild.fromString("38"),
            RefGetAccessionProvider.GenomeBuild.GRCH38);
        Assert.assertEquals(RefGetAccessionProvider.GenomeBuild.fromString("GRCh37"),
            RefGetAccessionProvider.GenomeBuild.GRCH37);
        Assert.assertEquals(RefGetAccessionProvider.GenomeBuild.fromString("hg19"),
            RefGetAccessionProvider.GenomeBuild.GRCH37);
        Assert.assertEquals(RefGetAccessionProvider.GenomeBuild.fromString("37"),
            RefGetAccessionProvider.GenomeBuild.GRCH37);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testGenomeBuildFromStringUnknownThrows() {
        RefGetAccessionProvider.GenomeBuild.fromString("hg17");
    }

    // -------------------------------------------------------------------------
    // Contig normalization: bare names and GRCh37 "MT" alias
    // -------------------------------------------------------------------------

    @Test
    public void testBareContigName_GRCh38() {
        Assert.assertEquals(RefGetAccessionProvider.getAccession("1"),
            RefGetAccessionProvider.getAccession("chr1"));
    }

    @Test
    public void testBareContigName_GRCh37() {
        Assert.assertEquals(
            RefGetAccessionProvider.getAccession("1", RefGetAccessionProvider.GenomeBuild.GRCH37),
            RefGetAccessionProvider.getAccession("chr1", RefGetAccessionProvider.GenomeBuild.GRCH37));
    }

    @Test
    public void testMtAlias_GRCh37() {
        // GRCh37 VCFs may use "MT" instead of "chrM"
        final String viaMT  = RefGetAccessionProvider.getAccession("MT",   RefGetAccessionProvider.GenomeBuild.GRCH37);
        final String viaChrM = RefGetAccessionProvider.getAccession("chrM", RefGetAccessionProvider.GenomeBuild.GRCH37);
        Assert.assertEquals(viaMT, viaChrM);
    }

    @Test
    public void testChrMtAlias() {
        // "chrMT" should also normalize to chrM
        final String viaChrMT = RefGetAccessionProvider.getAccession("chrMT", RefGetAccessionProvider.GenomeBuild.GRCH37);
        final String viaChrM  = RefGetAccessionProvider.getAccession("chrM",  RefGetAccessionProvider.GenomeBuild.GRCH37);
        Assert.assertEquals(viaChrMT, viaChrM);
    }

    // -------------------------------------------------------------------------
    // All 25 chromosomes valid for both builds
    // -------------------------------------------------------------------------

    @DataProvider(name = "allChromosomes")
    public Object[][] allChromosomes() {
        return new Object[][] {
            {"chr1"}, {"chr2"}, {"chr3"}, {"chr4"}, {"chr5"},
            {"chr6"}, {"chr7"}, {"chr8"}, {"chr9"}, {"chr10"},
            {"chr11"}, {"chr12"}, {"chr13"}, {"chr14"}, {"chr15"},
            {"chr16"}, {"chr17"}, {"chr18"}, {"chr19"}, {"chr20"},
            {"chr21"}, {"chr22"}, {"chrX"}, {"chrY"}, {"chrM"}
        };
    }

    @Test(dataProvider = "allChromosomes")
    public void testAllChromosomesReturnValidAccession_GRCh38(final String contig) {
        final String accession = RefGetAccessionProvider.getAccession(contig);
        Assert.assertNotNull(accession);
        Assert.assertTrue(accession.startsWith("SQ."),
            "Accession should start with 'SQ.' for " + contig);
        Assert.assertEquals(accession.length(), 35,
            "SQ.* accession should be 35 chars for " + contig);
    }

    @Test(dataProvider = "allChromosomes")
    public void testAllChromosomesReturnValidAccession_GRCh37(final String contig) {
        final String accession = RefGetAccessionProvider.getAccession(contig, RefGetAccessionProvider.GenomeBuild.GRCH37);
        Assert.assertNotNull(accession);
        Assert.assertTrue(accession.startsWith("SQ."),
            "Accession should start with 'SQ.' for GRCh37 " + contig);
        Assert.assertEquals(accession.length(), 35,
            "SQ.* accession should be 35 chars for GRCh37 " + contig);
    }

    @Test
    public void testAllChromosomesAreUnique_GRCh38() {
        final String[] contigs = {
            "chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
            "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20",
            "chr21","chr22","chrX","chrY","chrM"
        };
        final java.util.Set<String> seen = new java.util.HashSet<>();
        for (final String contig : contigs) {
            Assert.assertTrue(seen.add(RefGetAccessionProvider.getAccession(contig)),
                "Duplicate GRCh38 accession for " + contig);
        }
        Assert.assertEquals(seen.size(), 25);
    }

    @Test
    public void testAllChromosomesAreUnique_GRCh37() {
        final String[] contigs = {
            "chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
            "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20",
            "chr21","chr22","chrX","chrY"  // chrM shares rCRS with GRCh38, excluded here
        };
        final java.util.Set<String> seen = new java.util.HashSet<>();
        for (final String contig : contigs) {
            Assert.assertTrue(
                seen.add(RefGetAccessionProvider.getAccession(contig, RefGetAccessionProvider.GenomeBuild.GRCH37)),
                "Duplicate GRCh37 accession for " + contig);
        }
        Assert.assertEquals(seen.size(), 24);
    }

    // -------------------------------------------------------------------------
    // isKnownContig
    // -------------------------------------------------------------------------

    @Test
    public void testIsKnownContigTrue() {
        Assert.assertTrue(RefGetAccessionProvider.isKnownContig("chr20"));
        Assert.assertTrue(RefGetAccessionProvider.isKnownContig("20"));
        Assert.assertTrue(RefGetAccessionProvider.isKnownContig("chrM"));
        Assert.assertTrue(RefGetAccessionProvider.isKnownContig("chr20", RefGetAccessionProvider.GenomeBuild.GRCH37));
        Assert.assertTrue(RefGetAccessionProvider.isKnownContig("20",    RefGetAccessionProvider.GenomeBuild.GRCH37));
        Assert.assertTrue(RefGetAccessionProvider.isKnownContig("MT",    RefGetAccessionProvider.GenomeBuild.GRCH37));
    }

    @Test
    public void testIsKnownContigFalse() {
        Assert.assertFalse(RefGetAccessionProvider.isKnownContig("chrUn_gl000220"));
        Assert.assertFalse(RefGetAccessionProvider.isKnownContig(""));
    }

    // -------------------------------------------------------------------------
    // Error handling
    // -------------------------------------------------------------------------

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testUnknownContigThrows_GRCh38() {
        RefGetAccessionProvider.getAccession("chrUn_gl000220");
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testUnknownContigThrows_GRCh37() {
        RefGetAccessionProvider.getAccession("chrUn_gl000220", RefGetAccessionProvider.GenomeBuild.GRCH37);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testEmptyContigThrows() {
        RefGetAccessionProvider.getAccession("");
    }
}
