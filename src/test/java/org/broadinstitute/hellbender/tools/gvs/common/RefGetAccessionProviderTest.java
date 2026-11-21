package org.broadinstitute.hellbender.tools.gvs.common;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

public class RefGetAccessionProviderTest {

    // -------------------------------------------------------------------------
    // GRCh38 known-good values: digests of the RefSeq sequences (see RefGetAccessionProvider)
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
    // Contig names must match the hg38 FASTA exactly, so no other build's names resolve
    // -------------------------------------------------------------------------

    @DataProvider(name = "nonHg38ContigNames")
    public Object[][] nonHg38ContigNames() {
        return new Object[][] { {"1"}, {"20"}, {"X"}, {"MT"}, {"chrMT"} };
    }

    @Test(dataProvider = "nonHg38ContigNames")
    public void testNonHg38ContigNameIsUnknown(final String contig) {
        Assert.assertFalse(RefGetAccessionProvider.isKnownContig(contig), contig);
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

    // -------------------------------------------------------------------------
    // isKnownContig
    // -------------------------------------------------------------------------

    @Test
    public void testIsKnownContigTrue() {
        Assert.assertTrue(RefGetAccessionProvider.isKnownContig("chr20"));
        Assert.assertTrue(RefGetAccessionProvider.isKnownContig("chrM"));
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
    public void testUnknownContigThrows() {
        RefGetAccessionProvider.getAccession("chrUn_gl000220");
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testEmptyContigThrows() {
        RefGetAccessionProvider.getAccession("");
    }
}
