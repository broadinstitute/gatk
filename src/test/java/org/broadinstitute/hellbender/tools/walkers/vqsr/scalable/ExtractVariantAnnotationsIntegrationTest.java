package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable;

import org.broadinstitute.hdf5.HDF5File;
import org.broadinstitute.hdf5.HDF5LibException;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.tools.copynumber.utils.HDF5Utils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;

public final class ExtractVariantAnnotationsIntegrationTest extends CommandLineProgramTest {

    @Test
    public void test1kgp50ExomesSNPASExactMatch() {
        final File outputDir = createTempDir("extract-test");
        final String[] arguments = {
                "-L", "chr1",
                "-V", largeFileTestDir + "VQSR/chr1snippet.doctoredMQ.sites_only.vcf.gz",
                "-O", new File(outputDir, "test.snp.as").getAbsolutePath(),
                "--maximum-chunk-size", "100000",
                "--use-allele-specific-annotations",
                "-an", "AS_FS",
                "-an", "AS_ReadPosRankSum",
                "-an", "AS_MQRankSum",
                "-an", "AS_QD",
                "-an", "AS_SOR",
                "-an", "AS_MQ",
                "--trust-all-polymorphic",
                "-mode", "SNP",
                "--resource:hapmap,training=true,truth=true", "/mnt/4AB658D7B658C4DB/working/ref/hapmap_3.3.hg38.vcf.gz",
                "--resource:omni,training=true,truth=true", "/mnt/4AB658D7B658C4DB/working/ref/1000G_omni2.5.hg38.vcf.gz",
                "--resource:1000G,training=true,truth=false", "/mnt/4AB658D7B658C4DB/working/ref/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
                "--verbosity", "DEBUG"
        };
        runCommandLine(arguments);
        assertAnnotationFilesEqual(new File(outputDir, "test.snp.as.annot.hdf5"), new File("/home/slee/working/vqsr/scalable/extract-exact-match", "test.snp.as.annot.hdf5"));
    }

    @Test
    public void test1kgp50ExomesSNPExactMatch() {
        final File outputDir = createTempDir("extract-test");
        final String[] arguments = {
                "-L", "chr1",
                "-V", "/home/slee/working/vqsr/1kgp-50-exomes/resources/1kgp-50-exomes.sites_only.vcf.gz",
                "-O", new File(outputDir, "test.snp").getAbsolutePath(),
                "--maximum-chunk-size", "100000",
                "-an", "FS",
                "-an", "ReadPosRankSum",
                "-an", "MQRankSum",
                "-an", "QD",
                "-an", "SOR",
                "-an", "MQ",
                "--trust-all-polymorphic",
                "-mode", "SNP",
                "--resource:hapmap,training=true,truth=true", "/mnt/4AB658D7B658C4DB/working/ref/hapmap_3.3.hg38.vcf.gz",
                "--resource:omni,training=true,truth=true", "/mnt/4AB658D7B658C4DB/working/ref/1000G_omni2.5.hg38.vcf.gz",
                "--resource:1000G,training=true,truth=false", "/mnt/4AB658D7B658C4DB/working/ref/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
                "--verbosity", "DEBUG"
        };
        runCommandLine(arguments);
        assertAnnotationFilesEqual(new File(outputDir, "test.snp.annot.hdf5"), new File("/home/slee/working/vqsr/scalable/extract-exact-match", "test.snp.annot.hdf5"));
    }

    private static void assertAnnotationFilesEqual(final File file1,
                                                   final File file2) {
        try (final HDF5File annotationsHDF5File1 = new HDF5File(file1, HDF5File.OpenMode.READ_ONLY);
             final HDF5File annotationsHDF5File2 = new HDF5File(file2, HDF5File.OpenMode.READ_ONLY)) {
            Assert.assertEquals(
                    annotationsHDF5File1.readStringArray("/data/annotation_names"),
                    annotationsHDF5File2.readStringArray("/data/annotation_names"));
            for (final String label : Arrays.asList("is_training", "is_truth")) {
                Assert.assertEquals(
                        annotationsHDF5File1.readDoubleArray(String.format("/data/%s", label)),
                        annotationsHDF5File2.readDoubleArray(String.format("/data/%s", label)));
            }
            Assert.assertEquals(
                    HDF5Utils.readChunkedDoubleMatrix(annotationsHDF5File1, "/data/annotations"),
                    HDF5Utils.readChunkedDoubleMatrix(annotationsHDF5File2, "/data/annotations"));
        } catch (final HDF5LibException exception) {
            Assert.fail("Exception encountered during reading of annotations:", exception);
        }
    }

    @Test
    public void test1kgp50ExomesSNP() {
        final String[] arguments = {
                "-L", "chr1",
                "-V", "/home/slee/working/vqsr/1kgp-50-exomes/resources/1kgp-50-exomes.sites_only.vcf.gz",
                "-O", "/home/slee/working/vqsr/scalable/extract-test/test.snp",
                "--maximum-chunk-size", "100000",
                "-an", "FS",
                "-an", "ReadPosRankSum",
                "-an", "MQRankSum",
                "-an", "QD",
                "-an", "SOR",
                "-an", "MQ",
                "--trust-all-polymorphic",
                "-mode", "SNP",
                "--resource:hapmap,training=true,truth=true", "/mnt/4AB658D7B658C4DB/working/ref/hapmap_3.3.hg38.vcf.gz",
                "--resource:omni,training=true,truth=true", "/mnt/4AB658D7B658C4DB/working/ref/1000G_omni2.5.hg38.vcf.gz",
                "--resource:1000G,training=true,truth=false", "/mnt/4AB658D7B658C4DB/working/ref/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
                "--verbosity", "DEBUG"
        };
        runCommandLine(arguments);
    }

    @Test
    public void test1kgp50ExomesIndel() {
        final String[] arguments = {
                "-L", "chr1",
                "-V", "/home/slee/working/vqsr/1kgp-50-exomes/resources/1kgp-50-exomes.sites_only.vcf.gz",
                "-O", "/home/slee/working/vqsr/scalable/extract-test/test.indel",
                "-an", "FS",
                "-an", "ReadPosRankSum",
                "-an", "MQRankSum",
                "-an", "QD",
                "-an", "SOR",
                "-an", "MQ",
                "--trust-all-polymorphic",
                "-mode", "INDEL",
                "--resource:mills,training=true,truth=true", "/mnt/4AB658D7B658C4DB/working/ref/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
                "--verbosity", "DEBUG"
        };
        runCommandLine(arguments);
    }

    @Test
    public void testJbxSNP() {
        final String[] arguments = {
                "-V", "/home/slee/working/vqsr/scalable/jbx/resources/Test50Callset.annoated_pids.sites-only.vcf.gz",
                "-O", "/home/slee/working/vqsr/scalable/jbx/Test50Callset.snp.extract",
                "-an", "FS",
                "-an", "ReadPosRankSum",
                "-an", "MQRankSum",
                "-an", "QD",
                "-an", "SOR",
                "-an", "MQ",
                "-an", "COMBINED_TREE_SCORE",
                "--trust-all-polymorphic",
                "-mode", "SNP",
                "--resource:hapmap,training=true,truth=true", "/mnt/4AB658D7B658C4DB/working/ref/hapmap_3.3.hg38.vcf.gz",
                "--resource:omni,training=true,truth=true", "/mnt/4AB658D7B658C4DB/working/ref/1000G_omni2.5.hg38.vcf.gz",
                "--resource:1000G,training=true,truth=false", "/mnt/4AB658D7B658C4DB/working/ref/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
                "--verbosity", "DEBUG"
        };
        runCommandLine(arguments);
    }
}