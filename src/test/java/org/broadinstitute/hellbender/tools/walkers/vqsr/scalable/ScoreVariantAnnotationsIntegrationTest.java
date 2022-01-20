package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.testng.annotations.Test;

import static org.testng.Assert.*;

public class ScoreVariantAnnotationsIntegrationTest extends CommandLineProgramTest {

    @Test
    public void test1kgp50ExomesSNP() {
        final String[] arguments = {
                "-L", "chr1",
                "-V", "/home/slee/working/vqsr/1kgp-50-exomes/resources/1kgp-50-exomes.sites_only.vcf.gz",
                "-O", "/home/slee/working/vqsr/scalable/score-test/test.snp",
                "--python-script", "/home/slee/working/vqsr/scalable/score-test/isolation-forest-scoring.py",
                "--model-prefix", "/home/slee/working/vqsr/scalable/train-test/test.snp",
                "-an", "FS",
                "-an", "ReadPosRankSum",
                "-an", "MQRankSum",
                "-an", "QD",
                "-an", "SOR",
                "-an", "MQ",
                "--trust-all-polymorphic",
                "-mode", "SNP",
                "--resource:extracted-training,training=true,truth=false", "/home/slee/working/vqsr/scalable/extract-test/test.snp.vcf",
                "--resource:hapmap,training=false,truth=true", "/mnt/4AB658D7B658C4DB/working/ref/hapmap_3.3.hg38.vcf.gz",
                "--resource:omni,training=false,truth=true", "/mnt/4AB658D7B658C4DB/working/ref/1000G_omni2.5.hg38.vcf.gz",
                "--verbosity", "DEBUG"
        };
        runCommandLine(arguments);
    }

    @Test
    public void test1kgp50ExomesIndel() {
        final String[] arguments = {
                "-L", "chr1",
                "-V", "/home/slee/working/vqsr/1kgp-50-exomes/resources/1kgp-50-exomes.sites_only.vcf.gz",
                "-O", "/home/slee/working/vqsr/scalable/score-test/test.indel",
                "--python-script", "/home/slee/working/vqsr/scalable/score-test/isolation-forest-scoring.py",
                "--model-prefix", "/home/slee/working/vqsr/scalable/train-test/test.indel",
                "-an", "FS",
                "-an", "ReadPosRankSum",
                "-an", "MQRankSum",
                "-an", "QD",
                "-an", "SOR",
                "-an", "MQ",
                "--trust-all-polymorphic",
                "-mode", "INDEL",
                "--resource:extracted-training,training=true,truth=false", "/home/slee/working/vqsr/scalable/extract-test/test.indel.vcf",
                "--resource:mills,training=false,truth=true", "/mnt/4AB658D7B658C4DB/working/ref/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
                "--verbosity", "DEBUG"
        };
        runCommandLine(arguments);
    }
}