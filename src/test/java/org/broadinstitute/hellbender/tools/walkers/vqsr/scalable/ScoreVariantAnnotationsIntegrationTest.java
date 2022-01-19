package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable;

import org.broadinstitute.hellbender.cmdline.CommandLineProgramIntegrationTest;
import org.testng.annotations.Test;

import static org.testng.Assert.*;

public class ScoreVariantAnnotationsIntegrationTest extends CommandLineProgramIntegrationTest {
    @Test
    public void test1kgp50Exomes() {
        final String[] arguments = {
                "-L", "chr1",
                "-V", "/home/slee/working/vqsr/1kgp-50-exomes/resources/1kgp-50-exomes.sites_only.vcf.gz",
                "-O", "/home/slee/working/vqsr/scalable/score-test/test",
                "--python-script", "/home/slee/working/vqsr/scalable/score-test/isolation-forest-scoring.py",
                "--model-prefix", "/home/slee/working/vqsr/scalable/train-test/test",
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
}