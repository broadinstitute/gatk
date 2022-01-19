package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.testng.annotations.Test;

import static org.testng.Assert.*;

public class TrainVariantAnnotationModelIntegrationTest extends CommandLineProgramTest {

    @Test
    public void test1kgp50ExomesSNP() {
        final String[] arguments = {
                "--annotations-hdf5", "/home/slee/working/vqsr/scalable/extract-test/test.snp.annot.hdf5",
                "-O", "/home/slee/working/vqsr/scalable/train-test/test.snp",
                "--python-script", "/home/slee/working/vqsr/scalable/train-test/isolation-forest.py",
                "--hyperparameters-json", "/home/slee/working/vqsr/scalable/train-test/hyperparameters.json",
                "-mode", "SNP",
                "--verbosity", "DEBUG"
        };
        runCommandLine(arguments);
    }

    @Test
    public void test1kgp50ExomesIndel() {
        final String[] arguments = {
                "--annotations-hdf5", "/home/slee/working/vqsr/scalable/extract-test/test.indel.annot.hdf5",
                "-O", "/home/slee/working/vqsr/scalable/train-test/test.indel",
                "--python-script", "/home/slee/working/vqsr/scalable/train-test/isolation-forest.py",
                "--hyperparameters-json", "/home/slee/working/vqsr/scalable/train-test/hyperparameters.json",
                "-mode", "INDEL",
                "--verbosity", "DEBUG"
        };
        runCommandLine(arguments);
    }
}