package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable;

import org.broadinstitute.hellbender.cmdline.CommandLineProgramIntegrationTest;
import org.testng.annotations.Test;

import static org.testng.Assert.*;

public class TrainVariantAnnotationModelIntegrationTest extends CommandLineProgramIntegrationTest {

    @Test
    public void test1kgp50Exomes() {
        final String[] arguments = {
                "--annotations-hdf5", "/home/slee/working/vqsr/scalable/extract-test/test.annot.hdf5",
                "-O", "/home/slee/working/vqsr/scalable/train-test/test",
                "--python-script", "/home/slee/working/vqsr/scalable/train-test/isolation-forest.py",
                "--hyperparameters-json", "/home/slee/working/vqsr/scalable/train-test/hyperparameters.json",
                "--verbosity", "DEBUG"
        };
        runCommandLine(arguments);
    }
}