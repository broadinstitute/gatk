package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;

public class TrainVariantAnnotationsModelIntegrationTest extends CommandLineProgramTest {

    private static final String PYTHON_SCRIPT = packageMainResourcesDir + "tools/walkers/vqsr/scalable/isolation-forest.py";

    @Test
    public void test1kgp50ExomesAll() {
        final String[] arguments = {
                "--annotations-hdf5", "/home/slee/working/vqsr/scalable/extract-test/test.all.annot.hdf5",
                "-O", "/home/slee/working/vqsr/scalable/train-test/test.all",
                "--python-script", PYTHON_SCRIPT,
                "--hyperparameters-json", "/home/slee/working/vqsr/scalable/train-test/hyperparameters.json",
                "--mode", "SNP",
                "--mode", "INDEL",
                "--verbosity", "INFO"
        };
        runCommandLine(arguments);
    }

    @Test
    public void test1kgp50ExomesAllUnlabeled() {
        final String[] arguments = {
                "--annotations-hdf5", "/home/slee/working/vqsr/scalable/extract-test/test.all-unlabeled.annot.hdf5",
                "--unlabeled-annotations-hdf5", "/home/slee/working/vqsr/scalable/extract-test/test.all-unlabeled.unlabeled.annot.hdf5",
                "--truth-sensitivity-threshold", "0.95",
                "-O", "/home/slee/working/vqsr/scalable/train-test/test.all-unlabeled",
                "--python-script", PYTHON_SCRIPT,
                "--hyperparameters-json", "/home/slee/working/vqsr/scalable/train-test/hyperparameters.json",
                "--mode", "SNP",
                "--mode", "INDEL",
                "--verbosity", "INFO"
        };
        runCommandLine(arguments);

        runSystemCommand("h5diff /home/slee/working/vqsr/scalable/train-test/test.all-unlabeled.snp.trainingScores.hdf5 /home/slee/working/vqsr/scalable/train-test/expected/test.all-unlabeled.snp.trainingScores.hdf5");
        runSystemCommand("h5diff /home/slee/working/vqsr/scalable/train-test/test.all-unlabeled.snp.truthScores.hdf5 /home/slee/working/vqsr/scalable/train-test/expected/test.all-unlabeled.snp.truthScores.hdf5");
        runSystemCommand("h5diff /home/slee/working/vqsr/scalable/train-test/test.all-unlabeled.snp.unlabeledScores.hdf5 /home/slee/working/vqsr/scalable/train-test/expected/test.all-unlabeled.snp.unlabeledScores.hdf5");
        runSystemCommand("diff /home/slee/working/vqsr/scalable/train-test/test.all-unlabeled.snp.scorer.pkl /home/slee/working/vqsr/scalable/train-test/expected/test.all-unlabeled.snp.scorer.pkl");
        runSystemCommand("diff /home/slee/working/vqsr/scalable/train-test/test.all-unlabeled.snp.negative.scorer.pkl /home/slee/working/vqsr/scalable/train-test/expected/test.all-unlabeled.snp.negative.scorer.pkl");
        runSystemCommand("h5diff /home/slee/working/vqsr/scalable/train-test/test.all-unlabeled.indel.trainingScores.hdf5 /home/slee/working/vqsr/scalable/train-test/expected/test.all-unlabeled.indel.trainingScores.hdf5");
        runSystemCommand("h5diff /home/slee/working/vqsr/scalable/train-test/test.all-unlabeled.indel.truthScores.hdf5 /home/slee/working/vqsr/scalable/train-test/expected/test.all-unlabeled.indel.truthScores.hdf5");
        runSystemCommand("h5diff /home/slee/working/vqsr/scalable/train-test/test.all-unlabeled.indel.unlabeledScores.hdf5 /home/slee/working/vqsr/scalable/train-test/expected/test.all-unlabeled.indel.unlabeledScores.hdf5");
        runSystemCommand("diff /home/slee/working/vqsr/scalable/train-test/test.all-unlabeled.indel.scorer.pkl /home/slee/working/vqsr/scalable/train-test/expected/test.all-unlabeled.indel.scorer.pkl");
        runSystemCommand("diff /home/slee/working/vqsr/scalable/train-test/test.all-unlabeled.indel.negative.scorer.pkl /home/slee/working/vqsr/scalable/train-test/expected/test.all-unlabeled.indel.negative.scorer.pkl");
    }

    private static void runSystemCommand(final String command) {
        try {
            Process process = Runtime.getRuntime().exec(command);

            BufferedReader reader = new BufferedReader(
                    new InputStreamReader(process.getInputStream()));
            while (reader.readLine() != null) {
                Assert.fail(command);
            }

            reader.close();

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    @Test
    public void test1kgp50ExomesSNP() {
        final String[] arguments = {
                "--annotations-hdf5", "/home/slee/working/vqsr/scalable/extract-test/test.snp.annot.hdf5",
                "-O", "/home/slee/working/vqsr/scalable/train-test/test",
                "--python-script", PYTHON_SCRIPT,
                "--hyperparameters-json", "/home/slee/working/vqsr/scalable/train-test/hyperparameters.json",
                "--mode", "SNP",
                "--verbosity", "INFO"
        };
        runCommandLine(arguments);
    }

    @Test
    public void test1kgp50ExomesIndel() {
        final String[] arguments = {
                "--annotations-hdf5", "/home/slee/working/vqsr/scalable/extract-test/test.indel.annot.hdf5",
                "-O", "/home/slee/working/vqsr/scalable/train-test/test",
                "--python-script", PYTHON_SCRIPT,
                "--hyperparameters-json", "/home/slee/working/vqsr/scalable/train-test/hyperparameters.json",
                "--mode", "INDEL",
                "--verbosity", "INFO"
        };
        runCommandLine(arguments);
    }

    @Test
    public void testSNPAS() {
        final String[] arguments = {
                "--annotations-hdf5", "/home/slee/working/vqsr/scalable/extract-test/test.snp.as.annot.hdf5",
                "-O", "/home/slee/working/vqsr/scalable/train-test/test.snp.as",
                "--python-script", PYTHON_SCRIPT,
                "--hyperparameters-json", "/home/slee/working/vqsr/scalable/train-test/hyperparameters.json",
                "--mode", "SNP",
                "--verbosity", "INFO"
        };
        runCommandLine(arguments);
    }

    @Test
    public void test1kgp50ExomesBGMMAll() {
        final String[] arguments = {
                "--annotations-hdf5", "/home/slee/working/vqsr/scalable/extract-test/test.all.annot.hdf5",
                "-O", "/home/slee/working/vqsr/scalable/train-test/test.bgmm.all",
                "--hyperparameters-json", "/home/slee/working/vqsr/scalable/train-test/bgmm-hyperparameters.json",
                "--mode", "SNP",
                "--mode", "INDEL",
                "--verbosity", "INFO"
        };
        runCommandLine(arguments);
    }

    @Test
    public void test1kgp50ExomesBGMMSNP() {
        final String[] arguments = {
                "--annotations-hdf5", "/home/slee/working/vqsr/scalable/extract-test/test.snp.annot.hdf5",
                "-O", "/home/slee/working/vqsr/scalable/train-test/test.bgmm",
                "--hyperparameters-json", "/home/slee/working/vqsr/scalable/train-test/bgmm-hyperparameters.json",
                "--mode", "SNP",
                "--verbosity", "INFO"
        };
        runCommandLine(arguments);
    }

    @Test
    public void test1kgp50ExomesBGMMIndel() {
        final String[] arguments = {
                "--annotations-hdf5", "/home/slee/working/vqsr/scalable/extract-test/test.indel.annot.hdf5",
                "-O", "/home/slee/working/vqsr/scalable/train-test/test.bgmm",
                "--hyperparameters-json", "/home/slee/working/vqsr/scalable/train-test/bgmm-hyperparameters.json",
                "--mode", "INDEL",
                "--verbosity", "INFO"
        };
        runCommandLine(arguments);
    }

    @Test
    public void testJbxAll() {
        final String[] arguments = {
                "--annotations-hdf5", "/home/slee/working/vqsr/scalable/jbx/Test50Callset.all.extract.annot.hdf5",
                "-O", "/home/slee/working/vqsr/scalable/jbx/Test50Callset.all.train",
                "--python-script", PYTHON_SCRIPT,
                "--hyperparameters-json", "/home/slee/working/vqsr/scalable/jbx/hyperparameters.json",
                "--mode", "SNP",
                "--mode", "INDEL",
                "--verbosity", "INFO"
        };
        runCommandLine(arguments);
    }

    @Test
    public void testJbxAllUnlabeled() {
        final String[] arguments = {
                "--annotations-hdf5", "/home/slee/working/vqsr/scalable/jbx/Test50Callset.all-unlabeled.extract.annot.hdf5",
                "--unlabeled-annotations-hdf5", "/home/slee/working/vqsr/scalable/jbx/Test50Callset.all-unlabeled.extract.unlabeled.annot.hdf5",
                "--truth-sensitivity-threshold", "0.95",
                "-O", "/home/slee/working/vqsr/scalable/jbx/Test50Callset.all-unlabeled.train",
                "--python-script", PYTHON_SCRIPT,
                "--hyperparameters-json", "/home/slee/working/vqsr/scalable/jbx/hyperparameters.json",
                "--mode", "SNP",
                "--mode", "INDEL",
                "--verbosity", "INFO"
        };
        runCommandLine(arguments);
    }

    @Test
    public void testJbxSNP() {
        final String[] arguments = {
                "--annotations-hdf5", "/home/slee/working/vqsr/scalable/jbx/Test50Callset.snp.extract.annot.hdf5",
                "-O", "/home/slee/working/vqsr/scalable/jbx/Test50Callset.snp.train",
                "--python-script", PYTHON_SCRIPT,
                "--hyperparameters-json", "/home/slee/working/vqsr/scalable/jbx/hyperparameters.json",
                "--mode", "SNP",
                "--verbosity", "INFO"
        };
        runCommandLine(arguments);
    }

    @Test
    public void testJbxBGMMAll() {
        final String[] arguments = {
                "--annotations-hdf5", "/home/slee/working/vqsr/scalable/jbx/Test50Callset.all.extract.annot.hdf5",
                "-O", "/home/slee/working/vqsr/scalable/jbx/Test50Callset.bgmm.all.train",
                "--hyperparameters-json", "/home/slee/working/vqsr/scalable/train-test/bgmm-hyperparameters.json",
                "--mode", "SNP",
                "--mode", "INDEL",
                "--verbosity", "INFO"
        };
        runCommandLine(arguments);
    }
}