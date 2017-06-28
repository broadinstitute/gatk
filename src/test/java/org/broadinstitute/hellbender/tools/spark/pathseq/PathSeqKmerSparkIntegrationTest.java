package org.broadinstitute.hellbender.tools.spark.pathseq;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVKmer;
import org.broadinstitute.hellbender.tools.spark.utils.HopscotchSet;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.testng.annotations.Test;
import org.apache.commons.io.FileUtils;
import org.testng.Assert;
import java.io.File;
import java.util.Set;

public class PathSeqKmerSparkIntegrationTest extends CommandLineProgramTest {

    @Override
    public String getTestedClassName() {
        return PathSeqKmerSpark.class.getSimpleName();
    }

    @SuppressWarnings("unchecked")
    @Test(groups = "spark")
    public void test() throws Exception {
        final File expectedFile = getTestFile("kmer.hss");
        final File ref = getTestFile("hg19mini.fasta");
        final File output = createTempFile("test", ".hss");
        if (!output.delete()) {
            Assert.fail();
        }
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addFileArgument("reference", ref);
        args.addOutput(output);
        this.runCommandLine(args.getArgsArray());

        Input input_expected = new Input(FileUtils.openInputStream(expectedFile));
        Input input_test = new Input(FileUtils.openInputStream(output));

        Kryo kryo=new Kryo();
        kryo.setReferences(false);

        Set<SVKmer> expectedKmerLib = (HopscotchSet<SVKmer>)kryo.readClassAndObject(input_expected);
        Set<SVKmer> testKmerLib = (HopscotchSet<SVKmer>)kryo.readClassAndObject(input_test);

        Assert.assertEquals(expectedKmerLib,testKmerLib);
    }

}
