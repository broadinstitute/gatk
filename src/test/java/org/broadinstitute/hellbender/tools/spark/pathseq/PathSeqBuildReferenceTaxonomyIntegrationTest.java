package org.broadinstitute.hellbender.tools.spark.pathseq;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;

public class PathSeqBuildReferenceTaxonomyIntegrationTest extends CommandLineProgramTest {

    @Override
    public String getTestedClassName() {
        return PathSeqBuildReferenceTaxonomy.class.getSimpleName();
    }

    @Test
    public void test() throws Exception {

        final File expectedFile = getTestFile("tax.db");
        final File reference = getTestFile("test.fasta");
        final File catalog = getTestFile("catalog.txt.gz");
        final File taxdump = getTestFile("taxdump.tar.gz");
        final File output = createTempFile("test", ".db");

        ArgumentsBuilder args = new ArgumentsBuilder();
        args.addReference(reference);
        args.addFileArgument("refseqCatalogPath", catalog);
        args.addFileArgument("taxdumpPath", taxdump);
        args.addOutput(output);
        this.runCommandLine(args.getArgsArray());

        try {
            final Kryo kryo = new Kryo();
            kryo.setReferences(false);
            Input input = new Input(new FileInputStream(expectedFile));
            final PSTree tree_expected = (PSTree) kryo.readClassAndObject(input);

            input = new Input(new FileInputStream(output));
            final PSTree tree_actual = (PSTree) kryo.readClassAndObject(input);
            Assert.assertEquals(tree_actual, tree_expected);
        } catch (FileNotFoundException e) {
            throw new IOException("Error with Kryo IO", e);
        }

        final File genbank = getTestFile("genbank.txt.gz");

        args = new ArgumentsBuilder();
        args.addReference(reference);
        args.addFileArgument("genbankCatalogPath", genbank);
        args.addFileArgument("taxdumpPath", taxdump);
        args.addOutput(output);
        this.runCommandLine(args.getArgsArray());

        try {
            final Kryo kryo = new Kryo();
            kryo.setReferences(false);
            Input input = new Input(new FileInputStream(expectedFile));
            final PSTree tree_expected = (PSTree) kryo.readClassAndObject(input);

            input = new Input(new FileInputStream(output));
            final PSTree tree_actual = (PSTree) kryo.readClassAndObject(input);
            Assert.assertEquals(tree_actual, tree_expected);
        } catch (FileNotFoundException e) {
            throw new IOException("Error with Kryo IO", e);
        }

        //Try with invalid catalog file
        final File catalog_bad = getTestFile("catalog_that_does_not_exist.gz");
        args = new ArgumentsBuilder();
        args.addReference(reference);
        args.addFileArgument("refseqCatalogPath", catalog_bad);
        args.addFileArgument("taxdumpPath", taxdump);
        args.addOutput(output);
        try {
            this.runCommandLine(args.getArgsArray());
            Assert.fail("Did not throw UserException for non-existent catalog file");
        } catch (UserException e) {
        }

        //Try with invalid taxonomy dump file
        final File taxdump_bad = getTestFile("taxdump_that_does_not_exist.tar.gz");
        args = new ArgumentsBuilder();
        args.addReference(reference);
        args.addFileArgument("refseqCatalogPath", catalog);
        args.addFileArgument("taxdumpPath", taxdump_bad);
        args.addOutput(output);
        try {
            this.runCommandLine(args.getArgsArray());
            Assert.fail("Did not throw UserException for non-existent taxonomy dump file");
        } catch (UserException e) {
        }

    }

}