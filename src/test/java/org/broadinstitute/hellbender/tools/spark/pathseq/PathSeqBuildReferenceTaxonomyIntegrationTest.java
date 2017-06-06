package org.broadinstitute.hellbender.tools.spark.pathseq;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.testng.Assert;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;

public class PathSeqBuildReferenceTaxonomyIntegrationTest extends CommandLineProgramTest {

    private File expectedFile;
    private File referenceFile;
    private File catalogFile;
    private File taxdumpFile;
    private File outputFile;

    @Override
    public String getTestedClassName() {
        return PathSeqBuildReferenceTaxonomy.class.getSimpleName();
    }

    @BeforeTest
    public void before() {
        expectedFile = getTestFile("tax.db");
        referenceFile = getTestFile("test.fasta");
        catalogFile = getTestFile("catalog.txt.gz");
        taxdumpFile = getTestFile("taxdump.tar.gz");
        outputFile = createTempFile("test", ".db");
    }


    @Test
    public void test() throws Exception {

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addReference(referenceFile);
        args.addFileArgument("refseqCatalogPath", catalogFile);
        args.addFileArgument("taxdumpPath", taxdumpFile);
        args.addOutput(outputFile);
        this.runCommandLine(args.getArgsArray());

        try {
            final Kryo kryo = new Kryo();
            kryo.setReferences(false);
            Input input = new Input(new FileInputStream(expectedFile));
            final PSTaxonomyDatabase dbExpected = kryo.readObject(input, PSTaxonomyDatabase.class);

            input = new Input(new FileInputStream(outputFile));
            final PSTaxonomyDatabase dbTest = kryo.readObject(input, PSTaxonomyDatabase.class);
            Assert.assertEquals(dbTest.tree, dbExpected.tree);
            Assert.assertEquals(dbTest.accessionToTaxId, dbExpected.accessionToTaxId);
        } catch (FileNotFoundException e) {
            throw new IOException("Error with Kryo IO", e);
        }
    }

    @Test
    public void testGenBank() {

        //Test with Genbank catalog file
        final File genbank = getTestFile("genbank.txt.gz");
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addReference(referenceFile);
        args.addFileArgument("genbankCatalogPath", genbank);
        args.addFileArgument("taxdumpPath", taxdumpFile);
        args.addOutput(outputFile);
        this.runCommandLine(args.getArgsArray());
    }

    @Test(expectedExceptions = Exception.class)
    public void testBadCatalog() {

        //Try with invalid catalog file
        final File catalog_bad = getSafeNonExistentFile("bad.gz");
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addReference(referenceFile);
        args.addFileArgument("refseqCatalogPath", catalog_bad);
        args.addFileArgument("taxdumpPath", taxdumpFile);
        args.addOutput(outputFile);
        this.runCommandLine(args.getArgsArray());
    }

    @Test(expectedExceptions = Exception.class)
    public void testBadTaxdump() {

        //Try with invalid taxonomy dump file
        final File taxdump_bad = getSafeNonExistentFile("bad.tar.gz");
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addReference(referenceFile);
        args.addFileArgument("refseqCatalogPath", catalogFile);
        args.addFileArgument("taxdumpPath", taxdump_bad);
        args.addOutput(outputFile);
        this.runCommandLine(args.getArgsArray());

    }

}