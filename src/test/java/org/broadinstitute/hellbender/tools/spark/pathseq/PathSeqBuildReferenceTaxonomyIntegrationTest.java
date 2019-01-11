package org.broadinstitute.hellbender.tools.spark.pathseq;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.testng.Assert;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;

public class PathSeqBuildReferenceTaxonomyIntegrationTest extends CommandLineProgramTest {

    private File refseqReferenceFile;
    private File genbankReferenceFile;
    private File taxdumpFile;

    @Override
    public String getTestedClassName() {
        return PathSeqBuildReferenceTaxonomy.class.getSimpleName();
    }

    @BeforeTest
    public void before() { 
        refseqReferenceFile = getTestFile("test.fasta");
        genbankReferenceFile = getTestFile("genbank_test.fasta");
        taxdumpFile = getTestFile("taxdump.tar.gz");
    }

    @DataProvider(name = "buildTaxonomyData")
    public Object[][] getTaxonomyBuilderData() {
        return new Object[][]{
                {"catalog.txt.gz", null, 0, "tax.db"},
                {"catalog.txt.gz", null, 50, "tax_min50.db"},
                {null, "genbank.txt.gz", 0, "tax_genbank.db"}
        };
    }

    @Test(dataProvider = "buildTaxonomyData")
    public void test(final String refseqFilename, final String genbankFilename, final int minLength, final String expectedDatabaseFilename) throws Exception {

        final ArgumentsBuilder args = new ArgumentsBuilder();
        if (refseqFilename != null) {
            final File catalogFile = getTestFile(refseqFilename);
            args.addFileArgument(PathSeqBuildReferenceTaxonomy.REFSEQ_CATALOG_LONG_NAME, catalogFile);
            args.addReference(refseqReferenceFile);
        } else if (genbankFilename != null) {
            final File catalogFile = getTestFile(genbankFilename);
            args.addFileArgument(PathSeqBuildReferenceTaxonomy.GENBANK_CATALOG_LONG_NAME, catalogFile);
            args.addReference(genbankReferenceFile);
        }
        args.addFileArgument(PathSeqBuildReferenceTaxonomy.TAX_DUMP_LONG_NAME, taxdumpFile);
        args.addArgument(PathSeqBuildReferenceTaxonomy.MIN_NON_VIRUS_CONTIG_LENGTH_LONG_NAME, Integer.toString(minLength));
        final File outputFile = createTempFile("test_out", ".db");
        args.addOutput(outputFile);

        this.runCommandLine(args.getArgsArray());

        final File expectedFile = getTestFile(expectedDatabaseFilename);

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

    @Test(expectedExceptions = Exception.class)
    public void testBadCatalog() {

        //Try with invalid catalog file
        final File catalog_bad = getSafeNonExistentFile("bad.gz");
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addReference(refseqReferenceFile);
        args.addFileArgument(PathSeqBuildReferenceTaxonomy.REFSEQ_CATALOG_LONG_NAME, catalog_bad);
        args.addFileArgument(PathSeqBuildReferenceTaxonomy.TAX_DUMP_LONG_NAME, taxdumpFile);
        final File outputFile = createTempFile("test_bad_catalog", ".db");
        args.addOutput(outputFile);
        this.runCommandLine(args.getArgsArray());
    }

    @Test(expectedExceptions = Exception.class)
    public void testBadTaxdump() {

        //Try with invalid taxonomy dump file
        final File taxdump_bad = getSafeNonExistentFile("bad.tar.gz");
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addReference(refseqReferenceFile);
        args.addFileArgument(PathSeqBuildReferenceTaxonomy.REFSEQ_CATALOG_LONG_NAME, getTestFile("tax.db"));
        args.addFileArgument(PathSeqBuildReferenceTaxonomy.TAX_DUMP_LONG_NAME, taxdump_bad);
        final File outputFile = createTempFile("test_bad_dump", ".db");
        args.addOutput(outputFile);
        this.runCommandLine(args.getArgsArray());

    }

}