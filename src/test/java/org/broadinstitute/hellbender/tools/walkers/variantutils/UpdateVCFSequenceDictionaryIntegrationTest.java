package org.broadinstitute.hellbender.tools.walkers.variantutils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileNotFoundException;

public class UpdateVCFSequenceDictionaryIntegrationTest extends CommandLineProgramTest {
    private File testDir = new File(getTestDataDir(), "walkers/variantutils/UpdateVCFSequenceDictionary");

    @DataProvider(name="UpdateGoodSequenceDictionaryData")
    public Object[][] updateGoodSequenceDictionaryData() {
        return new Object[][]{
                new Object[]{ new File(testDir, "variantsNoDict.vcf"), new File(testDir, "variantsWithDict.vcf"), null, false},
                // pass a reference as a reference
                new Object[]{ new File(testDir, "variantsNoDict.vcf"), null, new File(testDir, "exampleFASTA.fasta"), false},
                // pass a reference as a source; we need to test both to ensure the user can bypass the framework sequence
                // dictionary validation that will occur if you use -R
                new Object[]{ new File(testDir, "variantsNoDict.vcf"), new File(testDir, "exampleFASTA.fasta"), null, false},
                new Object[]{ new File(testDir, "variantsNoDict.vcf"), new File(testDir, "exampleFASTA.dict"), null, false},
                // can't handle CRAM - see https://github.com/samtools/htsjdk/issues/731
                new Object[]{ new File(testDir, "variantsNoDict.vcf"), new File(testDir, "exampleBAM.bam"), null, false},
                // already has a dictionary - but force a replace
                new Object[]{ new File(testDir, "variantsWithDict.vcf"), new File(testDir, "exampleFASTA.dict"), null, true},
                // already has a dictionary - but force a replace
                new Object[]{ new File(testDir, "variantsWithDict.vcf"), new File(testDir, "exampleFASTA.dict"), null, true},
        };
    }

    @Test(dataProvider="UpdateGoodSequenceDictionaryData")
    private void testGoodUpdateSequenceDictionary(
            final File inputVariantsFile,
            final File inputSourceFile,
            final File inputReferenceFile,
            final boolean replace) throws FileNotFoundException {
        final SAMSequenceDictionary resultingDictionary =
                updateSequenceDictionary(inputVariantsFile, inputSourceFile, inputReferenceFile, replace);

        // get the original sequence dictionary from the source for comparison
        SAMSequenceDictionary sourceDictionary =
                SAMSequenceDictionaryExtractor.extractDictionary(
                        inputSourceFile == null ?
                                inputReferenceFile.toPath() : inputSourceFile.toPath()
                );

        // Some sequence dictionary sources will contain optional attributes (i.e., if the source is a .dict file,
        // or if only a reference is presented to the tool using -R, which will in turn cause the framework to retrieve
        // the dictionary from the .dict file accompanying the reference, the dictionary will likely include md5 and
        // UR attributes). However, htsjdk doesn't propagate these attributes to the VCF header properly
        // (https://github.com/samtools/htsjdk/issues/730), and many are stripped out. In order to ensure the
        // roundtrip comparison succeeds, roundtrip it through a VCF header to match what htsjdk will have written out.
        VCFHeader sourceVCFHeader = new VCFHeader();
        sourceVCFHeader.setSequenceDictionary(sourceDictionary);
        sourceDictionary = sourceVCFHeader.getSequenceDictionary();
        Assert.assertEquals(sourceDictionary, resultingDictionary);
    }

    @DataProvider(name="UpdateBadSequenceDictionaryData")
    public Object[][] updateBadSequenceDictionaryData() {
        return new Object[][]{
            // already has a dictionary
            new Object[]{ new File(testDir, "variantsWithDict.vcf"), new File(testDir, "variantsWithDict.vcf"), null, false},
            // source has no dictionary
            new Object[]{ new File(testDir, "variantsNoDict.vcf"), new File(testDir, "variantsNoDict.vcf"), null, false},
            // source dictionary is a mismatch for the variant records
            new Object[]{ new File(testDir, "variantsNoDictWithBadContig.vcf"), new File(testDir, "variantsWithDict.vcf"), null, false},
            new Object[]{ new File(testDir, "variantsNoDictWithBadContigLength.vcf"), new File(testDir, "variantsWithDict.vcf"), null, false},
        };
    }

    @Test(dataProvider="UpdateBadSequenceDictionaryData", expectedExceptions= CommandLineException.BadArgumentValue.class)
    private void testBadUpdateSequenceDictionary(
            final File inputVariantsFile,
            final File inputSourceFile,
            final File inputReferenceFile,
            final boolean replace) {
        updateSequenceDictionary(inputVariantsFile, inputSourceFile, inputReferenceFile, replace);
    }

    private SAMSequenceDictionary updateSequenceDictionary(
        final File inputVariantsFile,
        final File inputSourceFile,
        final File inputReferenceFile,
        final boolean replace)
    {
        ArgumentsBuilder argBuilder = new ArgumentsBuilder();

        argBuilder.addVCF(inputVariantsFile);
        if (inputReferenceFile != null) {
            argBuilder.addReference(inputReferenceFile);
        }
        if (inputSourceFile != null) {
            argBuilder.addFileArgument(UpdateVCFSequenceDictionary.DICTIONARY_ARGUMENT_NAME, inputSourceFile);
        }
        if (replace) {
            argBuilder.addArgument(UpdateVCFSequenceDictionary.REPLACE_ARGUMENT_NAME, Boolean.toString(replace));
        }

        File outFile = createTempFile("updateSequenceDictionary", ".vcf");
        argBuilder.addOutput(outFile);
        runCommandLine(argBuilder.getArgsList());

        // Don't require an index, since the framework doesn't create one if no input sequnce
        // dictionary is available via getBestAvailableSequenceDictionary.
        try (VCFFileReader vcfReader = new VCFFileReader(outFile, false)) {
            return vcfReader.getFileHeader().getSequenceDictionary();
        }
    }

}
