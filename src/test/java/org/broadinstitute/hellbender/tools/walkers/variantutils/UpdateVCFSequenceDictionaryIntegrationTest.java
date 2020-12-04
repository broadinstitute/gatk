package org.broadinstitute.hellbender.tools.walkers.variantutils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.net.URISyntaxException;
import java.nio.file.Paths;
import java.util.stream.Collectors;

public class UpdateVCFSequenceDictionaryIntegrationTest extends CommandLineProgramTest {
    private File testDir = new File(getTestDataDir(), "walkers/variantutils/UpdateVCFSequenceDictionary");

    @DataProvider(name="UpdateGoodSequenceDictionaryData")
    public Object[][] updateGoodSequenceDictionaryData() {
        return new Object[][]{
                new Object[]{ new File(testDir, "variantsNoDict.vcf"), new File(testDir, "variantsWithDict.vcf"), null, null, false, false},
                // pass a reference as a reference
                new Object[]{ new File(testDir, "variantsNoDict.vcf"), null, new File(testDir, "exampleFASTA.fasta"), null, false, false},
                // pass a reference as a source; we need to test both to ensure the user can bypass the framework sequence
                // dictionary validation that will occur if you use -R
                new Object[]{ new File(testDir, "variantsNoDict.vcf"), new File(testDir, "exampleFASTA.fasta"), null, null, false, false},
                new Object[]{ new File(testDir, "variantsNoDict.vcf"), new File(testDir, "exampleFASTA.dict"), null, null, false, false},
                // can't handle CRAM - see https://github.com/samtools/htsjdk/issues/731
                new Object[]{ new File(testDir, "variantsNoDict.vcf"), new File(testDir, "exampleBAM.bam"), null, null, false, false},
                // already has a dictionary - but force a replace; validation passes because vcf dictionary is subset of .dict
                new Object[]{ new File(testDir, "variantsWithSubsetDict.vcf"), new File(testDir, "exampleFASTA.dict"), null, null, true, false},
                // can force a replace with an invalid sequence dictionary if also disable sequence validation
                new Object[]{ new File(testDir, "variantsWithDictBadContigLength.vcf"), new File(testDir, "exampleFASTA.dict"), null, null, true, true}
        };
    }

    @Test(dataProvider="UpdateGoodSequenceDictionaryData")
    public void testGoodUpdateSequenceDictionary(
            final File inputVariantsFile,
            final File inputSourceFile,
            final File inputReferenceFile,
            final String masterSequenceDictionary,
            final boolean replace,
            final boolean disableSequenceDictionaryValidation) throws FileNotFoundException, URISyntaxException {
        final SAMSequenceDictionary resultingDictionary =
                updateSequenceDictionary(inputVariantsFile, inputSourceFile, inputReferenceFile, masterSequenceDictionary, replace, disableSequenceDictionaryValidation);

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
            new Object[]{ new File(testDir, "variantsWithDict.vcf"), new File(testDir, "variantsWithDict.vcf"), null, null, false, false},
            // source has no dictionary
            new Object[]{ new File(testDir, "variantsNoDict.vcf"), new File(testDir, "variantsNoDict.vcf"), null, null, false, false},
            // source dictionary is a mismatch for the variant records
            new Object[]{ new File(testDir, "variantsNoDictWithBadContig.vcf"), new File(testDir, "variantsWithDict.vcf"), null, null, false, false},
            new Object[]{ new File(testDir, "variantsNoDictWithBadContigLength.vcf"), new File(testDir, "variantsWithDict.vcf"), null, null, false, false},
        };
    }

    @Test(dataProvider="UpdateBadSequenceDictionaryData", expectedExceptions= CommandLineException.BadArgumentValue.class)
    public void testBadUpdateSequenceDictionary(
            final File inputVariantsFile,
            final File inputSourceFile,
            final File inputReferenceFile,
            final String masterSequenceDictionary,
            final boolean replace,
            final boolean disableSequenceDictionaryValidation) {
        updateSequenceDictionary(inputVariantsFile, inputSourceFile, inputReferenceFile, masterSequenceDictionary, replace, disableSequenceDictionaryValidation);
    }

    @Test
    public void testUseMasterDictionary() {
        final SAMSequenceDictionary actualSequenceDictionary = updateSequenceDictionary(
                new File(testDir, "variantsNoDict.vcf"),
                null,
                null,
                new File(testDir, "exampleFASTA.dict").getAbsolutePath(),
                false,
                false);
        final SAMSequenceDictionary expectedSequenceDictionary = SAMSequenceDictionaryExtractor.extractDictionary(
                Paths.get(new File(testDir, "exampleFASTA.dict").getAbsolutePath()));

        // verify only the sequence names and lengths, since other attributes such as MD/UR will have been updated
        Assert.assertEquals(actualSequenceDictionary.getSequences().stream().map(seq -> seq.getSequenceName()).collect(Collectors.toList()),
                expectedSequenceDictionary.getSequences().stream().map(seq -> seq.getSequenceName()).collect(Collectors.toList()));
        Assert.assertEquals(actualSequenceDictionary.getSequences().stream().map(seq -> seq.getSequenceLength()).collect(Collectors.toList()),
                expectedSequenceDictionary.getSequences().stream().map(seq -> seq.getSequenceLength()).collect(Collectors.toList()));
    }


    @Test(expectedExceptions = UserException.SequenceDictionaryIsMissingContigLengths.class)
    public void testBadContigLengthWithValidation() {
        // throw an error if trying to force a replace with an invalid sequence dictionary without disabling sequence validation
        updateSequenceDictionary(
                new File(testDir, "variantsWithDict.vcf"),
                new File(testDir, "variantsWithDictBadContigLength.vcf"),
                null,
                null,
                true,
                false);
    }

    @Test(expectedExceptions=CommandLineException.class)
    public void testMasterDictionaryAmbiguous() {
        // specifying both a source dictionary and a master dictionary is ambiguous
        updateSequenceDictionary(
                new File(testDir, "variantsNoDict.vcf"),
                new File(testDir, "exampleFASTA.dict"),
                null,
                new File(testDir, "exampleFASTA.dict").getAbsolutePath(),
                false,
                false);
    }

    private SAMSequenceDictionary updateSequenceDictionary(
        final File inputVariantsFile,
        final File inputSourceFile,
        final File inputReferenceFile,
        final String masterSequenceDictionary,
        final boolean replace,
        final boolean disableSequenceDictionaryValidation)
    {
        ArgumentsBuilder argBuilder = new ArgumentsBuilder();

        argBuilder.addVCF(inputVariantsFile);
        if (inputReferenceFile != null) {
            argBuilder.addReference(inputReferenceFile);
        }
        if (inputSourceFile != null) {
            argBuilder.add(UpdateVCFSequenceDictionary.DICTIONARY_ARGUMENT_NAME, inputSourceFile);
        }
        if (replace) {
            argBuilder.add(UpdateVCFSequenceDictionary.REPLACE_ARGUMENT_NAME, Boolean.toString(replace));
        }
        if (disableSequenceDictionaryValidation) {
            argBuilder.add(StandardArgumentDefinitions.DISABLE_SEQUENCE_DICT_VALIDATION_NAME, Boolean.toString(disableSequenceDictionaryValidation));
        }
        if (masterSequenceDictionary != null) {
            argBuilder.add(StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME, masterSequenceDictionary);
        }

        File outFile = createTempFile("updateSequenceDictionary", ".vcf.gz");
        argBuilder.addOutput(outFile);
        runCommandLine(argBuilder.getArgsList());

        // Don't require an index, since the framework doesn't create one if no input sequence
        // dictionary is available via getBestAvailableSequenceDictionary.
        try (VCFFileReader vcfReader = new VCFFileReader(outFile, false)) {
            return vcfReader.getFileHeader().getSequenceDictionary();
        }
    }

}
