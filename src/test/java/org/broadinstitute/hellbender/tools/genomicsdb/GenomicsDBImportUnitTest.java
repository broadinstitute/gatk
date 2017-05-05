package org.broadinstitute.hellbender.tools.genomicsdb;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.mockito.internal.util.io.IOUtil;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.HashMap;
import java.util.Map;

public class GenomicsDBImportUnitTest extends BaseTest{

    @DataProvider
    public Object[][] getBadTestFiles(){
        return new Object[][]{
                {"Sample1\tsamplePath\n"
                +"Sample1\tsamplePath"},  //duplicate sample name
                {""}, // empty file
                {"Sample1\tSample2\tFile"}, //various wrong numbers of inputs
                {"Sample1\t"},
                {"Sample1"}
        };
    }

    @Test(dataProvider = "getBadTestFiles", expectedExceptions = UserException.BadInput.class)
    public void testBadInputFiles(final String text){
        final File sampleFile = IOUtils.writeTempFile(text, "badSampleMapping", ".txt");
        GenomicsDBImport.loadSampleNameMapFile(sampleFile.toPath()  );
    }

    @Test
    public void testSampleNameToFileMap(){
        final File sampleFile = IOUtils.writeTempFile("Sample1\tfile1\nSample2\tfile2", "sampleMapping", ".txt");
        final Map<String, String> expected = new HashMap<>();
        expected.put("Sample1", "file1");
        expected.put("Sample2", "file2");
        Assert.assertEquals(GenomicsDBImport.loadSampleNameMapFile(sampleFile.toPath()), expected);
    }
}