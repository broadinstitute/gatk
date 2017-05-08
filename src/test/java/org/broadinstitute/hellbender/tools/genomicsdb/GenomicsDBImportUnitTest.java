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
import java.util.LinkedHashMap;
import java.util.Map;

public class GenomicsDBImportUnitTest extends BaseTest{

    @DataProvider
    public Object[][] getBadTestFiles(){
        return new Object[][]{
                {"Sample1\tsamplePath\n"
                +"Sample1\tsamplePath"},        // duplicate sample name
                {""},                           // empty file
                {"Sample1\tSample2\tFile"},     // 3 columns
                {"Sample1\t"},                  // 1 column
                {"Sample1"},                    // 1 column no delimiter
                {"\tfile"},                     // empty first token
                {" \tfile"},                    // first token only whitespace
                {"Sample1\tfile1\t"},           // extra tab
                {"Sample1\nfile"},              // newline instead of tab
                {"\t"},                         // only tab
                {"Sample1 file1"},              // non-tab whitespace
                {"Sample 1\tfile1"},            // white space in a token
        };
    }

    @Test(dataProvider = "getBadTestFiles", expectedExceptions = UserException.BadInput.class)
    public void testBadInputFiles(final String text){
        final File sampleFile = IOUtils.writeTempFile(text, "badSampleMapping", ".txt");
        GenomicsDBImport.loadSampleNameMapFile(sampleFile.toPath()  );
    }

    @Test
    public void testSampleNameToFileMap(){
        final String content =
                "Sample1\tfile1\n" +
                "Sample2\tfile2\n" +
                "Sample3\tfile3";
        final File sampleFile = IOUtils.writeTempFile(content, "sampleMapping", ".txt");
        final Map<String, String> expected = new LinkedHashMap<>();
        expected.put("Sample1", "file1");
        expected.put("Sample2", "file2");
        expected.put("Sample3", "file3");
        final LinkedHashMap<String, String> actual = GenomicsDBImport.loadSampleNameMapFile(sampleFile.toPath());
        Assert.assertEquals(actual, expected);
        Assert.assertEquals(actual.keySet().iterator().next(), "Sample1");
    }
}