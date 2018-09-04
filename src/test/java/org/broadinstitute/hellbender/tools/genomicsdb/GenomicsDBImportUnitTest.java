package org.broadinstitute.hellbender.tools.genomicsdb;

import com.intel.genomicsdb.importer.GenomicsDBImporter;
import htsjdk.tribble.FeatureReader;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.LinkedHashMap;
import java.util.Map;

public class GenomicsDBImportUnitTest extends GATKBaseTest {

    private static final String ORDERED_SAMPLE_MAP =    "Sample1\tfile1\n" +
                                                        "Sample2\tfile2\n" +
                                                        "Sample3\tfile3";

    private static final String UNORDERED_SAMPLE_MAP =  "Sample3\tfile3\n" +
                                                        "Sample2\tfile2\n" +
                                                        "Sample1\tfile1\n";

    @DataProvider
    public Object[][] getBadSampleNameMapFiles(){
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
                {"Sample1 file1"},              // 1 column
                {" name name\tfile1"},          // preceding whitespace
                {"name name \tfile1"},          // trailing whitespace
        };
    }

    @Test(dataProvider = "getBadSampleNameMapFiles", expectedExceptions = UserException.BadInput.class)
    public void testBadInputFiles(final String text){
        final File sampleFile = IOUtils.writeTempFile(text, "badSampleMapping", ".txt");
        GenomicsDBImport.loadSampleNameMapFile(sampleFile.toPath()  );
    }

    @DataProvider
    public Object[][] getGoodSampleNameMapFileSyntax(){
        return new Object[][]{
                // Note: none of these files are real, these are just valid files syntactically
                {"Sample1\tsamplePath1 \n"
                +"Sample2\tsamplePath2", new String[][] {{"Sample1","samplePath1"},{"Sample2","samplePath2"}}},     // normal sample names
                {"Sample1 001\tFile", new String[][] {{"Sample1 001","File"}}},          // sample names with whitespace
                {"name name\tfile1 ", new String[][] {{"name name","file1"}}},          // trailing whitespace second column
                {"name name\t file1 ", new String[][] {{"name name","file1"}}}        // leading and trailing whitespace second colum
                };
    }

    @Test(dataProvider = "getGoodSampleNameMapFileSyntax")
    public void testValidSampleFiles(final String text, final String[][] expectedEntries){
        final File sampleFile = IOUtils.writeTempFile(text, "goodSampleMapping", ".txt");
        final LinkedHashMap<String, Path> outputMap = GenomicsDBImport.loadSampleNameMapFile(sampleFile.toPath());
        Assert.assertEquals(outputMap.size(),expectedEntries.length);

        Arrays.stream(expectedEntries).forEach(s -> { Assert.assertTrue(outputMap.containsKey(s[0]));
                                                      Assert.assertEquals(outputMap.get(s[0]).toString(),s[1]);});
    }

    @Test
    public void testLoadSampleNameMapFilePreservesOrder(){
        final File sampleFile = IOUtils.writeTempFile(UNORDERED_SAMPLE_MAP, "badSampleMapping", ".txt");
        final LinkedHashMap<String, Path> unsortedMap = GenomicsDBImport.loadSampleNameMapFile(sampleFile.toPath());
        Assert.assertEquals(new ArrayList<>(unsortedMap.keySet()), Arrays.asList("Sample3", "Sample2", "Sample1"));
    }

    @DataProvider
    public Object[][] getSampleMaps(){
        return new Object[][]{
                {ORDERED_SAMPLE_MAP},
                {UNORDERED_SAMPLE_MAP}
        };
    }

    @Test(dataProvider = "getSampleMaps")
    public void testLoadSampleNameMapFileInSortedOrder(final String sampleMapText){
        final File sampleFile = IOUtils.writeTempFile(sampleMapText, "sampleMapping", ".txt");
        final Map<String, Path> expected = new LinkedHashMap<>();
        expected.put("Sample1", Paths.get("file1"));
        expected.put("Sample2", Paths.get("file2"));
        expected.put("Sample3", Paths.get("file3"));
        final Map<String, Path> actual = GenomicsDBImport.loadSampleNameMapFileInSortedOrder(sampleFile.toPath());
        Assert.assertEquals(actual, expected);
        Assert.assertEquals(actual.keySet().iterator().next(), "Sample1");
    }
}
