package org.broadinstitute.hellbender.tools.genomicsdb;

import com.intel.genomicsdb.GenomicsDBImporter;
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
                {"Sample1 file1"}               // non-tab whitespace
        };
    }

    @Test(dataProvider = "getBadTestFiles", expectedExceptions = UserException.BadInput.class)
    public void testBadInputFiles(final String text){
        final File sampleFile = IOUtils.writeTempFile(text, "badSampleMapping", ".txt");
        GenomicsDBImport.loadSampleNameMapFile(sampleFile.toPath()  );
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

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testNullFeatureReadersToFail() {
        final Map<String, FeatureReader<VariantContext>> sampleToReaderMap = new LinkedHashMap<>();
        sampleToReaderMap.put("Sample1", null);
        GenomicsDBImporter.generateSortedCallSetMap(sampleToReaderMap, true, true, 0L);
    }
}
