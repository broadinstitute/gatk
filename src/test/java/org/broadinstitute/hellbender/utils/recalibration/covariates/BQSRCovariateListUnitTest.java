package org.broadinstitute.hellbender.utils.recalibration.covariates;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationArgumentCollection;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Collections;
import java.util.List;

public final class BQSRCovariateListUnitTest extends GATKBaseTest {

    public BQSRCovariateList makeCovariateList() {
        return new BQSRCovariateList(new RecalibrationArgumentCollection(), Collections.singletonList("readGroup"));
    }

    @Test
    public void testSize() {
        BQSRCovariateList scl = makeCovariateList();
        Assert.assertEquals(scl.size(), 4);
    }

    @Test
    public void testCovariateNames() {
        BQSRCovariateList scl = makeCovariateList();
        Assert.assertEquals(scl.covariateNames(), "ReadGroupCovariate,QualityScoreCovariate,ContextCovariate,CycleCovariate");
    }

    @Test
    public void testIterator() {
        BQSRCovariateList scl = makeCovariateList();
        Assert.assertEquals(Utils.stream(scl).count(), 4);
    }

    @Test
    public void testGetCovariates() {
        BQSRCovariateList scl = makeCovariateList();
        Assert.assertEquals(scl.getReadGroupCovariate().parseNameForReport(), "ReadGroup");
        Assert.assertEquals(scl.getQualityScoreCovariate().parseNameForReport(), "QualityScore");
        Assert.assertEquals(scl.getAdditionalCovariates().get(0).parseNameForReport(), "Context");
        Assert.assertEquals(scl.getAdditionalCovariates().get(1).parseNameForReport(), "Cycle");
    }

    @Test
    public void testGetCovariatesByIndex() {
        BQSRCovariateList scl = makeCovariateList();
        Assert.assertEquals(scl.get(0).parseNameForReport(), "ReadGroup");
        Assert.assertEquals(scl.get(1).parseNameForReport(), "QualityScore");
        Assert.assertEquals(scl.get(2).parseNameForReport(), "Context");
        Assert.assertEquals(scl.get(3).parseNameForReport(), "Cycle");
    }

    @Test(expectedExceptions = IndexOutOfBoundsException.class)
    public void testGetCovariatesByIndexInvalid() {
        BQSRCovariateList scl = makeCovariateList();
        scl.get(4);
    }

    @Test
    public void testGetCovariatesByIndexClass() {
        BQSRCovariateList scl = makeCovariateList();
        Assert.assertEquals(scl.indexByClass(ReadGroupCovariate.class), 0);
        Assert.assertEquals(scl.indexByClass(QualityScoreCovariate.class), 1);
        Assert.assertEquals(scl.indexByClass(ContextCovariate.class), 2);
        Assert.assertEquals(scl.indexByClass(CycleCovariate.class), 3);

        //finally, test an anonymous subclass
        Assert.assertEquals(scl.indexByClass(new Covariate() {
            private static final long serialVersionUID = 1L;

            @Override
            public void initialize( RecalibrationArgumentCollection RAC, List<String> readGroups ) {
                
            }

            @Override
            public void recordValues(GATKRead read, SAMFileHeader header, PerReadCovariateMatrix values, boolean recordIndels) {

            }

            @Override
            public String formatKey(int key) {
                return null;
            }

            @Override
            public int keyFromValue(Object value) {
                return 0;
            }

            @Override
            public int maximumKeyValue() {
                return 0;
            }
        }.getClass()), -1);
    }

    @Test
    public void testGetCovariatesByParsedName() {
        BQSRCovariateList scl = makeCovariateList();
        final String[] parsedNames = {"ReadGroup", "QualityScore", "Context", "Cycle"};
        for (String parsedName : parsedNames) {
            Assert.assertEquals(scl.getCovariateByParsedName(parsedName).parseNameForReport(), parsedName);
        }
        Assert.assertEquals(scl.getCovariateByParsedName("fred"), null);
    }

    @Test
    public void testCovariateInitialize() {
        BQSRCovariateList scl = makeCovariateList();
        //this just tests non blowing up.
    }
}
