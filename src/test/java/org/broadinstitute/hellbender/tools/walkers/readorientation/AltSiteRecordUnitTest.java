package org.broadinstitute.hellbender.tools.walkers.readorientation;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.testutils.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * Created by tsato on 10/16/17.
 */
public class AltSiteRecordUnitTest extends GATKBaseTest {
    @Test
    public void test() throws IOException {
        File table = BaseTest.createTempFile("alt-site-design-matrix","tsv");
        AltSiteRecord.AltSiteRecordTableWriter writer = new AltSiteRecord.AltSiteRecordTableWriter(table);

        List<String> referenceContexts = Arrays.asList("ACT", "GAT", "CCA");
        final int numRowsPerContext = 100;
        final int refCount = 10;
        final int altCount = 10;
        final int depth = refCount + altCount;
        final int refF1R2 = 5;
        final int altF1R2 = 5;
        final Nucleotide altAllele = Nucleotide.T;

        for (int i = 0; i < numRowsPerContext; i++){
            for (String referenceContext : referenceContexts){
                AltSiteRecord record = new AltSiteRecord(referenceContext, refCount, altCount, refF1R2, altF1R2, altAllele);
                writer.writeRecord(record);
            }
        }

        writer.close();

        List<AltSiteRecord> altDesignMatrix = AltSiteRecord.readAltSiteRecords(table, numRowsPerContext * referenceContexts.size());
        final Map<String, List<AltSiteRecord>> recordsByContext = altDesignMatrix.stream()
                .collect(Collectors.groupingBy(AltSiteRecord::getReferenceContext));
        recordsByContext.values().forEach(recs -> Assert.assertEquals(recs.size(), numRowsPerContext));
        referenceContexts.forEach(context -> Assert.assertTrue(recordsByContext.containsKey(context)));

        for (String referenceContext : referenceContexts){
            List<AltSiteRecord> recs = recordsByContext.get(referenceContext);
            final int sumBaseCounts = recs.stream().mapToInt(rec -> rec.getDepth()).sum();
            Assert.assertEquals(sumBaseCounts, depth * numRowsPerContext);
            final int sumF1R2Counts = recs.stream().mapToInt(rec -> rec.getAltF1R2() + rec.getRefF1R2()).sum();
            Assert.assertEquals(sumF1R2Counts, (refF1R2 + altF1R2)*numRowsPerContext);
        }
    }

    @Test
    public void testGetReverseComplement() throws IOException {
        final String referenceContext = "AGC"; // G -> A
        final String revCompContext = "GCT"; // C -> T
        final Nucleotide revCompAlt = Nucleotide.T;
        final Nucleotide standardAlt = Nucleotide.A;
        final int refCount = 60;
        final int altCount = 10;
        final int depth = refCount + altCount;
        final int refF1R2 = 15;
        final int altF1R2 = 2;

        // Create a reverse comp version of an alt site record and change it into the standard form
        final AltSiteRecord revComp = new AltSiteRecord(revCompContext, refCount, altCount, refF1R2, altF1R2, revCompAlt);
        final AltSiteRecord standard = revComp.getReverseComplementOfRecord();

        Assert.assertEquals(standard.getReferenceContext(), referenceContext);
        Assert.assertEquals(standard.getAltAllele(), standardAlt);
        Assert.assertEquals(standard.getRefCount(), refCount);
        Assert.assertEquals(standard.getAltCount(), altCount);
        Assert.assertEquals(standard.getAltF1R2(), altCount - altF1R2);
        Assert.assertEquals(standard.getRefF1R2(), refCount - refF1R2);
        Assert.assertEquals(standard.getDepth(), depth);

    }

}