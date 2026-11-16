package org.broadinstitute.hellbender.tools.walkers.readorientation;

import java.nio.file.Path;
import org.apache.commons.lang3.tuple.Pair;
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
    // a table of several reference contexts, as F1R2CountsCollector writes before splitting it by context, has a column of contexts
    @Test
    public void test() throws IOException {
        Path table = BaseTest.createTempPath("alt-site-design-matrix","tsv");
        AltSiteRecord.AltSiteRecordTableWriter writer = new AltSiteRecord.AltSiteRecordTableWriter(table, "SAMPLE", true);

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
                AltSiteRecord record = new AltSiteRecord(referenceContext, refCount, altCount, refF1R2, altF1R2, altAllele, "1", i + 1);
                writer.writeRecord(record);
            }
        }

        writer.close();

        final List<AltSiteRecord> altDesignMatrix;
        try (AltSiteRecord.AltSiteRecordTableReader reader = new AltSiteRecord.AltSiteRecordTableReader(table)) {
            altDesignMatrix = reader.toList();
        }
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

    // a table of one reference context, as LearnReadOrientationModel reads, has the context in its name rather than a column
    @Test
    public void testTableOfOneContext() throws IOException {
        final Path table = BaseTest.createTempPath("alt-site-design-matrix", "tsv");
        final String referenceContext = "ACT";
        final int numRows = 10;
        try (AltSiteRecord.AltSiteRecordTableWriter writer = new AltSiteRecord.AltSiteRecordTableWriter(table, "SAMPLE")) {
            for (int i = 0; i < numRows; i++) {
                writer.writeRecord(new AltSiteRecord(referenceContext, 10 + i, 5, 4, 3, Nucleotide.G, "chr2", 1000 * (i + 1)));
            }
        }

        final Pair<String, List<AltSiteRecord>> sampleAndRecords = AltSiteRecord.readAltSiteRecords(table, referenceContext);
        Assert.assertEquals(sampleAndRecords.getLeft(), "SAMPLE");
        final List<AltSiteRecord> records = sampleAndRecords.getRight();
        Assert.assertEquals(records.size(), numRows);
        for (int i = 0; i < numRows; i++) {
            final AltSiteRecord record = records.get(i);
            Assert.assertEquals(record.getReferenceContext(), referenceContext);
            Assert.assertEquals(record.getRefCount(), 10 + i);
            Assert.assertEquals(record.getAltCount(), 5);
            Assert.assertEquals(record.getDepth(), 15 + i);
            Assert.assertEquals(record.getRefF1R2(), 4);
            Assert.assertEquals(record.getAltF1R2(), 3);
            Assert.assertEquals(record.getAltAllele(), Nucleotide.G);
            Assert.assertEquals(record.getContig(), "chr2");
            Assert.assertEquals(record.getPosition(), 1000 * (i + 1));
        }

        // the name of an alt table of the F1R2 counts gives its context
        Assert.assertEquals(F1R2CountsCollector.getAltTableContext(new File(F1R2CountsCollector.getAltTableName("tumor.1", "ACT"))), "ACT");
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
        final AltSiteRecord revComp = new AltSiteRecord(revCompContext, refCount, altCount, refF1R2, altF1R2, revCompAlt, "chr3", 12345);
        final AltSiteRecord standard = revComp.getReverseComplementOfRecord();

        Assert.assertEquals(standard.getReferenceContext(), referenceContext);
        Assert.assertEquals(standard.getAltAllele(), standardAlt);
        Assert.assertEquals(standard.getRefCount(), refCount);
        Assert.assertEquals(standard.getAltCount(), altCount);
        Assert.assertEquals(standard.getAltF1R2(), altCount - altF1R2);
        Assert.assertEquals(standard.getRefF1R2(), refCount - refF1R2);
        Assert.assertEquals(standard.getDepth(), depth);
        Assert.assertEquals(standard.getContig(), "chr3");
        Assert.assertEquals(standard.getPosition(), 12345);

    }

}
