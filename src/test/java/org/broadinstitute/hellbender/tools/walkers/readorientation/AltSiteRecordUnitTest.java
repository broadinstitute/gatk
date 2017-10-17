package org.broadinstitute.hellbender.tools.walkers.readorientation;

import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.testng.Assert;
import org.testng.annotations.Test;
import org.testng.internal.junit.ArrayAsserts;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import static org.testng.Assert.*;

/**
 * Created by tsato on 10/16/17.
 */
public class AltSiteRecordUnitTest {
    @Test
    public void test() throws IOException {
        File table = File.createTempFile("alt-site-design-matrix","tsv");
        AltSiteRecord.AltSiteRecordTableWriter writer = new AltSiteRecord.AltSiteRecordTableWriter(table);

        List<String> referenceContexts = Arrays.asList("ACT", "GAT", "CCA");
        final int numRowsPerContext = 100;
        final int depth = 200;
        final String contig = "1";
        final int position = 100_000_000;

        for (int i = 0; i < numRowsPerContext; i++){
            for (String referenceContext : referenceContexts){
                final int[] counts = new int[] { 0, 0, 0, i };
                final int[] f1r2Coutns = new int[] { 0, 0, 0, i };
                AltSiteRecord record = new AltSiteRecord(contig, position, referenceContext, counts,
                        f1r2Coutns, depth, Nucleotide.T);
                writer.writeRecord(record);
            }
        }

        writer.close();

        List<AltSiteRecord> altDesignMatrix = AltSiteRecord.readAltSiteRecords(table, numRowsPerContext * referenceContexts.size());
        final Map<String, List<AltSiteRecord>> recordsByContext = altDesignMatrix.stream()
                .collect(Collectors.groupingBy(AltSiteRecord::getReferenceContext));
        recordsByContext.values().forEach(recs -> Assert.assertEquals(recs.size(), numRowsPerContext));
        referenceContexts.forEach(cont -> Assert.assertTrue(recordsByContext.containsKey(cont)));

        for (String referenceContext : referenceContexts){
            List<AltSiteRecord> recs = recordsByContext.get(referenceContext);
            final long sumBaseCounts = recs.stream().mapToLong(rec -> MathUtils.sum(rec.getBaseCounts())).sum();
            Assert.assertEquals(sumBaseCounts, (numRowsPerContext-1)*numRowsPerContext/2); // use  1 + 2 + , ..., + n = n(n+1)/2
            final long sumF1R2Counts = recs.stream().mapToLong(rec -> MathUtils.sum(rec.getF1R2Counts())).sum();
            Assert.assertEquals(sumF1R2Counts, (numRowsPerContext-1)*numRowsPerContext/2);
        }
    }

    @Test
    public void testReverseComplement() throws IOException {
        final String contig = "1";
        final int position = 100_000;
        final int depth = 200;
        final String referenceContext = "AGC";
        final Nucleotide altAllele = Nucleotide.A; // G -> A transition
        final String revCompContext = "GCT";
        final int[] baseCounts = new int[]{ 10, 0, 60, 0 };
        final int[] f1r2Counts = new int[]{ 2, 0, 15, 0 };

        final AltSiteRecord record = new AltSiteRecord(contig, position, referenceContext, baseCounts, f1r2Counts, depth, altAllele);
        final AltSiteRecord revComp = record.getReverseComplementOfRecord();

        Assert.assertEquals(revComp.getReferenceContext(), revCompContext);
        ArrayAsserts.assertArrayEquals(revComp.getBaseCounts(), new int[]{ 0, 60, 0, 10 } );
        ArrayAsserts.assertArrayEquals(revComp.getF1R2Counts(), new int[]{ 0, 15, 0, 2 } );

    }

}