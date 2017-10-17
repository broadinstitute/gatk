package org.broadinstitute.hellbender.tools.walkers.readorientation;

import org.broadinstitute.hellbender.utils.Nucleotide;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;

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
        for (int i = 0; i < numRowsPerContext; i++){
            for (String referenceContext : referenceContexts){
                final int[] counts = new int[] { 0, 0, 0, i };
                final int[] f1r2Coutns = new int[] { 0, 0, 0, i/2 };
                AltSiteRecord record = new AltSiteRecord(referenceContext, counts, f1r2Coutns, depth, Nucleotide.T);
                writer.writeRecord(record);
            }
        }

        writer.close();

        List<AltSiteRecord> altDesignMatrix = AltSiteRecord.readAltSiteRecords(table, numRowsPerContext*3);
        Assert.assertEquals(altDesignMatrix.size(), numRowsPerContext * 3);
        Assert.assertEquals(altDesignMatrix.stream()
                .filter(rec -> rec.getReferenceContext().equals("ACT")).count(), numRowsPerContext);
        Assert.assertEquals(altDesignMatrix.stream()
                .filter(rec -> rec.getReferenceContext().equals("GAT")).count(), numRowsPerContext);
        Assert.assertEquals(altDesignMatrix.stream()
                .filter(rec -> rec.getReferenceContext().equals("CCA")).count(), numRowsPerContext);
    }

}