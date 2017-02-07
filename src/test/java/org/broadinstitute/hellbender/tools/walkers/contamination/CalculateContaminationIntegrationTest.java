package org.broadinstitute.hellbender.tools.walkers.contamination;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by David Benjamin on 2/16/17.
 */
public class CalculateContaminationIntegrationTest extends CommandLineProgramTest {
    @Test
    public void test() {
        final String contig = "chr1";
        final int spacing = 100000;
        final double contamination = 0.07;
        final double alleleFrequency = 0.2;
        final double refContamination = contamination * (1 - alleleFrequency);

        final int depth = 100;

        final List<PileupSummary> ps = new ArrayList<>();

        // if allele frequency is 0.2, then hets are about 10 times as common as hom var
        // thus we make nine hets for every hom var
        for (int n = 0; n < 1000; n++) {
            int position = n * spacing;
                ps.add(n % 10 == 0
                        ? new PileupSummary(contig, position, (int) (depth*refContamination), (int) (depth * (1 - refContamination)), 0, alleleFrequency)
                        : new PileupSummary(contig, position, depth / 2, depth/2, 1, alleleFrequency));
        }

        // same idea, but loss of heterozygosity
        // as above we make nine hets for every hom var, but these are LoH hets with a skewed allele fration of 0.8
        for (int n = 1000; n < 1100; n++) {
            int position = n * spacing;
            ps.add(n % 10 == 0
                    ? new PileupSummary(contig, position, (int) (depth*refContamination), (int) (depth * (1 - refContamination)), 0, alleleFrequency)
                    : new PileupSummary(contig, position, (int) 0.2 * depth, (int) (0.8 * depth), 1, alleleFrequency));
        }

        final File psTable = createTempFile("pileups", ".table");
        PileupSummary.writePileupSummaries(ps, psTable);
        final File contaminationTable = createTempFile("contamination", ".table");

        final String[] args = {
                "-I", psTable.getAbsolutePath(),
                "-O", contaminationTable.getAbsolutePath(),
        };
        runCommandLine(args);

        final double calculatedContamination = ContaminationRecord.readContaminationTable(contaminationTable).get(0).getContamination();
        Assert.assertEquals(calculatedContamination, contamination, 0.01);
    }

}