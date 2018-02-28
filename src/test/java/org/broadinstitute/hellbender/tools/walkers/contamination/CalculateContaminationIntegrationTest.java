package org.broadinstitute.hellbender.tools.walkers.contamination;

import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.distribution.UniformRealDistribution;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.RandomGeneratorFactory;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.stream.IntStream;

/**
 * Created by David Benjamin on 2/16/17.
 */
public class CalculateContaminationIntegrationTest extends CommandLineProgramTest {
    public static final File SPIKEIN_DATA_DIRECTORY = new File(getTestDataDir(), "calculatecontamination");
    public static final File NA12891_1_PCT_NA12892_99_PCT = new File(SPIKEIN_DATA_DIRECTORY, "NA12891_0.01_NA12892_0.99.table");
    public static final File NA12891_3_PCT_NA12892_97_PCT = new File(SPIKEIN_DATA_DIRECTORY, "NA12891_0.03_NA12892_0.97.table");
    public static final File NA12891_5_PCT_NA12892_95_PCT = new File(SPIKEIN_DATA_DIRECTORY, "NA12891_0.05_NA12892_0.95.table");
    public static final File NA12891_8_PCT_NA12892_92_PCT = new File(SPIKEIN_DATA_DIRECTORY, "NA12891_0.08_NA12892_0.92.table");
    public static final double BASELINE_CONTAMINATION_OF_NA12892 = 0.001;

    @Test
    public void testArtificialData() {
        final RandomGenerator randomGenerator = RandomGeneratorFactory.createRandomGenerator(new Random(111));

        final String contig = "chr1";
        final int spacing = 100000;
        final int numSites = 2000;
        final int lohStart = 1000;
        final int lohEnd = 1100;
        final double lohMinorAlleleFraction = 0.2;
        final double[] minorAlleleFractions = Collections.nCopies(numSites, 0.5).stream().mapToDouble(x -> x).toArray();
        IntStream.range(lohStart, lohEnd).forEach(n -> minorAlleleFractions[n] = lohMinorAlleleFraction);


        final int sampleDepth = 100;
        final int contaminantDepth = 7;
        final int totalDepth = sampleDepth + contaminantDepth;
        final double contamination = (double) contaminantDepth / (sampleDepth + contaminantDepth);

        final UniformRealDistribution alleleFrequencyDistribution = new UniformRealDistribution(randomGenerator, 0.05, 0.3);
        final List<PileupSummary> ps = new ArrayList<>();

        for (int n = 0; n < numSites; n++) {
            final int position = n * spacing;
            final double alleleFrequency = alleleFrequencyDistribution.sample();

            // get the contaminant "genotype"
            final int contaminantAltCount = new BinomialDistribution(randomGenerator, contaminantDepth, alleleFrequency).sample();
            int totalAltCount = contaminantAltCount;

            final double hetProbability = 2 * alleleFrequency * (1 - alleleFrequency);
            final double homAltProbability = alleleFrequency * alleleFrequency;
            final double x = randomGenerator.nextDouble();
            if (x < hetProbability) {    // het
                //draw alt allele fractions from mixture of binomials centered at maf and 1 - maf
                final double altFraction = randomGenerator.nextDouble() < 0.5 ? minorAlleleFractions[n] : 1 - minorAlleleFractions[n];
                totalAltCount += new BinomialDistribution(randomGenerator, sampleDepth, altFraction).sample();

            } else if (x < hetProbability + homAltProbability) {    //hom alt
                totalAltCount += sampleDepth;
            } else {    // hom ref
                //do nothing -- no alts from the sample
            }
            final int totalRefCount = totalDepth - totalAltCount;
            ps.add(new PileupSummary(contig, position, totalRefCount, totalAltCount, 0, alleleFrequency));
        }

        final File psTable = createTempFile("pileups", ".table");
        PileupSummary.writeToFile(ps, psTable);
        final File contaminationTable = createTempFile("contamination", ".table");
        final File segmentationsTable = createTempFile("segments", ".table");

        final String[] args = {
                "-I", psTable.getAbsolutePath(),
                "-O", contaminationTable.getAbsolutePath(),
                "-" + CalculateContamination.TUMOR_SEGMENTATION_SHORT_NAME, segmentationsTable.getAbsolutePath()
        };
        runCommandLine(args);

        final double calculatedContamination = ContaminationRecord.readFromFile(contaminationTable).get(0).getContamination();
        final List<MinorAlleleFractionRecord> mafRecords = MinorAlleleFractionRecord.readFromFile(segmentationsTable);

        Assert.assertEquals(calculatedContamination, contamination, 0.01);
        Assert.assertEquals(mafRecords.size(), 3);
        Assert.assertEquals(mafRecords.get(0).getMinorAlleleFraction(), 0.5, 0.05);
        Assert.assertEquals(mafRecords.get(1).getMinorAlleleFraction(), lohMinorAlleleFraction, 0.05);
        Assert.assertEquals(mafRecords.get(2).getMinorAlleleFraction(), 0.5, 0.05);
    }


    // spike-in where even at 0% spikein the bam file still had some baseline contamination
    @Test(dataProvider = "spikeInData")
    public void testSpikeIn(final File pileupSummary, final double spikeIn, final double baselineContamination) {
        final File contaminationTable = createTempFile("contamination", ".table");
        final double contamination = spikeIn + baselineContamination;

        final String[] args = {
                "-I", pileupSummary.getAbsolutePath(),
                "-O", contaminationTable.getAbsolutePath(),
        };
        runCommandLine(args);

        final double calculatedContamination = ContaminationRecord.readFromFile(contaminationTable).get(0).getContamination();
        Assert.assertEquals(calculatedContamination, contamination, 0.015);
    }

    // pileup summary table, spikein fraction, baseline contamination before spike-in
    @DataProvider(name = "spikeInData")
    public Object[][] spikeInData() {
        return new Object[][]{
                {NA12891_1_PCT_NA12892_99_PCT, 0.01, BASELINE_CONTAMINATION_OF_NA12892},
                {NA12891_3_PCT_NA12892_97_PCT, 0.03, BASELINE_CONTAMINATION_OF_NA12892},
                {NA12891_5_PCT_NA12892_95_PCT, 0.05, BASELINE_CONTAMINATION_OF_NA12892},
                {NA12891_8_PCT_NA12892_92_PCT, 0.08, BASELINE_CONTAMINATION_OF_NA12892}
        };
    }

    // Using 1% spike-in as a very close approximation to an uncontaminated matched normal.
    @Test
    public void testMatchedNormal() {
        final File normal = NA12891_1_PCT_NA12892_99_PCT;
        final File contaminated = NA12891_8_PCT_NA12892_92_PCT;
        final double baselineContamination = BASELINE_CONTAMINATION_OF_NA12892;
        final double contamination = 0.08 + baselineContamination;
        final File contaminationTable = createTempFile("contamination", ".table");


        final String[] args = {
                "-I", contaminated.getAbsolutePath(),
                "-" + CalculateContamination.MATCHED_NORMAL_SHORT_NAME, normal.getAbsolutePath(),
                "-O", contaminationTable.getAbsolutePath(),
        };
        runCommandLine(args);

        final double calculatedContamination = ContaminationRecord.readFromFile(contaminationTable).get(0).getContamination();
        Assert.assertEquals(calculatedContamination, contamination, 0.01);
    }
}