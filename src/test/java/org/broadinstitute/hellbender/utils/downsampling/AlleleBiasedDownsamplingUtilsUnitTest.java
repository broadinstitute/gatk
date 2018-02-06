package org.broadinstitute.hellbender.utils.downsampling;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.variant.variantcontext.Allele;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.LogManager;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.*;


/**
 * Basic unit test for AlleleBiasedDownsamplingUtils
 */
public class AlleleBiasedDownsamplingUtilsUnitTest extends GATKBaseTest {

    private static final org.apache.logging.log4j.Logger logger = LogManager.getLogger(AlleleBiasedDownsamplingUtilsUnitTest.class);

    /**
     * Directory where the testdata is located.
     */
    private static final File TEST_DATA_DIR = new File(CommandLineProgramTest.getTestDataDir(),"ArtificallyContaminatedBams/");

    @DataProvider(name = "oneCase")
    public Iterator<Object[]> oneCase() {

        final int[] idealHetAlleleCounts = {0, 50, 0, 50};
        final int[] idealHomAlleleCounts = {0, 100, 0, 0};

        final List<Object[]> data= new ArrayList<>();
        // no contamination, no removal
        data.add(new Object[]{0, 0, 0, 0, 0.1, 100, idealHetAlleleCounts, idealHetAlleleCounts});
        data.add(new Object[]{0, 0, 0, 0, 0.1, 100, idealHomAlleleCounts, idealHomAlleleCounts});

        // hom sample, het contaminant, different alleles
        data.add(new Object[]{5, 0, 0, 0, 0.1, 100, idealHomAlleleCounts, idealHomAlleleCounts});
        data.add(new Object[]{0, 0, 5, 0, 0.1, 100, idealHomAlleleCounts, idealHomAlleleCounts});
        data.add(new Object[]{0, 0, 0, 5, 0.1, 100, idealHomAlleleCounts, idealHomAlleleCounts});

        // hom sample, hom contaminant, different alleles
        data.add(new Object[]{10, 0, 0, 0, 0.1, 100, idealHomAlleleCounts, idealHomAlleleCounts});
        data.add(new Object[]{0, 0, 10, 0, 0.1, 100, idealHomAlleleCounts, idealHomAlleleCounts});
        data.add(new Object[]{0, 0, 0, 10, 0.1, 100, idealHomAlleleCounts, idealHomAlleleCounts});

        // het sample, het contaminant, different alleles
        data.add(new Object[]{5, 0, 0, 0, 0.1, 100, idealHetAlleleCounts, idealHetAlleleCounts});
        data.add(new Object[]{0, 0, 5, 0, 0.1, 100, idealHetAlleleCounts, idealHetAlleleCounts});

        // het sample, hom contaminant, different alleles
        data.add(new Object[]{10, 0, 0, 0, 0.1, 100, idealHetAlleleCounts, idealHetAlleleCounts});
        data.add(new Object[]{0, 0, 10, 0, 0.1, 100, idealHetAlleleCounts, idealHetAlleleCounts});

        // hom sample, het contaminant, overlapping alleles
        final int[] enhancedHomAlleleCounts = {0, 105, 0, 0};
        data.add(new Object[]{5, 5, 0, 0, 0.1, 100, idealHomAlleleCounts, enhancedHomAlleleCounts});
        data.add(new Object[]{0, 5, 5, 0, 0.1, 100, idealHomAlleleCounts, enhancedHomAlleleCounts});
        data.add(new Object[]{0, 5, 0, 5, 0.1, 100, idealHomAlleleCounts, enhancedHomAlleleCounts});

        // hom sample, hom contaminant, overlapping alleles
        data.add(new Object[]{0, 10, 0, 0, 0.1, 100, idealHomAlleleCounts, new int[]{0, 110, 0, 0}});

        // het sample, het contaminant, overlapping alleles
        data.add(new Object[]{5, 5, 0, 0, 0.1, 100, idealHetAlleleCounts, idealHetAlleleCounts});
        data.add(new Object[]{0, 5, 5, 0, 0.1, 100, idealHetAlleleCounts, idealHetAlleleCounts});
        data.add(new Object[]{0, 5, 0, 5, 0.1, 100, idealHetAlleleCounts, new int[]{0, 55, 0, 55}});
        data.add(new Object[]{5, 0, 0, 5, 0.1, 100, idealHetAlleleCounts, idealHetAlleleCounts});
        data.add(new Object[]{0, 0, 5, 5, 0.1, 100, idealHetAlleleCounts, idealHetAlleleCounts});

        // het sample, hom contaminant, overlapping alleles
        data.add(new Object[]{0, 10, 0, 0, 0.1, 100, idealHetAlleleCounts, idealHetAlleleCounts});
        data.add(new Object[]{0, 0, 0, 10, 0.1, 100, idealHetAlleleCounts, idealHetAlleleCounts});
        return data.iterator();
    }

    @Test(dataProvider= "oneCase")
    public void testOneCase(final int addA, final int addC, final int addG, final int addT, final double contaminationFraction,
                                    final int pileupSize, final int[] initialCounts, final int[] targetCounts) {

        final int[] actualCounts = initialCounts.clone();
        actualCounts[0] += addA;
        actualCounts[1] += addC;
        actualCounts[2] += addG;
        actualCounts[3] += addT;

        final int[] results = AlleleBiasedDownsamplingUtils.runSmartDownsampling(actualCounts, (int) (pileupSize * contaminationFraction));
        Assert.assertTrue(countsAreEqual(results, targetCounts));
    }

    @Test(dataProvider= "oneCase")
    public void testOneCaseUsingPublicAPI(final int addA, final int addC, final int addG, final int addT, final double contaminationFraction,
                                          final int pileupSize, final int[] initialCounts, final int[] targetCounts) {

        final int[] actualCounts = initialCounts.clone();
        actualCounts[0] += addA;
        actualCounts[1] += addC;
        actualCounts[2] += addG;
        actualCounts[3] += addT;

        final long actualCount = MathUtils.sum(actualCounts);
        final long expectedCount = MathUtils.sum(targetCounts);

        final SAMFileHeader header= ArtificialReadUtils.createArtificialSamHeader();

        final char[] bases= {'A', 'C', 'G', 'T'};   //note: hardwired to use same order as actualCounts
        final byte[] quals = {30};
        final Map<Allele, List<GATKRead>> readMap= new LinkedHashMap<>(bases.length);
        for (int idx = 0; idx < bases.length; idx++) {
            final Allele nonRefAllele = Allele.create(String.valueOf(bases[idx]));
            readMap.put(nonRefAllele, new ArrayList<>());
            for (int i = 0; i < actualCounts[idx]; i++) {
                final byte[] readBases = {(byte) bases[idx]};
                final GATKRead newRead = ArtificialReadUtils.createArtificialRead(header, "read" + i, 0, 1, readBases, quals, "1M");
                readMap.get(nonRefAllele).add(newRead);
            }
        }

        Assert.assertEquals(AlleleBiasedDownsamplingUtils.totalReads(readMap), actualCount);

        final List<GATKRead> results = AlleleBiasedDownsamplingUtils.selectAlleleBiasedReads(readMap, pileupSize, contaminationFraction);
        //TODO should add assertions to test that the reads are downsampled properly by allele.
        final long toRemove = actualCount - expectedCount;
        Assert.assertEquals(results.size(), toRemove);
    }

    private static boolean countsAreEqual(final int[] counts1, final int[] counts2) {
        for ( int i = 0; i < 4; i++ ) {
            if ( counts1[i] != counts2[i] )
                return false;
        }
        return true;
    }

    @Test
    public void testLoadContaminationFileDetails() throws IOException {
        final File ContamFile1=new File(TEST_DATA_DIR, "contamination.case.1.txt");

        final Map<String,Double> Contam1=new LinkedHashMap<>();
        final Set<String> Samples1= new LinkedHashSet<>();

        Contam1.put("NA11918",0.15);
        Samples1.addAll(Contam1.keySet());
        testLoadFile(ContamFile1,Samples1,Contam1,logger);

        Contam1.put("NA12842",0.13);
        Samples1.addAll(Contam1.keySet());
        testLoadFile(ContamFile1,Samples1,Contam1,logger);

        Samples1.add("DUMMY");
        testLoadFile(ContamFile1,Samples1,Contam1,logger);
   }

    private static void testLoadFile(final File file, final Set<String> Samples, final Map<String,Double> map, final Logger logger) throws IOException {
        final Map<String,Double> loadedMap = AlleleBiasedDownsamplingUtils.loadContaminationFile(file,0.0,Samples,logger);
        Assert.assertTrue(loadedMap.equals(map));
    }

    @DataProvider(name = "goodContaminationFiles")
    public Object[][] goodContaminationFiles() {
        return new Object[][]{
                {1, 2},
                {2, 3},
                {3, 2},
                {4, 2},
                {5, 3},
                {6, 2},
                {7, 2},
                {8, 2}
        };
    }

    @Test(dataProvider = "goodContaminationFiles")
    public void testLoadContaminationFile(final int artificalBAMnumber, final int numberOfSamples) throws IOException {
        final String ArtificialBAM = String.format("contamination.case.%d.txt", artificalBAMnumber);

        final File file = new File(TEST_DATA_DIR, ArtificialBAM);
        Assert.assertTrue(AlleleBiasedDownsamplingUtils.loadContaminationFile(file, 0.0, null, logger).size() == numberOfSamples);

    }


    @DataProvider(name = "badContaminationFiles")
    public Object[][] badContaminationFiles() {
        return new Object[][]{{1}, {2}, {3}, {4}, {5}};
    }

    @Test(dataProvider = "badContaminationFiles", expectedExceptions = UserException.MalformedFile.class)
    public void testLoadBrokenContaminationFile(final int i) throws IOException {
        final File file = new File(TEST_DATA_DIR, String.format("contamination.case.broken.%d.txt", i));
        AlleleBiasedDownsamplingUtils.loadContaminationFile(file, 0.0, null, logger);

    }


}
