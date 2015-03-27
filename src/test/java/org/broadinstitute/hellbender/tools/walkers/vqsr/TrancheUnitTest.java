package org.broadinstitute.hellbender.tools.walkers.vqsr;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.text.XReadLines;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public final class TrancheUnitTest extends BaseTest {
    private final String testDir = CommandLineProgramTest.getTestDataDir() + "/TrancheManagerUnitTest/";
    private final File QUAL_DATA = new File(testDir + "tranches.raw.dat");
    private final double[] TRUTH_SENSITIVITY_CUTS = new double[]{99.9, 99.0, 97.0, 95.0};
    private final File EXPECTED_TRANCHES_NEW = new File(testDir + "tranches.6.txt");
    private final File EXPECTED_TRANCHES_OLD = new File(testDir + "tranches.4.txt");

    private ArrayList<VariantDatum> readData() throws IOException{
        ArrayList<VariantDatum> vd = new ArrayList<>();
        try (XReadLines xrl = new XReadLines(QUAL_DATA, true)){
            for ( String line : xrl ) {
                String[] parts = line.split("\t");
                // QUAL,TRANSITION,ID,LOD,FILTER
                if ( ! parts[0].equals("QUAL") ) {
                    VariantDatum datum = new VariantDatum();
                    datum.lod = Double.valueOf(parts[3]);
                    datum.isTransition = parts[1].equals("1");
                    datum.isKnown = ! parts[2].equals(".");
                    datum.isSNP = true;
                    datum.atTruthSite = datum.isKnown;
                    vd.add(datum);
                }
            }
        }

        return vd;
    }

    @Test(expectedExceptions = {UserException.MalformedFile.class})
    public final void readBadFormat() throws IOException{
        Tranche.readTranches(QUAL_DATA);
    }

    @Test
    public final void readNewFormat() throws IOException{
        read(EXPECTED_TRANCHES_NEW);
    }

    @Test(expectedExceptions = {UserException.MalformedFile.class})
    public final void readOldFormat() throws IOException {
        read(EXPECTED_TRANCHES_OLD);
    }

    public final List<Tranche> read(File f) throws IOException{
        return Tranche.readTranches(f);
    }

    private static void assertTranchesAreTheSame(List<Tranche> newFormat, List<Tranche> oldFormat) {
        Assert.assertEquals(oldFormat.size(), newFormat.size());
        for ( int i = 0; i < newFormat.size(); i++ ) {
            Tranche n = newFormat.get(i);
            Tranche o = oldFormat.get(i);
            Assert.assertEquals(n.targetTruthSensitivity, o.targetTruthSensitivity, 1e-3);
            Assert.assertEquals(n.numNovel, o.numNovel);
            Assert.assertEquals(n.novelTiTv, o.novelTiTv, 1e-3);
            Assert.assertEquals(n.numKnown, o.numKnown);
            Assert.assertEquals(n.knownTiTv, o.knownTiTv, 1e-3);
        }
    }

    private static List<Tranche> findMyTranches(ArrayList<VariantDatum> vd, double[] tranches) {
        final int nCallsAtTruth = VariantDatum.countCallsAtTruth( vd, Double.NEGATIVE_INFINITY );
        final Tranche.TruthSensitivityMetric metric = new Tranche.TruthSensitivityMetric( nCallsAtTruth );
        return Tranche.findTranches(vd, tranches, metric, VariantRecalibratorArgumentCollection.Mode.SNP);
    }

    @Test
    public final void testFindTranches1() throws IOException {
        ArrayList<VariantDatum> vd = readData();
        List<Tranche> tranches = findMyTranches(vd, TRUTH_SENSITIVITY_CUTS);
        tranches.sort(Tranche.TRUTH_SENSITIVITY_ORDER);
        assertTranchesAreTheSame(read(EXPECTED_TRANCHES_NEW), tranches);
    }

    @Test(expectedExceptions = {UserException.class})
    public final void testBadFDR() throws IOException {
        ArrayList<VariantDatum> vd = readData();
        findMyTranches(vd, new double[]{-1});
    }
}
