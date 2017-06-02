package org.broadinstitute.hellbender.tools.coveragemodel;

import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.hellbender.tools.exome.sexgenotyper.SexGenotypeData;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

/**
 * Unit tests for {@link CoverageModelCopyRatioEmissionProbabilityCalculator}
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public class CoverageModelCopyRatioEmissionProbabilityCalculatorUnitTest extends BaseTest {

    private static final double[] truthReadCountsTable = {1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8};
    private static final double[] truthMuTable = {-20.0741, 22.4479, 9.7571, 16.9622, -19.2745};
    private static final double[] truthPsiTable = {0.529822, 0.523619, 0.181005, 0.308569, 0.121021};

    /**
     * log normalization constant for (mu, psi)
     * these are calculated accurately using Mathematica via brute-force summation
     */
    private static final double[][] truthLogProbabilityMassTable = {
            {-132.836567457, -133.370722005, -171.607050746, -155.027162677, -180.709888717},
            {1.1614683835e-10, 1.1578717714e-10, 9.755796876e-11, 1.0398313589e-10, 9.467541922e-11},
            {0.0000377381250187, 0.0000376210877616, 0.0000316917381709, 0.0000337810883095, 0.0000307545504235},
            {2.8018462e-8, 2.7931699e-8, 2.3534208e-8, 2.5084171e-8, 2.2838841e-8},
            {-122.552919532, -123.045208201, -158.286065831, -143.004754512, -166.676111059}
    };

    /**
     * log probability density for (read count, mu, psi)
     * these are calculated accurately using Mathematica
     */
    private static final double[][][] truthNormalizedLogProbabilityDensity = {
            {
                {0.0, 0.0, 0.0, 0.0, 0.0},
                {-165.826976262, -166.495438446, -214.341374329, -193.595593167, -225.730883739},
                {-32.2465932017, -32.3712344664, -41.3072080963, -37.4294149773, -43.4378275750},
                {-95.1669811037, -95.5477785606, -122.811526495, -110.988398813, -129.303364166},
                {0.0, 0.0, 0.0, 0.0, 0.0}
            }, {
                {-267.660662513, -271.075357709, -721.91921978, -460.516092201, -954.49618457},
                {-325.172485776, -328.372118637, -724.699537888, -499.427803360, -920.559974210},
                {-47.1058387174, -47.5396837102, -101.463601762, -70.7793365351, -128.178566238},
                {-173.596705671, -175.288691414, -384.969957945, -265.768992624, -488.628206408},
                {-250.042203165, -253.221236184, -672.702699190, -429.526528418, -889.01985508}
            }, {
                {-436.516664404, -442.534298939, -1427.46383024, -805.86765233, -2148.11972699},
                {-300.094481897, -303.516390718, -838.087299200, -504.630196675, -1219.44743552},
                {-29.8002426402, -30.0802340536, -74.1771253018, -46.6108897736, -105.798775784},
                {-146.647005245, -148.285231165, -404.411474633, -244.610218179, -587.223241368},
                {-410.835543900, -416.476987634, -1339.14058520, -756.94700352, -2013.97383281}
            }, {
                {-560.422879828, -567.990789769, -1835.37073100, -1028.07642395, -2809.25963892},
                {-234.984519227, -237.668181256, -670.409082464, -397.293784577, -996.350191741},
                {-15.1574059035, -15.2419468439, -29.2784928360, -20.3534284649, -40.0429229288},
                {-102.731350666, -103.851320574, -284.689867462, -170.517607539, -421.013296838},
                {-530.663375260, -537.799701841, -1731.90499911, -971.43633323, -2649.09492089}
            }, {
                {-686.131797723, -695.175240178, -2205.29705035, -1243.67162040, -3368.55658632},
                {-175.150913661, -177.103269366, -493.064836081, -293.395211830, -732.459650602},
                {-10.0938681174, -10.0913216056, -10.1003042586, -10.0258166199, -10.3079013560},
                {-66.5095293450, -67.1751588692, -175.175184076, -106.880059435, -257.135795835},
                {-652.829786859, -661.398875626, -2091.08440718, -1180.86632050, -3191.89596631}
            }, {
                {-820.84178953, -831.45535476, -2595.93141636, -1473.48762657, -3952.52642283},
                {-124.955951844, -126.286767125, -341.865261052, -205.593780903, -505.360079590},
                {-15.0236497667, -15.0522251523, -20.0929724658, -16.8393718830, -24.1122036100},
                {-40.1365601831, -40.4626193929, -93.5986225762, -59.9584106279, -134.048464038},
                {-784.05728628, -794.15515609, -2471.48376795, -1404.69375571, -3760.51119557}
            }, {
                {-965.43957054, -977.73896112, -3014.83635749, -1720.13452602, -4578.02318849},
                {-84.7409249182, -85.5680981255, -219.725393967, -134.894858840, -321.552099611},
                {-29.9604732155, -30.1387148021, -59.3775040899, -40.8353400111, -81.7272653271},
                {-23.7609444547, -23.8657467806, -41.2309477208, -30.1906157715, -54.5863287779},
                {-925.17912002, -936.92162940, -2880.20975711, -1645.37137721, -4370.77870242}
            }, {
                {-1120.03038409, -1134.13381126, -3462.91352441, -1983.92262176, -5247.06324355},
                {-54.5309722786, -54.9730016502, -126.860921297, -81.3726338101, -181.518188265},
                {-54.9042389829, -55.3506895279, -127.953434738, -82.0135124012, -183.152196755},
                {-17.3917807348, -17.3938572490, -18.1504710745, -17.6035796694, -18.9246662952},
                {-1076.29469407, -1089.80007114, -3318.11403600, -1903.19228164, -5024.60306319}
            }, {
                {-1284.62653105, -1300.65249932, -3940.26834414, -2264.88818689, -5959.88243320},
                {-34.3278397751, -34.5032652389, -63.2868395916, -45.0322612933, -85.2919065157},
                {-89.8549273734, -90.6881292471, -225.820634452, -140.373839455, -328.386722080},
                {-21.0295384778, -21.0474315229, -24.3612535523, -22.1986948101, -27.0725754712},
                {-1237.41567762, -1252.80242870, -3785.29661965, -2178.19087977, -5722.20801694}
            }
    };

    /* By setting the switching threshold read count to 0, this instance always performs calculations
     * in the Laplace approximation */
    private static final CoverageModelCopyRatioEmissionProbabilityCalculator emissionCalculatorLaplace =
            new CoverageModelCopyRatioEmissionProbabilityCalculator(0, false);

    @Test(dataProvider = "truthLogProbabilityMassTestData")
    public void testGetLogProbabilityMass(final double mu, final double psi,
                                          final double expectedLogProbabilityMass) {
        final double absAccuracy = FastMath.max(CoverageModelCopyRatioEmissionProbabilityCalculator.ABSOLUTE_ACCURACY,
                FastMath.abs(expectedLogProbabilityMass) * CoverageModelCopyRatioEmissionProbabilityCalculator.RELATIVE_ACCURACY);
        Assert.assertEquals(CoverageModelCopyRatioEmissionProbabilityCalculator
                .getLogProbabilityLaplaceApproximationNormalizationConstant(mu, psi), expectedLogProbabilityMass,
                absAccuracy);
    }

    @Test(dataProvider = "truthNormalizedLogProbabilityDensity")
    public void testLogLikelihood(final double readCounts, final double mu, final double psi,
                                  final double expectedLogLikelihood) {
        final double absAccuracy = FastMath.max(CoverageModelCopyRatioEmissionProbabilityCalculator.ABSOLUTE_ACCURACY,
                FastMath.abs(expectedLogLikelihood) * CoverageModelCopyRatioEmissionProbabilityCalculator.RELATIVE_ACCURACY);

        final CoverageModelCopyRatioEmissionData emissionData =
                new CoverageModelCopyRatioEmissionData(mu, psi, (int)readCounts, 0);
        emissionData.setCopyRatioCallingMetadata(CopyRatioCallingMetadata.builder()
                .sampleCoverageDepth(1.0)
                .sampleSexGenotypeData(new SexGenotypeData("N/A", "N/A", null, null))
                .emissionCalculationStrategy(CoverageModelCopyRatioEmissionProbabilityCalculator.EmissionCalculationStrategy.HYBRID_POISSON_GAUSSIAN)
                .build());
        Assert.assertEquals(emissionCalculatorLaplace.logLikelihood(emissionData, 1.0, null), expectedLogLikelihood, absAccuracy);
    }

    @DataProvider(name = "truthLogProbabilityMassTestData")
    public Object[][] getTruthLogProbabilityMassTestData() {
        final Double[][] truthData = new Double[truthMuTable.length * truthPsiTable.length][];
        for(int i = 0; i < truthMuTable.length; i++) {
            for (int j = 0; j < truthPsiTable.length; j++) {
                truthData[truthPsiTable.length * i + j] = new Double[] {
                        truthMuTable[i], truthPsiTable[j], truthLogProbabilityMassTable[i][j]};
            }
        }
        return truthData;
    }

    @DataProvider(name = "truthNormalizedLogProbabilityDensity")
    public Object[][] getTruthNormalizedLogProbabilityDensity() {
        final Double[][] truthData = new Double[truthReadCountsTable.length * truthMuTable.length *
                truthPsiTable.length][];
        int counter = 0;
        for (int k = 0; k < truthReadCountsTable.length; k++) {
            for (int i = 0; i < truthMuTable.length; i++) {
                for (int j = 0; j < truthPsiTable.length; j++) {
                    truthData[counter] = new Double[]{
                            truthReadCountsTable[k], truthMuTable[i], truthPsiTable[j],
                            truthNormalizedLogProbabilityDensity[k][i][j]};
                    counter++;
                }
            }
        }
        return truthData;
    }

}
