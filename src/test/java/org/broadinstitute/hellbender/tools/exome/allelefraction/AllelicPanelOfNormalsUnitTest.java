package org.broadinstitute.hellbender.tools.exome.allelefraction;

import htsjdk.samtools.util.Log;
import org.broadinstitute.hellbender.utils.LoggingUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;

/**
 * Tests for the {@link AllelicPanelOfNormals}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public class AllelicPanelOfNormalsUnitTest extends BaseTest {

    private static final String TEST_SUB_DIR = publicTestDir + "org/broadinstitute/hellbender/tools/exome/";

    // test data is a PON generated from 50 normals simulated from the allele-fraction model with alpha = 65 and beta = 60
    private static final File ALLELIC_PON_NORMAL_FILE = new File(TEST_SUB_DIR, "allelic-pon-test-pon-normal.tsv");

    private static final double MLE_ALPHA_EXPECTED = 65;
    private static final double MLE_BETA_EXPECTED = 60;
    private static final double MLE_MEAN_BIAS_EXPECTED = meanBias(MLE_ALPHA_EXPECTED, MLE_BETA_EXPECTED);
    private static final double MLE_BIAS_VARIANCE_EXPECTED = biasVariance(MLE_ALPHA_EXPECTED, MLE_BETA_EXPECTED);
    private static final double ALPHA_EXPECTED_AT_FIRST_SITE = 698.6;
    private static final double BETA_EXPECTED_AT_FIRST_SITE = 645.9;
    private static final double DELTA = 3.;

    @Test
    public void testPoNHyperparameterInitialization() {
        LoggingUtils.setLoggingLevel(Log.LogLevel.INFO);
        final AllelicPanelOfNormals allelicPON = new AllelicPanelOfNormals(ALLELIC_PON_NORMAL_FILE);

        final SimpleInterval firstSite = new SimpleInterval("1", 1, 1);  //all sites in PON are from chr1
        final SimpleInterval siteNotInPON = new SimpleInterval("2", 1, 1);  //all sites in PON are from chr1

        // test initialization of hyperparameters for first site in PON (a = 1218, r = 1317)
        final double alphaAtFirstSite = allelicPON.getAlpha(firstSite);
        final double betaAtFirstSite = allelicPON.getBeta(firstSite);

        Assert.assertEquals(alphaAtFirstSite, ALPHA_EXPECTED_AT_FIRST_SITE, DELTA);
        Assert.assertEquals(betaAtFirstSite, BETA_EXPECTED_AT_FIRST_SITE, DELTA);

        // test initialization of MLE hyperparameters (which are default values for sites not in PON)
        final double alphaNotInPON = allelicPON.getAlpha(siteNotInPON);
        final double betaNotInPON = allelicPON.getBeta(siteNotInPON);
        final double meanBias = allelicPON.getMLEMeanBias();
        final double biasVariance = allelicPON.getMLEBiasVariance();

        Assert.assertEquals(alphaNotInPON, MLE_ALPHA_EXPECTED, DELTA);
        Assert.assertEquals(betaNotInPON, MLE_BETA_EXPECTED, DELTA);
        Assert.assertEquals(meanBias, MLE_MEAN_BIAS_EXPECTED, DELTA);
        Assert.assertEquals(biasVariance, MLE_BIAS_VARIANCE_EXPECTED, DELTA);
    }

    private static double meanBias(final double alpha, final double beta) {
        return alpha / beta;
    }

    private static double biasVariance(final double alpha, final double beta) {
        return alpha / (beta * beta);
    }
}