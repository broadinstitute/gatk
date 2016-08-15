package org.broadinstitute.hellbender.tools.pon.allelic;

import htsjdk.samtools.util.Log;
import org.broadinstitute.hellbender.tools.exome.allelefraction.AlleleFractionLikelihoods;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCountCollection;
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
public final class AllelicPanelOfNormalsUnitTest extends BaseTest {
    private static final String TEST_SUB_DIR = publicTestDir + "org/broadinstitute/hellbender/tools/exome/";

    // test data is a PoN with counts generated from 50 normals simulated from the allele-fraction model with alpha = 65 and beta = 60
    private static final File ALLELIC_PON_NORMAL_COUNTS_FILE = new File(TEST_SUB_DIR, "allelic-pon-test-pon-counts-normal.tsv");

    private static final double MLE_ALPHA_EXPECTED = 65;
    private static final double MLE_BETA_EXPECTED = 60;
    private static final double MLE_MEAN_BIAS_EXPECTED = AlleleFractionLikelihoods.meanBias(MLE_ALPHA_EXPECTED, MLE_BETA_EXPECTED);
    private static final double MLE_BIAS_VARIANCE_EXPECTED = AlleleFractionLikelihoods.biasVariance(MLE_ALPHA_EXPECTED, MLE_BETA_EXPECTED);
    private static final double ALPHA_EXPECTED_AT_FIRST_SITE = 698.6;
    private static final double BETA_EXPECTED_AT_FIRST_SITE = 645.9;
    private static final double DELTA = 3.;

    @Test
    public void testPoNHyperparameterInitialization() {
        LoggingUtils.setLoggingLevel(Log.LogLevel.INFO);
        final AllelicPanelOfNormals allelicPoN = new AllelicPanelOfNormals(new AllelicCountCollection(ALLELIC_PON_NORMAL_COUNTS_FILE));

        final SimpleInterval firstSite = new SimpleInterval("1", 1, 1);
        final SimpleInterval siteNotInPoN = new SimpleInterval("2", 1, 1);  //all sites in PoN are from chr1

        // test initialization of hyperparameters for first site in PoN (a = 1218, r = 1317)
        final double alphaAtFirstSite = allelicPoN.getAlpha(firstSite);
        final double betaAtFirstSite = allelicPoN.getBeta(firstSite);

        Assert.assertEquals(alphaAtFirstSite, ALPHA_EXPECTED_AT_FIRST_SITE, DELTA);
        Assert.assertEquals(betaAtFirstSite, BETA_EXPECTED_AT_FIRST_SITE, DELTA);

        // test initialization of MLE hyperparameters (which are default values for sites not in PoN)
        final double alphaNotInPoN = allelicPoN.getAlpha(siteNotInPoN);
        final double betaNotInPoN = allelicPoN.getBeta(siteNotInPoN);
        final double meanBias = allelicPoN.getGlobalMeanBias();
        final double biasVariance = allelicPoN.getGlobalBiasVariance();

        Assert.assertEquals(alphaNotInPoN, MLE_ALPHA_EXPECTED, DELTA);
        Assert.assertEquals(betaNotInPoN, MLE_BETA_EXPECTED, DELTA);
        Assert.assertEquals(meanBias, MLE_MEAN_BIAS_EXPECTED, DELTA);
        Assert.assertEquals(biasVariance, MLE_BIAS_VARIANCE_EXPECTED, DELTA);
    }
}