package org.broadinstitute.hellbender.tools.walkers.validation;

import org.broadinstitute.hellbender.tools.walkers.validation.MixingFraction;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.Utils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;

/**
 * Created by David Benjamin on 1/31/17.
 */
public class MixingFractionUnitTest {

    @Test
    public void testMixingfraction() {
        final String sample = "SAMPLE";
        final double fraction = 0.15;
        final MixingFraction mixingFraction = new MixingFraction(sample, fraction);
        Assert.assertEquals(mixingFraction.getSample(), sample);
        Assert.assertEquals(mixingFraction.getMixingFraction(), fraction);
    }

    @Test
    public void testIO() throws IOException {
        Utils.resetRandomGenerator();
        final File file = File.createTempFile("mixing_fractions", ".table");

        final String sample1 = "SAMPLE1";
        final double fraction1 = 0.15;
        final String sample2 = "SAMPLE2";
        final double fraction2 = 0.17;
        final List<MixingFraction> original = Arrays.asList(new MixingFraction(sample1, fraction1), new MixingFraction(sample2, fraction2));

        MixingFraction.writeMixingFractions(original, file);
        final List<MixingFraction> copy = MixingFraction.readMixingFractions(file);

        Assert.assertEquals(original.size(), copy.size());
        new IndexRange(0, original.size()).forEach(n -> {
                    Assert.assertEquals(original.get(n).getSample(), copy.get(n).getSample());
                    Assert.assertEquals(original.get(n).getMixingFraction(), copy.get(n).getMixingFraction());
                });
    }
}