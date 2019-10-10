package org.broadinstitute.hellbender.engine.filters;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.tools.walkers.variantutils.ValidateVariants;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.annotations.Test;

public class UMIReadFilterTest extends GATKBaseTest {
    @Test
    public void test(){
        

        int d = 3;

        final ArgumentsBuilder validateVariantsArgs = new ArgumentsBuilder()
                .addArgument("R", MITO_REF.getAbsolutePath())
                .addArgument("V", standardVcf.getAbsolutePath())
                .addArgument("L", IntervalUtils.locatableToString(new SimpleInterval("chrM:1-1000")))
                .add("-gvcf");
        runCommandLine(validateVariantsArgs, ValidateVariants.class.getSimpleName());
    }
}