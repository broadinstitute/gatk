package org.broadinstitute.hellbender.tools.copynumber;

import org.testng.annotations.Test;

/**
 * Integration test for {@link CollectAllelicCountsSpark}.  Uses a BAM with sites generated from hg19mini using wgsim.
 */
@Test(groups = "spark")
public final class CollectAllelicCountsSparkIntegrationTest extends CollectAllelicCountsIntegrationTest {
}