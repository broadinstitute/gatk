package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVFastqUtils;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAligner;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndexCache;
import org.broadinstitute.hellbender.utils.fermi.FermiLiteAssembler;
import org.broadinstitute.hellbender.utils.fermi.FermiLiteAssembly;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.List;

public class ReviseAssemblyUnitTest extends GATKBaseTest {
    @Test(groups = "sv")
    void testReviseAssembly() {
        final String fastqFile = getToolTestDataDir() + "/test.fastq";
        final List<SVFastqUtils.FastqRead> readList = SVFastqUtils.readFastqFile(fastqFile);
        final FermiLiteAssembly initialAssembly = new FermiLiteAssembler().createAssembly(readList);
        Assert.assertEquals(initialAssembly.getNContigs(), 9);
        final FermiLiteAssembly revisedAssembly =
                FermiLiteAssemblyHandler.reviseAssembly(initialAssembly, true, true);
        Assert.assertEquals(revisedAssembly.getNContigs(), 6);
    }
}
