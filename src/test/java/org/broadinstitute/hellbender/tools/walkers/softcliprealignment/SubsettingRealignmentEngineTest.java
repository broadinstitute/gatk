package org.broadinstitute.hellbender.tools.walkers.softcliprealignment;

import com.google.api.client.util.Lists;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.tools.spark.bwa.BwaSpark;
import org.broadinstitute.hellbender.tools.spark.bwa.BwaSparkIntegrationTest;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Iterator;
import java.util.List;
import java.util.stream.Collectors;

public class SubsettingRealignmentEngineTest extends CommandLineProgramTest {

    @Override
    public String getTestedClassName() {
        return SubsettingRealignmentEngine.class.getSimpleName();
    }

    private static final String BWA_IMAGE_PATH = "src/test/resources/" + BwaSparkIntegrationTest.class.getPackage().getName().replace(".","/") +"/" + BwaSpark.class.getSimpleName() + "/ref.fa.img";
    private final File BAM_FILE = getTestFile("test.bam");
    @Test
    public void test() {
        final SamReader reader = SamReaderFactory.make().open(BAM_FILE);
        final SAMFileHeader inputHeader = reader.getFileHeader();
        try (final SubsettingRealignmentEngine engine = new SubsettingRealignmentEngine(BWA_IMAGE_PATH, inputHeader, 100, 1)) {
            final Iterator<SAMRecord> iter = reader.iterator();
            while (iter.hasNext()) {
                final GATKRead read = new SAMRecordToGATKReadAdapter(iter.next());
                engine.addRead(read, r -> r.getMappingQuality() < 30);
            }
            final List<GATKRead> output = Lists.newArrayList(engine.alignAndMerge());
            final List<GATKRead> realignedOutput = output.stream().filter(SubsettingRealignmentEngine::checkRealignedTag).collect(Collectors.toList());
            final List<GATKRead> notRealignedOutput = output.stream().filter(r -> !SubsettingRealignmentEngine.checkRealignedTag(r)).collect(Collectors.toList());
            Assert.assertEquals(realignedOutput.size(), 89);
            Assert.assertEquals(notRealignedOutput.size(), 474);
            final List<GATKRead> properlyPairedOutput = output.stream().filter(GATKRead::isProperlyPaired).collect(Collectors.toList());
            Assert.assertEquals(properlyPairedOutput.size(), 0);
        }
    }

}