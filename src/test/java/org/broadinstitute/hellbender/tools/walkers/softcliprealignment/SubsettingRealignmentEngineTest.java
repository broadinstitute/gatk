package org.broadinstitute.hellbender.tools.walkers.softcliprealignment;

import com.google.api.client.util.Lists;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Iterator;
import java.util.List;

public class SubsettingRealignmentEngineTest extends CommandLineProgramTest {

    @Override
    public String getTestedClassName() {
        return SubsettingRealignmentEngine.class.getSimpleName();
    }

    private static final String BWA_IMAGE_PATH = "src/test/resources/large/SubsettingRealignmentEngine/hg38_chr22_test.fasta.img";
    private final String BAM_FILE = "src/test/resources/large/SubsettingRealignmentEngine/test.realigned.bam";

    @Test
    public void test() {
        final SamReader reader = SamReaderFactory.make().open(new File(BAM_FILE));
        final SAMFileHeader inputHeader = reader.getFileHeader();
        try (final SubsettingRealignmentEngine engine = new SubsettingRealignmentEngine(BWA_IMAGE_PATH, inputHeader, 10, 1)) {
            final Iterator<SAMRecord> iter = reader.iterator();
            while (iter.hasNext()) {
                final GATKRead read = new SAMRecordToGATKReadAdapter(iter.next());
                engine.addRead(read, r -> r.getMappingQuality() < 30);
            }
            final List<GATKRead> output = Lists.newArrayList(engine.alignAndMerge());

            Assert.assertEquals(engine.getSelectedReadsCount(), 103);
            Assert.assertEquals(engine.getNonselectedReadsCount(), 7713);
            Assert.assertEquals(engine.getPairedAlignmentReadsCount(), 34);
            Assert.assertEquals(engine.getUnpairedAlignmentReadsCount(), 69);

            final long realignedOutput = output.stream().filter(SubsettingRealignmentEngine::checkRealignedTag).count();
            Assert.assertEquals(realignedOutput, 105);

            final long notRealignedOutput = output.stream().filter(r -> !SubsettingRealignmentEngine.checkRealignedTag(r)).count();
            Assert.assertEquals(notRealignedOutput, 7713);

            final long properlyPairedOutput = output.stream().filter(GATKRead::isProperlyPaired).count();
            Assert.assertEquals(properlyPairedOutput, 7562);

            final long hasReadGroupOutput = output.stream().map(GATKRead::getReadGroup).filter(rg -> rg != null).count();
            Assert.assertEquals(hasReadGroupOutput, 7818);
        }
    }

}