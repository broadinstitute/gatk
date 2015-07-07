package org.broadinstitute.hellbender.tools.exome;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Unit tests for {@link SampleCollection}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class SampleCollectionUnitTest extends BaseTest {

    private static SAMFileHeader emptyHeader;
    private static SAMFileHeader singleSampleReadGroupHeader;
    private static SAMFileHeader singleSampleLessReadGroupHeader;
    private static SAMFileHeader[] multiSampleHeaders;

    private final static int minMultiSampleCount = 2;
    private final static int maxMultiSampleCount = 101;
    private final static int minReadGroupPerSample = 1;
    private final static int maxReadGroupPerSample = 101;
    private final static int minOrphanReadGroupCount = 1;
    private final static int maxOrphanReadGroupCount = 10;

    private final static int randomSeed = 131113;
    private final static int multiSampleWithOrphanReadGroupCount = 2;
    private final static int multiSampleWithoutOrphanReadGroupCount = 10;


    @BeforeClass
    public void setUp() {
        final Random rnd = new Random(randomSeed);
        emptyHeader = ArtificialReadUtils.createArtificialSamHeader();
        final SAMReadGroupRecord rg1_1 = new SAMReadGroupRecord("RG1");
        rg1_1.setSample("SM1");
        final SAMReadGroupRecord rg1_0 = new SAMReadGroupRecord("RG0");
        singleSampleReadGroupHeader = ArtificialReadUtils.createArtificialSamHeader();
        singleSampleReadGroupHeader.addReadGroup(rg1_1);
        singleSampleLessReadGroupHeader = ArtificialReadUtils.createArtificialSamHeader();
        singleSampleLessReadGroupHeader.addReadGroup(rg1_0);
        multiSampleHeaders = new SAMFileHeader[multiSampleWithOrphanReadGroupCount + multiSampleWithoutOrphanReadGroupCount];
        for (int i = 0; i < multiSampleWithOrphanReadGroupCount; i++) {
            final int orphanCount = rnd.nextInt(maxOrphanReadGroupCount - minOrphanReadGroupCount) + minOrphanReadGroupCount;
            final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();
            for (int j = 0; j < orphanCount; j++) {
                header.addReadGroup(new SAMReadGroupRecord("RG0_" + j));
            }
            addSampleReadGroups(rnd, header);
            multiSampleHeaders[i] = header;
        }

        for (int i = 0; i < multiSampleWithoutOrphanReadGroupCount; i++) {
            final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();
            addSampleReadGroups(rnd, header);
            multiSampleHeaders[i + multiSampleWithOrphanReadGroupCount] = header;
        }
    }

    private void addSampleReadGroups(Random rnd, SAMFileHeader header) {
        final int sampleCount = rnd.nextInt(maxMultiSampleCount - minMultiSampleCount) + minMultiSampleCount;

        for (int j = 0; j < sampleCount; j++) {
            final int readGroupCount = rnd.nextInt(maxReadGroupPerSample - minReadGroupPerSample) + minReadGroupPerSample;
            for (int k = 0; k < readGroupCount; k++) {
                final SAMReadGroupRecord rg = new SAMReadGroupRecord("RG" + (j + 1) + "_" + (k + 1));
                rg.setSample("SM" + (j+1));
                header.addReadGroup(rg);
            }
        }
    }

    @Test(dataProvider = "samFileHeaderData")
    public void testCreation(final SAMFileHeader header) {
        new SampleCollection(header);
    }

    @Test(dataProvider = "samFileHeaderData")
    public void testSampleCount(final SAMFileHeader header) {
        final SampleCollection sampleCollection = new SampleCollection(header);
        final Set<String> sampleIds = header.getReadGroups().stream().map(SAMReadGroupRecord::getSample).filter(Objects::nonNull).collect(Collectors.toSet());
        Assert.assertEquals(sampleCollection.sampleCount(), sampleIds.size(), header.getReadGroups().toString() + " " + sampleIds);
    }

    @Test(dataProvider = "samFileHeaderData", dependsOnMethods = "testCreation")
    public void testSamples(final SAMFileHeader header) {
        final SampleCollection sampleCollection = new SampleCollection(header);
        final List<Sample> samples = sampleCollection.samples();
        final List<String> sampleIds = header.getReadGroups().stream().map(SAMReadGroupRecord::getSample).filter(Objects::nonNull).sorted().distinct().collect(Collectors.toList());
        Assert.assertEquals(samples.size(),sampleIds.size());
        for (int i = 0; i < samples.size();i++) {
            Assert.assertEquals(samples.get(i).getId(),sampleIds.get(i));
            final int idx = i;
            final List<String> expectedGroups = header.getReadGroups().stream().filter(g -> sampleIds.get(idx).equals(g.getSample())).map(SAMReadGroupRecord::getId).sorted().collect(Collectors.toList());
            final Sample sample = samples.get(i);
            final List<String> observedGroups =  new ArrayList<>(sample.readGroups());
            Assert.assertEquals(expectedGroups,observedGroups);
        }
    }

    @Test(dataProvider = "samFileHeaderData", dependsOnMethods = "testCreation")
    public void testSampleById(final SAMFileHeader header) {
        final SampleCollection sampleCollection = new SampleCollection(header);
        final List<Sample> samples = sampleCollection.samples();
        for (int i = 0; i < samples.size();i++) {
            final Sample sample = samples.get(i);
            Assert.assertEquals(sampleCollection.sample(sample.getId()),sample);
            Assert.assertNull(sampleCollection.sample(sample.getId() + "_XXX"));
        }
    }

    @Test(dataProvider = "samFileHeaderData")
    public void testSampleIds(final SAMFileHeader header) {
        final SampleCollection sampleCollection = new SampleCollection(header);
        final Set<String> sampleIds = header.getReadGroups().stream().map(SAMReadGroupRecord::getSample).filter(Objects::nonNull).collect(Collectors.toSet());
        for (final SAMReadGroupRecord rg : header.getReadGroups()) {
            final String sampleId = rg.getSample();
            if (sampleId == null) {
                continue;
            }
            sampleIds.add(rg.getSample());
        }
        final List<String> observedSampleIds = sampleCollection.sampleIds();
        Assert.assertEquals(observedSampleIds.size(),sampleIds.size(), Arrays.toString(observedSampleIds.toArray()));
        for (final String sampleId : observedSampleIds) {
            Assert.assertTrue(sampleIds.contains(sampleId), "Missing expected sample id: " + sampleId);
        }
        final List<String> sortedSampleIds = new ArrayList<>(sampleIds);
        Collections.sort(sortedSampleIds);
        Assert.assertEquals(observedSampleIds,sortedSampleIds, Arrays.toString(observedSampleIds.toArray()) + " " + Arrays.toString(sortedSampleIds.toArray()));
    }

    @Test(dataProvider = "samFileHeaderData", dependsOnMethods = "testCreation")
    public void testReadGroupCount(final SAMFileHeader header) {
        final SampleCollection sampleCollection = new SampleCollection(header);
        Assert.assertEquals(sampleCollection.readGroupCount(),header.getReadGroups().size());
    }

    @Test(dataProvider = "samFileHeaderData", dependsOnMethods = "testCreation")
    public void testReadGroupIndexByReadWithNullReadGroup(final SAMFileHeader header) {
        final SampleCollection sampleCollection = new SampleCollection(header);
        final GATKRead read = ArtificialReadUtils.createRandomRead(header, 100);
        read.setReadGroup(null);
        Assert.assertEquals(sampleCollection.readGroupIndexByRead(read), -1);
    }

    @Test(dataProvider = "samFileHeaderData", dependsOnMethods = "testCreation", expectedExceptions = IllegalArgumentException.class)
    public void testReadGroupIndexByReadWithNonExistentReadGroup(final SAMFileHeader header) {
        final SampleCollection sampleCollection = new SampleCollection(header);
        final GATKRead read = ArtificialReadUtils.createRandomRead(header, 100);
        read.setReadGroup("NON_EXISTENT");
        Assert.assertEquals(sampleCollection.readGroupIndexByRead(read), -1);
    }

    @Test(dataProvider = "samFileHeaderData", dependsOnMethods = "testCreation")
    public void testReadGroupIndexByReadWithReadGroupWithSample(final SAMFileHeader header) {
        final SampleCollection sampleCollection = new SampleCollection(header);
        final GATKRead read = ArtificialReadUtils.createRandomRead(header, 100);
        if (sampleCollection.sampleCount() > 0) {
            final String readGroup = sampleCollection.samples().get(0).readGroups().stream().findFirst().get();
            read.setReadGroup(readGroup);
            Assert.assertEquals(sampleCollection.readGroupIndexByRead(read), sampleCollection.readGroupIndexById(readGroup));
        }
    }

    @Test(dataProvider = "samFileHeaderData", dependsOnMethods = "testCreation")
    public void testReadGroupIndexByReadWithReadGroupWithoutSample(final SAMFileHeader header) {
        final String readGroup = header.getReadGroups().stream().filter(rg -> rg.getSample() == null).map(SAMReadGroupRecord::getId).findFirst().orElse(null);
        if (readGroup == null)
            return;
        final SampleCollection sampleCollection = new SampleCollection(header);
        final GATKRead read = ArtificialReadUtils.createRandomRead(header, 100);
        read.setReadGroup(readGroup);
        Assert.assertEquals(sampleCollection.readGroupIndexByRead(read), sampleCollection.readGroupIndexById(readGroup));

    }

    @Test(dataProvider = "samFileHeaderData", dependsOnMethods = "testCreation")
    public void testSampleIndexByReadWithNullReadGroup(final SAMFileHeader header) {
        final SampleCollection sampleCollection = new SampleCollection(header);
        final GATKRead read = ArtificialReadUtils.createRandomRead(header, 100);
        read.setReadGroup(null);
        Assert.assertEquals(sampleCollection.sampleIndexByRead(read), -1);
    }

    @Test(dataProvider = "samFileHeaderData", dependsOnMethods = "testCreation", expectedExceptions = IllegalArgumentException.class)
    public void testSampleIndexByReadWithNonExistentReadGroup(final SAMFileHeader header) {
        final SampleCollection sampleCollection = new SampleCollection(header);
        final GATKRead read = ArtificialReadUtils.createRandomRead(header, 100);
        read.setReadGroup("NON_EXISTENT");
        Assert.assertEquals(sampleCollection.sampleIndexByRead(read), -1);
    }

    @Test(dataProvider = "samFileHeaderData", dependsOnMethods = "testCreation")
    public void testSampleIndexByReadWithReadGroupWithSample(final SAMFileHeader header) {
        final SampleCollection sampleCollection = new SampleCollection(header);
        final GATKRead read = ArtificialReadUtils.createRandomRead(header, 100);
        if (sampleCollection.sampleCount() > 0) {
            read.setReadGroup(sampleCollection.samples().get(0).readGroups().stream().findFirst().get());
            Assert.assertEquals(sampleCollection.sampleIndexByRead(read), 0);
        }
    }

    @Test(dataProvider = "samFileHeaderData", dependsOnMethods = "testCreation")
    public void testSampleIndexByReadWithReadGroupWithoutSample(final SAMFileHeader header) {
        final String readGroup = header.getReadGroups().stream().filter(rg -> rg.getSample() == null).map(SAMReadGroupRecord::getId).findFirst().orElse(null);
        if (readGroup == null)
            return;
        final SampleCollection sampleCollection = new SampleCollection(header);
        final GATKRead read = ArtificialReadUtils.createRandomRead(header, 100);
        if (sampleCollection.sampleCount() > 0) {
            read.setReadGroup(readGroup);
            Assert.assertEquals(sampleCollection.sampleIndexByRead(read), -1);
        }
    }

    @Test(dataProvider = "samFileHeaderData", dependsOnMethods = "testCreation")
    public void testReadGroupIds(final SAMFileHeader header) {
        final SampleCollection sampleCollection = new SampleCollection(header);
        final List<String> readGroupIds = sampleCollection.readGroups();
        Assert.assertEquals(readGroupIds.size(),header.getReadGroups().size());

        final List<String> expectedReadGroupIds = new ArrayList<>();
        header.getReadGroups().forEach(rg -> expectedReadGroupIds.add(rg.getId()));
        Collections.sort(expectedReadGroupIds, (a,b) -> {
            final SAMReadGroupRecord arg = header.getReadGroup(a);
            final SAMReadGroupRecord brg = header.getReadGroup(b);
            final String aSampleId = arg.getSample();
            final String bSampleId = brg.getSample();

            if (Objects.equals(aSampleId,bSampleId)) {
                return arg.getId().compareTo(brg.getId());
            } else if (aSampleId == null) {
                return 1;
            } else if (bSampleId == null) {
                return -1;
            } else {
                return aSampleId.compareTo(bSampleId);
            }
        });

        Assert.assertEquals(readGroupIds,expectedReadGroupIds);
    }


    @Test(dataProvider = "samFileHeaderData", dependsOnMethods = {"testSampleIds"})
    public void testSampleIndexByGroupId(final SAMFileHeader header) {
        final SampleCollection sampleCollection = new SampleCollection(header);
        final List<String> sampleIds = sampleCollection.sampleIds();

        for (final SAMReadGroupRecord rg : header.getReadGroups()) {
            final String sampleId = rg.getSample();
            if (sampleId == null) {
                Assert.assertEquals(sampleCollection.sampleIndexByGroupId(rg.getId()), -1);
            } else {
                Assert.assertEquals(sampleCollection.sampleIndexByGroupId(rg.getId()), sampleIds.indexOf(sampleId));
            }
            Assert.assertEquals(sampleCollection.sampleIndexByGroupId(rg.getId() + "random_extension"),-1);
        }
    }

    @Test(dataProvider = "samFileHeaderData", dependsOnMethods = {"testReadGroupIds"})
    public void testReadGroupIndexByGroupId(final SAMFileHeader header) {
        final SampleCollection sampleCollection = new SampleCollection(header);
        final List<String> readGroupIds = sampleCollection.readGroups();

        for (final SAMReadGroupRecord rg : header.getReadGroups()) {
            Assert.assertTrue(readGroupIds.indexOf(rg.getId()) >= 0);
            Assert.assertEquals(sampleCollection.readGroupIndexById(rg.getId()), readGroupIds.indexOf(rg.getId()));
            Assert.assertEquals(sampleCollection.readGroupIndexById(rg.getId() + "_random_extension"), -1);
        }

    }

    @Test(expectedExceptions = {IllegalArgumentException.class},
          dataProvider = "wrongCreationParamData")
    public void testWrongCreation(final SAMFileHeader header) {
        new SampleCollection(header);
    }

    @DataProvider(name = "samFileHeaderData")
    public Object[][] samFileHeaderData() {
       final List<Object[]> result = new ArrayList<>(multiSampleHeaders.length + 3);
       result.add(new Object[] { emptyHeader });
       result.add(new Object[] { singleSampleLessReadGroupHeader});
       result.add(new Object[] { singleSampleReadGroupHeader });
       for (final SAMFileHeader header : multiSampleHeaders) {
           result.add(new Object[]{header});
       }
       return result.toArray(new Object[result.size()][]);
    }

    @DataProvider(name = "wrongCreationParamData")
    public Object[][] wrongCreationParamData() {
        return new Object[][] {
                { null},
        };
    }
}
