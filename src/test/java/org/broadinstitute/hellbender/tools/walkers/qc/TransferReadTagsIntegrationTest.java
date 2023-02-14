package org.broadinstitute.hellbender.tools.walkers.qc;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.util.RuntimeIOException;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.engine.ReadsPathDataSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.tools.walkers.mutect.M2TestingUtils;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.read.SAMFileGATKReadWriter;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class TransferReadTagsIntegrationTest extends CommandLineProgramTest {

    public File getSamFile(final String prefix, final boolean aligned, final SAMFileHeader.SortOrder sortOrder, final Map<String, String> readNamesAndUmis) {
        final int readLength = 30;
        final byte[] bases = new byte[readLength];
        final byte[] quals = new byte[readLength];
        Arrays.fill(bases, (byte)'A');
        Arrays.fill(quals, (byte)30);

        final File samFile = createTempFile(prefix, ".bam");
        final SAMFileHeader samHeader = M2TestingUtils.createSamHeader("sample");
        if (sortOrder != null) {
            samHeader.setSortOrder(sortOrder);
        }

        try (final SAMFileGATKReadWriter samWriter = new SAMFileGATKReadWriter(ReadUtils.createCommonSAMWriter(samFile, null, samHeader, true, false, false))) {
            for (final Map.Entry<String, String> readNameAndUMI : readNamesAndUmis.entrySet()) {
                final GATKRead read = aligned? ArtificialReadUtils.createArtificialRead(samHeader, bases, quals, "30M") :
                        ArtificialReadUtils.createArtificialUnmappedRead(samHeader, bases, quals);
                read.setName(readNameAndUMI.getKey());
                final String umi = readNameAndUMI.getValue();
                if (umi != null) {
                    read.setAttribute("RX", umi);
                }
                samWriter.addRead(read);
            }
        } catch (RuntimeIOException e){
            throw new UserException("Failed to open SAM writers", e);
        }

        return samFile;
    }


    @DataProvider(name = "testSAMPairs")
    public Object[][] getTestSAMPairs(){
        final List<String> umis = Arrays.asList("TCC-ATG", "CTA-GGA", "CCT-AGC", "GCA-AAA");
        final LinkedHashMap<String, String> alignedUMIsMap = new LinkedHashMap<>();
        final LinkedHashMap<String, String> unmappedUMIsMap = new LinkedHashMap<>();
        final LinkedHashMap<String, String> alignedWithExistingTagsUMIsMap = new LinkedHashMap<>();
        final LinkedHashMap<String, String> expectedUMIsMap = new LinkedHashMap<>();

        for (int i = 0; i < umis.size(); i++) {
            if (i != 1) {
                alignedUMIsMap.put("read" + i, null);
                alignedWithExistingTagsUMIsMap.put("read"+ i, "ZZZ");
                expectedUMIsMap.put("read"+ i, umis.get(i));
            }
            unmappedUMIsMap.put("read" + i, umis.get(i));

        }
        final File alignedSAMFile = getSamFile("aligned", true, SAMFileHeader.SortOrder.queryname, alignedUMIsMap);
        final File alignedSAMWithExistingTagFile = getSamFile("aligned_with_existing_tag", true, SAMFileHeader.SortOrder.queryname, alignedWithExistingTagsUMIsMap);
        final File unmappedSAMFile = getSamFile("unmapped", false, null, unmappedUMIsMap);
        final File emptySAMFile = getSamFile("empty", false, SAMFileHeader.SortOrder.queryname, Collections.emptyMap());

        return new Object[][]{{alignedSAMFile, unmappedSAMFile, expectedUMIsMap },
                {alignedSAMWithExistingTagFile, unmappedSAMFile, expectedUMIsMap},
                {emptySAMFile, emptySAMFile, Collections.emptyMap()},
                {emptySAMFile, unmappedSAMFile, Collections.emptyMap()}
        };
    }

    @Test(dataProvider = "testSAMPairs")
    public void test(final File alignedSAMFile, final File unmappedSAMFile, final Map<String, String> expectedReadNamesAndUmis) {
        final File outputFile = createTempFile("output", ".bam");
        final ArgumentsBuilder args = new ArgumentsBuilder()
                .add("I", alignedSAMFile.getAbsolutePath())
                .add("unmapped-sam", unmappedSAMFile.getAbsolutePath())
                .add("read-tags", "RX")
                .add("O", outputFile.getAbsolutePath());
        runCommandLine(args);

        final ReadsPathDataSource outputReadSource = new ReadsPathDataSource(outputFile.toPath());

        final Set<String> observedReadName = new HashSet<>();
        for (GATKRead outputRead : outputReadSource) {
            Assert.assertTrue(expectedReadNamesAndUmis.containsKey(outputRead.getName()));
            observedReadName.add(outputRead.getName());
            Assert.assertEquals(outputRead.getAttributeAsString("RX"), expectedReadNamesAndUmis.get(outputRead.getName()));
        }
        Assert.assertTrue(observedReadName.containsAll(expectedReadNamesAndUmis.keySet()));
    }

    @Test(expectedExceptions = UserException.class, expectedExceptionsMessageRegExp = "Unmapped sam iterator is empty and aligned sam iterator is not.")
    public void testEmptyUnmapped() {
        final File alignedSAMFile = getSamFile("aligned", true, SAMFileHeader.SortOrder.queryname, Collections.singletonMap("read0", null));
        final File emptySAMFile = getSamFile("empty", false, SAMFileHeader.SortOrder.queryname, Collections.emptyMap());
        final File outputFile = createTempFile("output", ".bam");
        final ArgumentsBuilder args = new ArgumentsBuilder()
                .add("I", alignedSAMFile.getAbsolutePath())
                .add("unmapped-sam", emptySAMFile.getAbsolutePath())
                .add("read-tags", "RX")
                .add("O", outputFile.getAbsolutePath());
        runCommandLine(args);
    }
}