package org.broadinstitute.hellbender.engine.dataflow;

import com.google.api.services.genomics.model.Read;
import com.google.appengine.repackaged.com.google.common.collect.Lists;
import com.google.cloud.genomics.dataflow.readers.bam.ReadConverter;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.read.ArtificialSAMUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.List;
import java.util.stream.Collectors;

public class SAMSerializableFunctionTest extends CommandLineProgramTest{

    @Test
    public void testApply(){
        final SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader();
        final List<Integer> lengths = Lists.newArrayList(1, 2, 100, 1000);
        List<Read> reads = lengths.stream()
                .map(ArtificialSAMUtils::createRandomRead)
                .map(ReadConverter::makeRead)
                .collect(Collectors.toList());


        SAMSerializableFunction<Integer> f = new SAMSerializableFunction<>(header.toString() ,SAMRecord::getReadLength);

        final List<Integer> recoveredLengths = reads.stream()
                .map(f::apply)
                .collect(Collectors.toList());

        Assert.assertEquals(lengths, recoveredLengths);
    }
}