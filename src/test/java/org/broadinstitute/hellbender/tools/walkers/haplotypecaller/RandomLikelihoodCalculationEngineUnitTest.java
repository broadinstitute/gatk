package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import com.google.common.base.Strings;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public final class RandomLikelihoodCalculationEngineUnitTest {
    @Test
    public void testComputeLikelihoods(){
        final ReadLikelihoodCalculationEngine lce = new RandomLikelihoodCalculationEngine();

        final Map<String, List<GATKRead>> perSampleReadList= new HashMap<>();
        final int n = 10 ;
        final GATKRead read1= ArtificialReadUtils.createArtificialRead(TextCigarCodec.decode(n + "M"));
        read1.setMappingQuality(60);
        final String sample1 = "sample1";
        perSampleReadList.put(sample1, Arrays.asList(read1));

        final SampleList samples = new IndexedSampleList(sample1);

        final AssemblyResultSet assemblyResultSet = new AssemblyResultSet();
        final byte[] bases = Strings.repeat("A", n + 1).getBytes();
        final Haplotype hap1 = new Haplotype(bases, true);
        hap1.setGenomeLocation(read1);
        assemblyResultSet.add(hap1);

        final byte[] basesModified= bases;
        basesModified[5] = 'C';//different bases
        final Haplotype hap2 = new Haplotype(basesModified, false);
        hap2.setGenomeLocation(read1);//use same loc
        assemblyResultSet.add(hap2);


        final ReadLikelihoods<Haplotype> likes = lce.computeReadLikelihoods(assemblyResultSet, samples, perSampleReadList);
        final LikelihoodMatrix<Haplotype> mtx = likes.sampleMatrix(0);

        Assert.assertEquals(mtx.numberOfAlleles(), 2);
        Assert.assertEquals(mtx.numberOfReads(), 1);
        final double v1 = mtx.get(0, 0);
        final double v2 = mtx.get(1, 0);

        Assert.assertTrue(v1 < 0);
        Assert.assertTrue(v2 < 0);
        lce.close();
    }
}
