package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.samtools.*;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.FlowBasedAlignmentArgumentCollection;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.nio.file.FileSystems;
import java.nio.file.Path;
import java.util.*;
import java.util.stream.Collectors;

public class FlowBasedAlignmentLikelihoodEngineUnitTest extends GATKBaseTest {

    final static String    testResourceDir = publicTestDir + "org/broadinstitute/hellbender/utils/read/flow/reads/";
    final static String    inputDir = testResourceDir + "/input/";
    final Random           random = new Random(0);

    @DataProvider(name = "testData")
    public Object[][] getTestData() {

        // create engine, reads
        final Path                          inputFile = FileSystems.getDefault().getPath(inputDir, "sample.bam");
        final SamReader                     reader = SamReaderFactory.makeDefault().open(new File(inputFile.toString()));
        final List<GATKRead>                reads = getTestReads(reader);

        // build test data
        final Object[][]                    testData = new Object[reads.size()][];
        for ( int i = 0 ; i < reads.size() ; i++ ) {
            testData[i] = new Object[2];
            testData[i][0] = reads.get(i);
            testData[i][1] = reader.getFileHeader();
        }

        return testData;
    }

    @Test(dataProvider = "testData")
    void testComputeReadLikelihoodsWithRandomlyModifiedHmers(final GATKRead read, final SAMFileHeader fileHeader) {

        // create engine
        final FlowBasedAlignmentLikelihoodEngine engine = getTestAlignmentEngine();

        // for each read, create extrapolated reads, and haplotypes. then invoke engine
        final List<GATKRead>  extrapolatedReads = getTestExtrapolatedReads(read);
        final List<Haplotype> haplotypes = getTestHaplotypes(extrapolatedReads);

        // run the engine
        final AlleleLikelihoods<GATKRead, Haplotype> result = FlowBasedAlignmentLikelihoodEngineTestUtils.computeReadLikelihoods(haplotypes, extrapolatedReads, false, fileHeader, engine);
        Assert.assertNotNull(result);

        // check results (sanity)
        LikelihoodMatrix<GATKRead, Haplotype> matrix = result.sampleMatrix(0);
        Assert.assertNotNull(matrix);
        Assert.assertEquals(matrix.evidence().size(), extrapolatedReads.size());
        Assert.assertEquals(matrix.alleles().size(), haplotypes.size());

        // check results, first ref, second unmodified, rest modified
        double refScore = matrix.get(0, 0);
        for (int alleleIndex = 1; alleleIndex < matrix.numberOfAlleles(); alleleIndex++) {
            double score = matrix.get(alleleIndex, 0);
            if (alleleIndex == 1) {
                Assert.assertEquals(refScore, score);
            } else {
                Assert.assertTrue(refScore > score);
            }
        }
    }

    // create a new engine (to be used in testing code)
    private FlowBasedAlignmentLikelihoodEngine getTestAlignmentEngine() {

        // create (and possibly initialise) arguments
        FlowBasedAlignmentArgumentCollection args = new FlowBasedAlignmentArgumentCollection();

        // create engine
        return new FlowBasedAlignmentLikelihoodEngine(args,-5, 0.02, false, PairHMMLikelihoodCalculationEngine.DEFAULT_DYNAMIC_DISQUALIFICATION_SCALE_FACTOR);
    }

    // get test reads from a test file
    private List<GATKRead> getTestReads(SamReader reader) {

        return reader
                .iterator()
                .stream().map(samRecord -> new SAMRecordToGATKReadAdapter(samRecord))
                .collect(Collectors.toList());
    }

    // for now, we simply use the read itself (no extrapolated reads added)
    private List<GATKRead> getTestExtrapolatedReads(GATKRead read) {
        return Collections.singletonList(read);
    }

    // get test haplotypes from a set of reads
    private List<Haplotype> getTestHaplotypes(final List<GATKRead> reads) {

        final List<Haplotype>       haplotypes = new LinkedList<>();

        // this code builds haplotypes from the first read
        GATKRead                    read = reads.get(0);
        for ( int i = 0 ; i < 20 ; i++ ) {

            // first haplotype is reference, others are messed up randomly, except for first
            byte[]      bases = read.getBases();
            if ( i > 1 ) {
                final int     ofs = 10 + i;      // this ensures that no two haplotypes will be alike.
                                                 // also skipping initial hmers since they might be problematic
                final int     mod = random.nextInt(3);  // 0-SNP, 1-DEL, 2-INS
                if ( mod == 0 ) {
                    bases[ofs] = modifyBase(bases[ofs]);
                } else if ( mod == 1 ) {
                    final byte[]    delBases = new byte[bases.length - 1];
                    System.arraycopy(bases, 0, delBases, 0, ofs);
                    System.arraycopy(bases, ofs + 1, delBases, ofs, bases.length - ofs - 1);
                    bases = delBases;
                } else if ( mod == 2 ) {
                    final byte[]    insBases = new byte[bases.length + 1];
                    System.arraycopy(bases, 0, insBases, 0, ofs + 1);
                    insBases[ofs + 1] = insBases[ofs];
                    System.arraycopy(bases, ofs, insBases, ofs + 2, bases.length - ofs - 1);
                    bases = insBases;
                }
            }

            final Haplotype haplotype = new Haplotype(bases, i == 0);
            haplotype.setCigar(new Cigar(Collections.singletonList(new CigarElement(haplotype.length(), CigarOperator.M))));
            haplotype.setUniquenessValue(i);
            haplotype.setGenomeLocation(read);

            haplotypes.add(haplotype);
        }

        return haplotypes;
    }

    private byte modifyBase(byte b) {
        switch (b) {
            case 'A' : return 'C';
            case 'C' : return 'T';
            case 'T' : return 'G';
            case 'G' : return 'A';
            default: return 'N';
        }
    }
}
