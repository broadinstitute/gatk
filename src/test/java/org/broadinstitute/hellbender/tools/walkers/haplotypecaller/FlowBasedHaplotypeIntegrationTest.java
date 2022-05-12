package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.FlowBasedAlignmentArgumentCollection;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.genomicsdb.importer.model.ChromosomeInterval;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import org.broadinstitute.hellbender.utils.haplotype.FlowBasedHaplotype;
import org.broadinstitute.hellbender.utils.read.FlowBasedRead;
import org.broadinstitute.hellbender.tools.FlowBasedArgumentCollection;

import java.io.File;
import java.util.ArrayList;

public class FlowBasedHaplotypeIntegrationTest extends GATKBaseTest {
    private final Logger logger = LogManager.getLogger(this.getClass());

    @DataProvider(name = "haplotypeMatchingDataProvider")
    public Object[][] gethaplotypeMatchingTestData() {
        final Object[][]        testData = {

                { publicTestDir + "/large/FlowBasedHaplotype_HC_flow_chr9.part.bam", null, 100 }
        };

        return testData;

    }

    @Test(dataProvider = "haplotypeMatchingDataProvider")
    public void testHaplotypeMatching(final String readsFile, final String CRValue, final int limit) {

        // bypass until we have a real file
        if ( readsFile == null )
            return;

        // read reads
        final SamReader reader = SamReaderFactory.makeDefault().open(new File(readsFile));
        final ArrayList<GATKRead> reads = new ArrayList<>();
        final ArrayList<Haplotype> haplotypes = new ArrayList<>();
        int count_reads = 0;
        int count_haplotypes = 0;
        for ( final SAMRecord rec : reader ) {
            final SAMRecordToGATKReadAdapter tmp = new SAMRecordToGATKReadAdapter(rec);
            if (!tmp.getAttributeAsString("RG").startsWith("ArtificialHaplotype")) {
                if ((CRValue == null ) || tmp.getAttributeAsString(("CR")).equals(CRValue)) {
                    if ( count_reads < limit ) {
                        reads.add(tmp);
                        count_reads++;
                    }
                }
            } else {
                if ((CRValue == null) || tmp.getAttributeAsString(("CR")).equals(CRValue) ) {
                    ChromosomeInterval gl = new ChromosomeInterval(tmp.getContig(), tmp.getStart(), tmp.getEnd());
                    if ( count_haplotypes < limit ) {
                        Haplotype hap = new Haplotype(tmp.getBases(), gl);
                        hap.setCigar(tmp.getCigar());
                        haplotypes.add(hap);
                        count_haplotypes++;
                    }
                }
            }

        }

        logger.debug(String.format("%d reads %d haplotypes", count_reads, count_haplotypes));
        final ArrayList<FlowBasedRead> fbrs = new ArrayList<>();
        final FlowBasedArgumentCollection fbargs = new FlowBasedArgumentCollection();
        for (final GATKRead r : reads) {
            fbrs.add(new FlowBasedRead(r, "TACG", 8, fbargs));
        }
        for (FlowBasedRead fbr : fbrs ) {
            fbr.applyAlignment();
        }

        final ArrayList<FlowBasedHaplotype> fbhs = new ArrayList<>();
        for (final Haplotype hap : haplotypes) {
            fbhs.add(new FlowBasedHaplotype(hap, "TACG"));
        }

        final FlowBasedAlignmentLikelihoodEngine fbe = new FlowBasedAlignmentLikelihoodEngine(new FlowBasedAlignmentArgumentCollection(), -5, 0.02, false, PairHMMLikelihoodCalculationEngine.DEFAULT_DYNAMIC_DISQUALIFICATION_SCALE_FACTOR);

        final AlleleLikelihoods<GATKRead, Haplotype> haplotypeReadLikelihoods =
                FlowBasedAlignmentLikelihoodEngineTestUtils.computeReadLikelihoods(haplotypes, reads, true, reader.getFileHeader(), fbe);
        logger.debug("haplotypeReadLikelihoods: " + haplotypeReadLikelihoods);
    }
}
