package org.broadinstitute.hellbender.utils.fermi;

import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.IOException;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

public class FermiLiteAssemblerTest extends BaseTest {
    @Test
    public void testFermiLiteAssembly() {
        final List<GATKRead> readsList = new ArrayList<>(1200);
        final ReadsDataSource reads = new ReadsDataSource(Paths.get(getToolTestDataDir()+"x.bam"));
        for ( GATKRead read : reads ) {
            readsList.add(read);
        }
        final Assembly assembly = new FermiLiteAssembler().createAssembly(readsList);
        Assert.assertEquals(assembly.getNContigs(), 1);
        try( final ReferenceSequenceFile refSeqFile =
                     ReferenceSequenceFileFactory.getReferenceSequenceFile(Paths.get(getToolTestDataDir()+"x.fa")) ) {
            final ReferenceSequence seq = refSeqFile.nextSequence();
            Assert.assertEquals(assembly.getContig(0).getSequence(), seq.getBases());
        }
        catch ( IOException ioe ) {
            throw new GATKException("can't read expected contig sequence", ioe);
        }
    }

}
