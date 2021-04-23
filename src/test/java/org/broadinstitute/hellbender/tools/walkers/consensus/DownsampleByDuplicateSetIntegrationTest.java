package org.broadinstitute.hellbender.tools.walkers.consensus;

import htsjdk.samtools.SAMTag;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.engine.ReadsPathDataSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.nio.file.Paths;
import java.util.*;

public class DownsampleByDuplicateSetIntegrationTest extends CommandLineProgramTest {
    private static final String downsampleByDuplicateSetTestDir = publicTestDir + "org/broadinstitute/hellbender/tools/downsampleByDuplicateSet/";
    private static final String NA12878_GROUPED = downsampleByDuplicateSetTestDir + "NA12878.grouped.bam";

    @Test
    public void testMatesAreTogether() {
        final File out = createTempFile("downsampled", "bam");
        final ArgumentsBuilder args = new ArgumentsBuilder()
                .add("I", NA12878_GROUPED)
                .add("O", out.getAbsolutePath())
                .add(DownsampleByDuplicateSet.FRACTION_TO_KEEP_NAME, "1.0");
        runCommandLine(args, DownsampleByDuplicateSet.class.getSimpleName());

        try (final ReadsDataSource readsDataSource = new ReadsPathDataSource(Paths.get(out.getAbsolutePath()))) {
            final Iterator<GATKRead> iterator = readsDataSource.iterator();
            while (iterator.hasNext()) {
                // Make sure that the read and its mate are next to each other in the file
                final GATKRead read1 = iterator.next();
                final GATKRead read2 = iterator.next();
                Assert.assertEquals(read1.getName(), read2.getName());
            }
        } catch (Exception e) {
            throw new UserException("Encountered an IOException writing to " + out, e);
        }
    }

    /**
     * Test that the downsampling rate corresponds to the reduction in the number of duplicates in the output
     * file up to sampling noise.
     */
    @Test
    public void testDownsampleFraction(){
        final File out = createTempFile("downasampled", "bam");
        for (double downsampleRate : Arrays.asList(0.0, 0.1, 0.3, 0.5, 1.0)){
            final ArgumentsBuilder args = new ArgumentsBuilder()
                    .add("I", NA12878_GROUPED)
                    .add(DownsampleByDuplicateSet.FRACTION_TO_KEEP_NAME, Double.toString(downsampleRate))
                    .add("O", out.getAbsolutePath());
            runCommandLine(args, DownsampleByDuplicateSet.class.getSimpleName());

            try(final ReadsDataSource originalBam = new ReadsPathDataSource(Paths.get(NA12878_GROUPED));
                final ReadsDataSource downsampledBam = new ReadsPathDataSource(Paths.get(out.getAbsolutePath()))){

                final int originalMoleculeCount = countDuplicateSets(originalBam);
                final int downsampledMoleculeCount = countDuplicateSets(downsampledBam);

                if (downsampleRate == 0.0){ // Keep none of the reads
                    Assert.assertEquals(downsampledMoleculeCount, 0.0);
                    continue;
                }

                if (downsampleRate == 1.0){ // Keep all of the reads
                    Assert.assertEquals(downsampledMoleculeCount, originalMoleculeCount);
                    continue;
                }

                // In addition to the stochastic noise, there's an additional source of deviation from the mean due to
                // the fact that we are downsampling by molecule, whose family size may vary and therefore is not the same thing as
                // dropping 5% of the reads, for example.
                final double noise = originalMoleculeCount*0.05; // allow for a 5% deviation from the expected
                final double deviationFromExpected = Math.abs(downsampleRate * originalMoleculeCount - downsampledMoleculeCount);
                Assert.assertTrue(deviationFromExpected < noise);
            } catch (Exception e) {
                throw new UserException("Encountered an IO error", e);
            }
        }
    }

    private int countDuplicateSets(final ReadsDataSource readsDataSource){
        int count = 0;
        String currentMoleculeId = ""; // Note we are duplex aware: 12/A different from 12/B

        final Iterator<GATKRead> iterator = readsDataSource.iterator();
        while (iterator.hasNext()){
            final GATKRead read = iterator.next();
            final String moleculeID = read.getAttributeAsString(SAMTag.MI.name());
            if (!moleculeID.equals(currentMoleculeId)){
                count++;
                currentMoleculeId = moleculeID;
            }
        }

        return count;
    }
}