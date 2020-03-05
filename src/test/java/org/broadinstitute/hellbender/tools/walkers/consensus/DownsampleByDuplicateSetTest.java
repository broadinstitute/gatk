package org.broadinstitute.hellbender.tools.walkers.consensus;

import htsjdk.samtools.SamReader;
import org.apache.commons.lang3.mutable.Mutable;
import org.apache.commons.lang3.mutable.MutableInt;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.iterators.SamReaderQueryingIterator;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.runtime.ProcessController;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.Reader;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;


import static org.testng.Assert.*;

public class DownsampleByDuplicateSetTest extends CommandLineProgramTest {
    private final String hg19 = "/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta";
    private final String medium = "/dsde/working/tsato/consensus/tp53/test/bams/medium/Jonna_Grimsby_A04_denovo_bloodbiopsy_1pct_rep1.tp53.CTG-TTC.grouped.bam";
    private final String countScript = "/dsde/working/tsato/consensus/tp53/test/bams/count_MIs.sh";

    @Test
    public void testForReal(){
        final File out = new File("/Users/tsato/Downloads/downsampled.bam");
        final String cloud = "gs://broad-dsde-methods/cromwell-execution-39/SpikeinNA12878/99a73fe8-b4fa-40d2-a866-0c95dc56162c/call-GroupCoffee/Jonna_Grimsby_A05_denovo_bloodbiopsy_100pct_HD78_rep1.fgbio.groupByUmi.bam";
        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addArgument("R", hg19)
                .addArgument("I", cloud)
                .addArgument("O", out.getAbsolutePath())
                .addArgument("DS", "0.5");
        runCommandLine(args, DownsampleByDuplicateSet.class.getSimpleName());
    }

    @Test
    public void testDownsampleFraction(){
        // final File out = createTempFile("downsampled", "bam");
        final File out = new File("/dsde/working/tsato/gatk/test.bam");
        for (double downsampleRate : Arrays.asList(0.1, 0.3, 0.5)){
            final ArgumentsBuilder args = new ArgumentsBuilder()
                    .addArgument("R", hg19)
                    .addArgument("I", medium)
                    .addArgument("DS", Double.toString(downsampleRate))
                    .addArgument("O", out.getAbsolutePath());
            runCommandLine(args, DownsampleByDuplicateSet.class.getSimpleName());

            final ReadsDataSource originalBam = new ReadsDataSource(Paths.get(medium));
            final int originalMoleculeCount = countDuplicateSets(originalBam);

            final ReadsDataSource downsampledBam = new ReadsDataSource(Paths.get(out.getAbsolutePath()));
            final int downsampledMoleculeCount = countDuplicateSets(downsampledBam);

            final int noise = 10; // TODO: Can I do better?
            Assert.assertTrue(Math.abs(downsampleRate * originalMoleculeCount - downsampledMoleculeCount) < noise);
        }

    }

    private int countDuplicateSets(final ReadsDataSource readsDataSource){
        int count = 0;
        String currentMolecularId = ""; // Note we are duplex aware: 12/A different from 12/B
        final Iterator<GATKRead> iterator = readsDataSource.iterator();
        while (iterator.hasNext()){
            final GATKRead read = iterator.next();
            final String molecularID = read.getAttributeAsString("MI");
            if (!molecularID.equals(currentMolecularId)){
                count++;
                currentMolecularId = molecularID;
            }
        }

        return count;
    }


}