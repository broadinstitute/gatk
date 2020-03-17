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
    public void testMatesKeptTogether(){
        final boolean cloud = true;
        final File out;
        final String input;
        if (cloud){
            input = "gs://broad-dsde-methods-takuto/liquid-biopsy/tmp/Jonna_Grimsby_A05_denovo_bloodbiopsy_100pct_HD78_rep1.fgbio.groupByUmi.abbrv.bam";
            out = new File("/Users/tsato/workspace/gatk/tmp/cloud.bam");
        } else {
            input = "/Users/tsato/workspace/gatk/debug/march12/Jonna_Grimsby_A05_denovo_bloodbiopsy_100pct_HD78_rep1.fgbio.groupByUmi.abbrv.bam";
            out = new File("/Users/tsato/workspace/gatk/tmp/local.bam");


        }
        // final String cloud = "gs://broad-dsde-methods/cromwell-execution-39/SpikeinNA12878/31af0e21-c493-4091-9cc6-8782b83df606/call-GroupCoffee/Jonna_Grimsby_A05_denovo_bloodbiopsy_100pct_HD78_rep1.fgbio.groupByUmi.bam";
        final ArgumentsBuilder args = new ArgumentsBuilder()
                .add("R", hg19)
                .add("I", input)
                .add("O", out.getAbsolutePath())
                .add("DS", "1.0");
        runCommandLine(args, DownsampleByDuplicateSet.class.getSimpleName());
    }

    @Test
    public void testForReal(){
        final File out = new File("/Users/tsato/Downloads/downsampled.bam");
        final String cloud = "gs://broad-dsde-methods/cromwell-execution-39/SpikeinNA12878/99a73fe8-b4fa-40d2-a866-0c95dc56162c/call-GroupCoffee/Jonna_Grimsby_A05_denovo_bloodbiopsy_100pct_HD78_rep1.fgbio.groupByUmi.bam";
        final ArgumentsBuilder args = new ArgumentsBuilder()
                .add("R", hg19)
                .add("I", cloud)
                .add("O", out.getAbsolutePath())
                .add("DS", "0.5");
        runCommandLine(args, DownsampleByDuplicateSet.class.getSimpleName());
    }

    @Test
    public void testDownsampleMaintainsCounts(){
         final File out = new File("/dsde/working/tsato/gatk/test.bam");
        final double downsampleRate = 1.0;
        final ArgumentsBuilder args = new ArgumentsBuilder()
                .add("R", hg19)
                .add("I", medium)
                .add("DS", Double.toString(downsampleRate))
                .add("O", out.getAbsolutePath());
        runCommandLine(args, DownsampleByDuplicateSet.class.getSimpleName());

        final ReadsDataSource originalBam = new ReadsDataSource(Paths.get(medium));
        final Map<String, MutableInt> originalMoleculeCounts = molecularIDsAndCounts(originalBam);

        final ReadsDataSource downsampledBam = new ReadsDataSource(Paths.get(out.getAbsolutePath()));
        final Map<String, MutableInt> downsampledMoleculeCounts = molecularIDsAndCounts(downsampledBam);

        for (Map.Entry<String, MutableInt> originalIDAndCount : originalMoleculeCounts.entrySet()){
            final String originalID = originalIDAndCount.getKey();
            final int originalCount = originalIDAndCount.getValue().intValue();
            Assert.assertTrue(originalCount == downsampledMoleculeCounts.get(originalID).intValue());
        }
    }
    @Test
    public void testDownsampleFraction(){
        // final File out = createTempFile("downsampled", "bam");
        final File out = new File("/dsde/working/tsato/gatk/test.bam");
        for (double downsampleRate : Arrays.asList(0.1, 0.3, 0.5)){
            final ArgumentsBuilder args = new ArgumentsBuilder()
                    .add("R", hg19)
                    .add("I", medium)
                    .add("DS", Double.toString(downsampleRate))
                    .add("O", out.getAbsolutePath());
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

    private Map<String, MutableInt> molecularIDsAndCounts(final ReadsDataSource readsDataSource){
        final Map<String, MutableInt> map = new TreeMap<>();
        final Iterator<GATKRead> iterator = readsDataSource.iterator();
        while (iterator.hasNext()){
            final GATKRead read = iterator.next();
            final String molecularID = read.getAttributeAsString("MI"); // Note we are duplex aware: 12/A different from 12/B
            if (map.containsKey(molecularID)){
                map.get(molecularID).increment();
            } else {
                map.put(molecularID, new MutableInt(0));
            }
        }

        return map;
    }


}