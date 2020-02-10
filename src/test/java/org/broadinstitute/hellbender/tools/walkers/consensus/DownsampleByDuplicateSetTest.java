package org.broadinstitute.hellbender.tools.walkers.consensus;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.runtime.ProcessController;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.Reader;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.stream.Collectors;
import java.util.stream.Stream;


import static org.testng.Assert.*;

public class DownsampleByDuplicateSetTest extends CommandLineProgramTest {
    private final String hg19 = "/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta";
    private final String medium = "/dsde/working/tsato/consensus/tp53/test/bams/medium/Jonna_Grimsby_A04_denovo_bloodbiopsy_1pct_rep1.tp53.CTG-TTC.grouped.bam";
    private final String countScript = "/dsde/working/tsato/consensus/tp53/test/bams/count_MIs.sh";

    @Test
    public void testTemp() {
        final File out = new File("/dsde/working/tsato/consensus/tp53/test/tmp/test.bam");
        final String cloud = "gs://broad-dsde-methods/cromwell-execution-39/SpikeinNA12878/d6425c0d-4282-4a7b-bee9-f59134e500aa/call-GroupCoffee/Jonna_Grimsby_A05_denovo_bloodbiopsy_100pct_HD78_rep1.fgbio.groupByUmi.bam";

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addArgument("R", hg19)
                .addArgument("I", cloud)
                .addArgument("DS", "0.5")
                .addArgument("O", out.getAbsolutePath());
        runCommandLine(args, DownsampleByDuplicateSet.class.getSimpleName());
    }

    @Test
    public void test(){
        final File out = new File("/dsde/working/tsato/consensus/tp53/test/tmp/test.bam");
        final File countTableBefore = createTempFile("before", "csv");
        runProcess(new ProcessController(), new String[]{ countScript, medium, countTableBefore.getAbsolutePath()});
        final String cloud = "gs://broad-dsde-methods/cromwell-execution-39/SpikeinNA12878/d6425c0d-4282-4a7b-bee9-f59134e500aa/call-GroupCoffee/Jonna_Grimsby_A05_denovo_bloodbiopsy_100pct_HD78_rep1.fgbio.groupByUmi.bam";

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addArgument("R", hg19)
                .addArgument("I", medium)
                .addArgument("DS", "0.5")
                .addArgument("O", out.getAbsolutePath());
        runCommandLine(args, DownsampleByDuplicateSet.class.getSimpleName());

        final File countTableAfter = createTempFile("after", "csv");
        runProcess(new ProcessController(), new String[]{ countScript, out.getAbsolutePath(), countTableAfter.getAbsolutePath()});


        final List<Integer> idCountBefore = new ArrayList<>();
        final Map<Integer, Integer> idCountMap = new TreeMap<>();

        // IO code so unwieldly!
        try (Stream<String> stream = Files.lines(Paths.get(countTableBefore.getAbsolutePath()))) {
            stream.forEach(line -> {
                final String[] str = line.split(",");
                idCountBefore.add(Integer.parseInt(str[1]));
            });
        } catch (IOException e) {
            throw new UserException("", e);
        }

        try (Stream<String> stream = Files.lines(Paths.get(countTableAfter.getAbsolutePath()))) {
            stream.forEach(line -> {
                final String[] str = line.split(",");
                idCountMap.put(Integer.parseInt(str[0]), Integer.parseInt(str[1]));
            });
        } catch (IOException e) {
            throw new UserException("", e);
        }

        for (Map.Entry<Integer, Integer> entry : idCountMap.entrySet()){
            final int id = entry.getKey();
            final int count = entry.getValue();
            if (count > 0){
                Assert.assertEquals(count, (int) idCountBefore.get(id)); // why do I need to cast here?
            }
        }
    }

}