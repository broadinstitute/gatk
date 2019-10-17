package org.broadinstitute.hellbender.tools.walkers.contamination;


import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.utils.reference.ReferenceUtils;
import picard.cmdline.programgroups.DiagnosticsAndQCProgramGroup;

import java.io.File;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

@CommandLineProgramProperties(
        summary="Combine output files from GetPileupSummary in the order defined by a sequence dictionary",
        oneLineSummary = "Combine output files from GetPileupSummary in the order defined by a sequence dictionary",
        programGroup = DiagnosticsAndQCProgramGroup.class
)
public class GatherPileupSummaries extends CommandLineProgram {
    @Argument(fullName = StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME, doc = "sequence dictionary file")
    final File sequenceDictionaryFile = null;

    @Argument(fullName = StandardArgumentDefinitions.INPUT_SHORT_NAME, doc = "an output of PileupSummaryTable")
    final List<File> input = null;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "output")
    final File output = null;

    SAMSequenceDictionary sequenceDictionary = null;

    @Override
    protected void onStartup(){
        sequenceDictionary = ReferenceUtils.loadFastaDictionary(sequenceDictionaryFile);
    }

    @Override
    protected Object doWork() {
        final List<File> nonEmptyFiles = removeEmptyFiles(input);
        Collections.sort(nonEmptyFiles, new PileupSummaryFileComparator(sequenceDictionary));
        PileupSummary.writeToFile(nonEmptyFiles, output);
        return String.format("Successfully merged %d samples", nonEmptyFiles.size());
    }

    private List<File> removeEmptyFiles(final List<File> list){
        final List<File> nonEmptyList = list.stream().filter(f -> PileupSummary.readFromFile(f).getRight().size() > 0)
                .collect(Collectors.toList());
        if (nonEmptyList.size() < list.size()){
            logger.info(String.format("Removed %d empty samples", list.size() - nonEmptyList.size()));
        }
        return nonEmptyList;
    }

    /**
     * Compare two PilupSummary files under the assumption that
     *   1. PileupSummaries are already sorted within each file
     *   2. Files do not overlap
     */
    private class PileupSummaryFileComparator implements Comparator<File> {
        final SAMSequenceDictionary sequenceDictionary;

        private PileupSummaryFileComparator(final SAMSequenceDictionary sequenceDictionary){
            this.sequenceDictionary = sequenceDictionary;
        }

        @Override
        public int compare(File file1, File file2) {
            final PileupSummary ps1 = PileupSummary.readFromFile(file1).getRight().get(0);
            final PileupSummary ps2 = PileupSummary.readFromFile(file2).getRight().get(0);

            return new PileupSummary.PileupSummaryComparator(sequenceDictionary).compare(ps1, ps2);
        }
    }
}
