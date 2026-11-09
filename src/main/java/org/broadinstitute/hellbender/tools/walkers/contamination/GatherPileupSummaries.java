package org.broadinstitute.hellbender.tools.walkers.contamination;


import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.utils.reference.ReferenceUtils;
import picard.cmdline.programgroups.DiagnosticsAndQCProgramGroup;

import java.io.File;
import java.util.ArrayList;
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
        // Each file's first site, if it has any, is enough to order the files, under the assumption that
        //   1. PileupSummaries are already sorted within each file
        //   2. Files do not overlap
        final List<Pair<File, PileupSummary>> filesAndFirstSites = new ArrayList<>();
        for (final File file : input) {
            PileupSummary.readFirstSite(file).ifPresent(firstSite -> filesAndFirstSites.add(Pair.of(file, firstSite)));
        }
        if (filesAndFirstSites.size() < input.size()){
            logger.info(String.format("Removed %d empty samples", input.size() - filesAndFirstSites.size()));
        }

        final List<File> nonEmptyFiles = filesAndFirstSites.stream()
                .sorted(Comparator.comparing(Pair::getRight, new PileupSummary.PileupSummaryComparator(sequenceDictionary)))
                .map(Pair::getLeft)
                .collect(Collectors.toList());
        PileupSummary.writeToFile(nonEmptyFiles, output);
        return String.format("Successfully merged %d samples", nonEmptyFiles.size());
    }
}
