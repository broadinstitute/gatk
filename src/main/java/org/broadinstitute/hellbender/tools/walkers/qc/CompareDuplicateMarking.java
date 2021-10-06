package org.broadinstitute.hellbender.tools.walkers.qc;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

// Define an interface with methods match()
public class CompareDuplicateMarking implements QuerynameSetComparison {
    long genomeDupTranscrNot = 0;
    long bothDup = 0;
    long genomeNotTranscrDup = 0;
    long neitherDup;
    long genomeDupNotFoundInTranscr = 0;
    long genomeAbsentTranscriptDup = 0;
    long notFoundInTranscr = 0;
    long notFoundInGenome = 0;
    long numReadPairsGenome = 0;
    long numReadPairsTranscr = 0;

    public CompareDuplicateMarking(){}

    /**
     *
     * The queryname set was found in the genome, but not in the transcriptome.
     *
     * @param input1ReadPair
     */
    @Override
    public void processInput1(final ReadPair input1ReadPair){
        if (input1ReadPair.isDuplicateMarked()){
            genomeDupNotFoundInTranscr += 1;
        } else {
            notFoundInTranscr += 1;
        }
    }

    @Override
    public void processInput2(final ReadPair input2ReadPair){
        if (input2ReadPair.isDuplicateMarked()){
            genomeAbsentTranscriptDup += 1;
        } else {
            notFoundInGenome += 1;
        }
    }

    @Override
    public void processMatchingQuerynameSets(final ReadPair input1ReadPair, final ReadPair input2ReadPair){
        // This code chunk is a good candidate for Scala's match, which would improve readability
        if (input1ReadPair.isDuplicateMarked()){
            if (input2ReadPair.isDuplicateMarked()){
                bothDup += 1;
            } else {
                genomeDupTranscrNot += 1;
            }
        } else {
            if (input2ReadPair.isDuplicateMarked()){
                genomeNotTranscrDup += 1;
            } else {
                neitherDup += 1;
            }
        }
    }

    @Override
    public void writeSummary(final File outputTable, final SAMFileHeader header){
        try (PrintWriter pw = new PrintWriter(outputTable)){
            final String sample = header.getReadGroups().get(0).getSample();
            pw.println("#sample=" + sample);
            pw.println("genomeDupTranscrNot," + genomeDupTranscrNot);
            pw.println("bothDup," + bothDup);
            pw.println("genomeNotTranscrDup," + genomeNotTranscrDup);
            pw.println("neitherDup," + neitherDup);
            pw.println("genomeDupNotFoundInTranscr," + genomeDupNotFoundInTranscr);
            pw.println("notFoundInTranscr," + notFoundInTranscr);
            pw.println("notFoundInGenome," + notFoundInGenome);
            pw.println("numReadPairsGenome," + numReadPairsGenome);
            pw.println("numReadPairsTranscr," + numReadPairsTranscr);
        } catch (IOException e){
            throw new UserException("Could not write to " + outputTable, e);
        }
    }
}
