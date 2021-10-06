package org.broadinstitute.hellbender.tools.walkers.qc;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.List;
import java.util.StringJoiner;

/**
 * Recommended preprocessing.
 * 1. PrintReads -L Ribosome the (hg19-Twist) bam. Call the output (hg19-Twist-Ribo-Pre)
 * 2. samtools view -N (hg19-Twist-Ribo-Pre) (hg38-Twist). (Filters the hg38-Twist bam to the queryname sets of interest)
 *    samtools view -N (hg19-Twist-Ribo-Pre) (hg19-Twist). (This recovers the mates)
 * 3. samtools sort (hg19-Twist-Ribosome) by queryname.
 *    samtools sort (hg38-Twist-Ribosome) by queryname.
 * 4. Feed the resulting file to this tool.
 *
 */

public class CompareRibosomalReads implements QuerynameSetComparison {
    // These will be added after getting basic functionaly in place.
    private File hg19RibosomalIntervals = new File("");
    private List<SimpleInterval> hg19RibosomalInterval;

    PrintWriter tableWriter;

    // These should always be 0.
    int onlyInHg19 = 0;
    int onlyInHg38 = 0;

    public CompareRibosomalReads(File outputFile) throws FileNotFoundException {
        tableWriter = new PrintWriter(outputFile);
        tableWriter.println("hg19Contig\thg19Start\thg19End\thg38Contig\thg38Start\thg38End");
    }

    @Override
    public void processInput1(ReadPair input1ReadPair) {
        onlyInHg19++;
    }

    @Override
    public void processInput2(ReadPair input2ReadPair) {
        onlyInHg38++;
    }

    @Override
    public void processMatchingQuerynameSets(ReadPair input1ReadPair, ReadPair input2ReadPair) {
        final GATKRead hg19FirstOfPair = input1ReadPair.getFirstOfPair();
        // There has to be a way to get a simple interval of the insert (or read) from the read;
        final String hg19Contig = hg19FirstOfPair.getContig();
        final int hg19Start = hg19FirstOfPair.getStart();
        final int hg19End = hg19FirstOfPair.getEnd();

        final GATKRead hg38FirstOfPair = input2ReadPair.getFirstOfPair();
        // There has to be a way to get a simple interval of the insert (or read) from the read;
        final String hg38Contig = hg38FirstOfPair.getContig();
        final int hg38Start = hg38FirstOfPair.getStart();
        final int hg38End = hg38FirstOfPair.getEnd();

        StringJoiner sj = new StringJoiner("\t");
        // Consider moving the Integer.toString() to where those variables are obtained
        sj.add(hg19Contig).add(Integer.toString(hg19Start)).add(Integer.toString(hg19End))
                .add(hg38Contig).add(Integer.toString(hg38Start)).add(Integer.toString(hg38End));
        tableWriter.println(sj);
    }

    @Override
    public void writeSummary(File outputTable, SAMFileHeader header) {
        try (PrintWriter pw = new PrintWriter(outputTable)){
            final String sample = header.getReadGroups().get(0).getSample();
            pw.println("#sample=" + sample);
            pw.println("onlyInHg19," + onlyInHg19);
            pw.println("onlyInHg38," + onlyInHg38);
        } catch (IOException e){
            throw new UserException("Could not write to " + outputTable, e);
        }

    }
}
