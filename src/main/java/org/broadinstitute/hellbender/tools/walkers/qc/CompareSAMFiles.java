package org.broadinstitute.hellbender.tools.walkers.qc;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.util.PeekableIterator;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.ShortVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadQueryNameComparator;

import java.io.File;
import java.net.UnknownServiceException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * This should eventually be turned into a biread walker, a subclass of WalkerBase.
 */

@CommandLineProgramProperties(
        summary = "Given a query name sorted bam, call it read1, and another query sorted bam, call it read2. " +
                "Traverse read1 with the information about matching read2 reads",
        oneLineSummary = "s",
        programGroup = ShortVariantDiscoveryProgramGroup.class // Sato: change
)
public class CompareSAMFiles extends GATKTool {

    @Argument(fullName = "read2", doc = "read2")
    public GATKPath read2;

    PeekableIterator<GATKRead> read1Iterator;
    PeekableIterator<GATKRead> read2Iterator;

    GATKRead currentRead1;
    GATKRead currentRead2;
    ReadQueryNameComparator queryNameComparator;

    long genomeDupTranscrNot = 0;
    long bothDup = 0;
    long genomeNotTranscrDup = 0;
    long neitherDup;
    long genomeDupNotFoundInTranscr = 0;
    long notFoundInTranscr = 0;
    long notFoundInGenome = 0;
    long numReadPairsGenome = 0;
    long numReadPairsTranscr = 0;

    ReadsDataSource reads2;
    public void onTraversalStart(){
        queryNameComparator = new ReadQueryNameComparator();
        read1Iterator = new PeekableIterator<>(directlyAccessEngineReadsDataSource().iterator());
        SAMFileHeader.SortOrder so1 = directlyAccessEngineReadsDataSource().getHeader().getSortOrder();

        reads2 = new ReadsPathDataSource(read2.toPath());
        SAMFileHeader.SortOrder so2 = reads2.getHeader().getSortOrder();
        read2Iterator = new PeekableIterator<>(reads2.iterator());
    }

    boolean read1Behind;

    private class ReadPair {
        private GATKRead firstOfPair;
        private GATKRead secondOfPair;
        private List<GATKRead> secondaryAlignments = new ArrayList<>(10);
        private List<GATKRead> supplementaryAlignments = new ArrayList<>(10); // Finally understand the difference
        private ReadPair(){}

        private void add(GATKRead read){
            if (read.isFirstOfPair()){
                this.firstOfPair = read;
            } else if (read.isSecondOfPair()){
                this.secondOfPair = read;
            } else if (read.isSupplementaryAlignment()){
                this.secondaryAlignments.add(read);
            } else {
                int d = 3;
                throw new UserException("Unknown read type");
            }
        }

        public boolean isDuplicateMarked(){
            // Doing some investigation
            if (firstOfPair.isDuplicate()){
                // Make sure the rest is duplicate-marked
                if (!secondOfPair.isDuplicate() || secondaryAlignments.stream().anyMatch(r -> ! r.isDuplicate())){
                    throw new UserException("First of pair a duplicate but the rest is not" + secondOfPair.getName());
                }
            } else {
                // Make sure the rest is not duplicate-marked
                if (secondOfPair.isDuplicate() || secondaryAlignments.stream().anyMatch(r -> r.isDuplicate())){
                    throw new UserException("First of pair a not duplicate but the rest is " + secondOfPair.getName());
                }
            }
            return firstOfPair.isDuplicate();
        }
    }

    @Override
    public void traverse() {
        currentRead1 = read1Iterator.next();
        currentRead2 = read2Iterator.next(); // Check for hasNext()

        while (read1Iterator.hasNext()){
            final ReadPair input1ReadPair = new ReadPair();
            final ReadPair input2ReadPair = new ReadPair();
            final GATKRead firstInput1Read = currentRead1;
            final GATKRead firstInput2Read = currentRead2;

            int diff = queryNameComparator.compareReadNames(currentRead1, currentRead2);
            if (diff == 0){
                // fast-forward and gather all the reads that share the same query name
                // i.e. the mate and secondary and supplementary alignments
                input1ReadPair.add(firstInput1Read);
                input2ReadPair.add(firstInput2Read);

                // collectQueryNameGroup() advances the iterator and updates the ReadPair object by side effect
                // maybe avoid that, but OK for now.
                collectQueryNameGroup(read1Iterator, input1ReadPair, firstInput1Read);
                collectQueryNameGroup(read2Iterator, input2ReadPair, firstInput2Read);

                // Process the info. Abstract this section and put the ribosome part in.
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

                currentRead1 = read1Iterator.next();
                currentRead2 = read2Iterator.next();
            } else if (diff < 0){
                // Read1 is behind. This queryname group is not in input2.
                int d = 3;
                input1ReadPair.add(firstInput1Read);
                collectQueryNameGroup(read1Iterator, input1ReadPair, firstInput1Read);
                // Extract a method, like "handleWhenDiffLessThan0"
                if (input1ReadPair.isDuplicateMarked()){
                    genomeDupNotFoundInTranscr += 1;
                } else {
                    notFoundInTranscr += 1;
                }
                currentRead1 = read1Iterator.next();
            } else {
                int d = 3;
                // read2 not found in read1---what to do?
                currentRead2 = read2Iterator.next();

            }
        }
    }

    @Override
    public Object onTraversalSuccess(){
        File out = new File("/Users/tsato/workspace/gatk/tmp/duplicates.txt");

        System.out.println("genomeDupTranscrNot," + genomeDupTranscrNot);
        System.out.println("bothDup," + bothDup);
        System.out.println("genomeNotTranscrDup," + genomeNotTranscrDup);
        System.out.println("neitherDup," + neitherDup);
        System.out.println("genomeDupNotFoundInTranscr," + genomeDupNotFoundInTranscr);
        System.out.println("notFoundInTranscr," + notFoundInTranscr);
        System.out.println("notFoundInGenome," + notFoundInGenome);
        System.out.println("numReadPairsGenome," + numReadPairsGenome);
        System.out.println("numReadPairsTranscr," + numReadPairsTranscr);

        return "SUCCESS";
    }

    // Assumes that templateRead has already been added to the ReadPair object
    private void collectQueryNameGroup(PeekableIterator<GATKRead> iterator, ReadPair readPair, GATKRead templateRead) {
        GATKRead currentRead = templateRead;
        GATKRead nextRead = iterator.peek();
        while (iterator.hasNext() && templateRead.getName().equals(nextRead.getName())){ // Names match
            currentRead = iterator.next();
            readPair.add(currentRead);
            nextRead = iterator.peek();
        }
    }
}
