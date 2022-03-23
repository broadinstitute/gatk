package org.broadinstitute.hellbender.tools.walkers.qc;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.util.PeekableIterator;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.ShortVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadQueryNameComparator;
import org.broadinstitute.hellbender.utils.read.SAMFileGATKReadWriter;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;

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
    @Argument(fullName = "output", shortName = "O")
    public File outputTable = new File("duplicate_comaprison_counts.csv");

    @Argument(fullName = "read2", doc = "read2")
    public GATKPath read2;

    @Argument(fullName = "mode")
    public String mode = "RIBOSOME";

    @Argument(fullName = "out1")
    public File outSamFile1;

    @Argument(fullName = "out2")
    public File outSamFile2;

    PeekableIterator<GATKRead> read1Iterator;
    PeekableIterator<GATKRead> read2Iterator;

    GATKRead currentRead1;
    GATKRead currentRead2;
    ReadQueryNameComparator queryNameComparator;
    QuerynameSetComparison queryNameSetComparison;

    ReadsDataSource reads2;

    SAMFileGATKReadWriter writer1;
    SAMFileGATKReadWriter writer2;

    public void onTraversalStart(){
        queryNameComparator = new ReadQueryNameComparator();
        read1Iterator = new PeekableIterator<>(directlyAccessEngineReadsDataSource().iterator());
        SAMFileHeader.SortOrder so1 = directlyAccessEngineReadsDataSource().getHeader().getSortOrder();

        reads2 = new ReadsPathDataSource(read2.toPath());
        SAMFileHeader.SortOrder so2 = reads2.getHeader().getSortOrder();
        read2Iterator = new PeekableIterator<>(reads2.iterator());

        // Clean this Unhandled Exception mess somehow
        try {
            queryNameSetComparison = getQuerynameComparison(mode);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        if (outSamFile1 != null){
            Utils.nonNull(outSamFile2, "If out1 is specified, out2 must both be non-null");
            writer1 = createSAMWriter(new GATKPath(outSamFile1.getAbsolutePath()), false);
            writer2 = createSAMWriter(new GATKPath(outSamFile2.getAbsolutePath()), false);
        }

    }

    private QuerynameSetComparison getQuerynameComparison(String mode) throws FileNotFoundException {
        if (mode.equals("RIBOSOME")){
            return new CompareRibosomalReads(new File("ribosome.tsv"));
        } else if (mode.equals("DUPLICATE")){
            return new CompareDuplicateMarking();
        } else if (mode.equals("DIFF")) {
            Utils.nonNull(outSamFile1, "out1 must be specified under the DIFF mode");
            Utils.nonNull(outSamFile2, "out2 must be specified under the DIFF mode");
            return new SamDiff(writer1, writer2);
        } else {
            throw new UserException("Unrecognizable mode " + mode);
        }
    }

    boolean read1Behind;

    @Override
    public void traverse() {
        currentRead1 = read1Iterator.next();
        currentRead2 = read2Iterator.next(); // Check for hasNext()

        while (read1Iterator.hasNext() && read2Iterator.hasNext()){
            final ReadPair input1ReadPair = new ReadPair();
            final ReadPair input2ReadPair = new ReadPair();

            // First in the sense that these are the first read encountered in a querynameGroup
            final GATKRead firstInput1Read = currentRead1;
            final GATKRead firstInput2Read = currentRead2;

            if (currentRead1 == null || currentRead2 == null){
                logger.info("currentRead1 is null");
                logger.info("read1Iterator has next = " + read1Iterator.hasNext());
                logger.info("read2Iterator has next = " + read2Iterator.hasNext());
            }

            if (currentRead1 == null || currentRead2 == null){
                logger.info("currentRead2 is null");
                logger.info("read1Iterator has next = " + read1Iterator.hasNext());
                logger.info("read2Iterator has next = " + read2Iterator.hasNext());
            }

            int diff = queryNameComparator.compareReadNames(currentRead1, currentRead2);

            if (diff == 0){
                // The query names match:
                // Fast-forward and gather all the reads that share the same query name (queryname set)
                // i.e. the mate and secondary and supplementary alignments
                input1ReadPair.add(firstInput1Read);
                input2ReadPair.add(firstInput2Read);

                // collectQueryNameGroup() advances the iterator and updates the ReadPair object by side effect
                // maybe avoid that, but OK for now.
                collectQueryNameGroup(read1Iterator, input1ReadPair, firstInput1Read);
                collectQueryNameGroup(read2Iterator, input2ReadPair, firstInput2Read);

                // Process the info. Abstract this section and put the ribosome part in.
                queryNameSetComparison.processMatchingQuerynameSets(input1ReadPair, input2ReadPair);

                if (read1Iterator.hasNext()){
                    currentRead1 = read1Iterator.next();
                }

                if (read2Iterator.hasNext()){
                    currentRead2 = read2Iterator.next();
                }
            } else if (diff < 0){
                // Read1 is behind. This queryname group is not in input2.
                input1ReadPair.add(firstInput1Read);
                collectQueryNameGroup(read1Iterator, input1ReadPair, firstInput1Read);
                // Extract a method, like "handleWhenDiffLessThan0"
                queryNameSetComparison.processInput1(input1ReadPair);

                if (read1Iterator.hasNext()){
                    currentRead1 = read1Iterator.next();
                }
            } else {
                logger.info("Found transcriptome reads not in the genome");
                logger.info(currentRead2.getName());

                input2ReadPair.add(firstInput2Read);
                collectQueryNameGroup(read2Iterator, input2ReadPair, firstInput2Read);
                // Ditto above. Extract this method.
                queryNameSetComparison.processInput2(input2ReadPair);

                if (read2Iterator.hasNext()){
                    currentRead2 = read2Iterator.next();
                }
            }
        }


    }

    @Override
    public Object onTraversalSuccess(){
        if (read1Iterator.hasNext()){
            logger.info("read1Iterator has next");
            while (read1Iterator.hasNext()){
                final ReadPair input1ReadPair = new ReadPair();
                final GATKRead firstInput1Read = currentRead1;
                input1ReadPair.add(firstInput1Read);
                collectQueryNameGroup(read1Iterator, input1ReadPair, firstInput1Read);
                queryNameSetComparison.processInput1(input1ReadPair);
                currentRead1 = read1Iterator.next();
            }
        }

        if (read2Iterator.hasNext()){
            logger.info("read2Iterator has next");
            while (read2Iterator.hasNext()){ // CODE DUPLICATION
                final ReadPair input2ReadPair = new ReadPair();
                final GATKRead firstInput2Read = currentRead2;
                input2ReadPair.add(firstInput2Read);
                collectQueryNameGroup(read2Iterator, input2ReadPair, firstInput2Read);
                queryNameSetComparison.processInput2(input2ReadPair);
                currentRead2 = read2Iterator.next();
            }
        }

        queryNameSetComparison.writeSummary(outputTable, getHeaderForReads());

        if (writer1 != null) {
            writer1.close();
        } if (writer2 != null){
            writer2.close();
        }

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
