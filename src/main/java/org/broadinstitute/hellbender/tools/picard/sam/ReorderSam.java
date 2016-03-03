package org.broadinstitute.hellbender.tools.picard.sam;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.PicardCommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.io.File;
import java.util.HashMap;
import java.util.Map;

/**
 * Reorders a SAM/BAM input file according to the order of contigs in a second reference sequence
 *
 * @author mdepristo
 */
@CommandLineProgramProperties(
        summary = "Not to be confused with SortSam which sorts a SAM/BAM file with a valid sequence dictionary, " +
                "ReorderSam reorders reads in a SAM/BAM file to match the contig ordering in a provided reference file, " +
                "as determined by exact name matching of contigs.  Reads mapped to contigs absent in the new " +
                "reference are dropped. Runs substantially faster if the input is an indexed BAM file.",
        oneLineSummary = "Reorders reads in a SAM/BAM file to match ordering in reference",
        programGroup = ReadProgramGroup.class
)
public final class ReorderSam extends PicardCommandLineProgram {

    @Argument(fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME, doc = "Input file SAM.BAM to extract reads from.")
    public File INPUT;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "Output file SAM/BAM to write extracted reads to.")
    public File OUTPUT;

    @Argument(shortName = "S", doc = "If true, then allows only a partial overlap of the BAM contigs with the new reference " +
            "sequence contigs.  By default, this tool requires a corresponding contig in the new " +
            "reference for each read contig")
    public boolean ALLOW_INCOMPLETE_DICT_CONCORDANCE = false;

    @Argument(shortName = "U", doc = "If true, then permits mapping from a read contig to a new reference contig with the " +
            "same name but a different length.  Highly dangerous, only use if you know what you " +
            "are doing.")
    public boolean ALLOW_CONTIG_LENGTH_DISCORDANCE = false;

    @Override
    protected Object doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsReadable(REFERENCE_SEQUENCE);
        IOUtil.assertFileIsWritable(OUTPUT);

        final SamReader in = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(INPUT);

        ReferenceSequenceFile reference = ReferenceSequenceFileFactory.getReferenceSequenceFile(REFERENCE_SEQUENCE);
        SAMSequenceDictionary refDict = reference.getSequenceDictionary();

        if (refDict == null) {
            CloserUtil.close(in);
            throw new UserException("No reference sequence dictionary found. Aborting. " +
                    "You can create a sequence dictionary for the reference fasta using CreateSequenceDictionary.jar.");
        }

        printDictionary("SAM/BAM file", in.getFileHeader().getSequenceDictionary());
        printDictionary("Reference", refDict);
        Map<Integer, Integer> newOrder = buildSequenceDictionaryMap(refDict, in.getFileHeader().getSequenceDictionary());

        // has to be after we create the newOrder
        SAMFileHeader outHeader = ReadUtils.cloneSAMFileHeader(in.getFileHeader());
        outHeader.setSequenceDictionary(refDict);

        logger.info("Writing reads...");
        if (in.hasIndex()) {
            try (final SAMFileWriter out = createSAMWriter(OUTPUT, REFERENCE_SEQUENCE, outHeader, true)) {

                // write the reads in contig order
                for (final SAMSequenceRecord contig : refDict.getSequences()) {
                    final SAMRecordIterator it = in.query(contig.getSequenceName(), 0, 0, false);
                    writeReads(out, it, newOrder, contig.getSequenceName());
                }
                // don't forget the unmapped reads
                writeReads(out, in.queryUnmapped(), newOrder, "unmapped");
            }
        } else {
            try (final SAMFileWriter out = createSAMWriter(OUTPUT, REFERENCE_SEQUENCE, outHeader, false)) {
                writeReads(out, in.iterator(), newOrder, "All reads");
            }
        }

        // cleanup
        CloserUtil.close(in);
        return null;
    }

    /**
     * Low-level helper function that returns the new reference index for oldIndex according to the
     * ordering map newOrder.  Read is provided in case an error occurs, so that an informative message
     * can be made.
     */
    private int newOrderIndex(SAMRecord read, int oldIndex, Map<Integer, Integer> newOrder) {
        if (oldIndex == -1)
            return -1; // unmapped read
        else {
            final Integer n = newOrder.get(oldIndex);
            if (n == null) throw new UserException("BUG: no mapping found for read " + read.getSAMString());
            else return n;
        }
    }

    /**
     * Helper function that writes reads from iterator it into writer out, updating each SAMRecord along the way
     * according to the newOrder mapping from dictionary index -> index.  Name is used for printing only.
     */
    private void writeReads(final SAMFileWriter out,
                            final SAMRecordIterator it,
                            final Map<Integer, Integer> newOrder,
                            final String name) {
        long counter = 0;
        logger.info("  Processing " + name);

        while (it.hasNext()) {
            counter++;
            final SAMRecord read = it.next();
            final int oldRefIndex = read.getReferenceIndex();
            final int oldMateIndex = read.getMateReferenceIndex();
            final int newRefIndex = newOrderIndex(read, oldRefIndex, newOrder);

            read.setHeaderStrict(out.getFileHeader());
            read.setReferenceIndex(newRefIndex);

            final int newMateIndex = newOrderIndex(read, oldMateIndex, newOrder);
            if (oldMateIndex != -1 && newMateIndex == -1) { // becoming unmapped
                read.setMateAlignmentStart(0);
                read.setMateUnmappedFlag(true);
                read.setAttribute(SAMTag.MC.name(), null);      // Set the Mate Cigar String to null
            }
            read.setMateReferenceIndex(newMateIndex);

            out.addAlignment(read);
        }

        it.close();
        logger.info("Wrote " + counter + " reads");
    }

    /**
     * Constructs a mapping from read sequence records index -> new sequence dictionary index for use in
     * reordering the reference index and mate reference index in each read.  -1 means unmapped.
     */
    private Map<Integer, Integer> buildSequenceDictionaryMap(final SAMSequenceDictionary refDict,
                                                             final SAMSequenceDictionary readsDict) {
        Map<Integer, Integer> newOrder = new HashMap<>();

        logger.info("Reordering SAM/BAM file:");
        for (final SAMSequenceRecord refRec : refDict.getSequences()) {
            final SAMSequenceRecord readsRec = readsDict.getSequence(refRec.getSequenceName());

            if (readsRec != null) {
                if (refRec.getSequenceLength() != readsRec.getSequenceLength()) {
                    String msg = String.format("Discordant contig lengths: read %s LN=%d, ref %s LN=%d",
                            readsRec.getSequenceName(), readsRec.getSequenceLength(),
                            refRec.getSequenceName(), refRec.getSequenceLength());
                    if (ALLOW_CONTIG_LENGTH_DISCORDANCE) {
                        logger.warn(msg);
                    } else {
                        throw new UserException(msg);
                    }
                }

                logger.info(String.format("  Reordering read contig %s [index=%d] to => ref contig %s [index=%d]%n",
                        readsRec.getSequenceName(), readsRec.getSequenceIndex(),
                        refRec.getSequenceName(), refRec.getSequenceIndex()));
                newOrder.put(readsRec.getSequenceIndex(), refRec.getSequenceIndex());
            }
        }

        for (SAMSequenceRecord readsRec : readsDict.getSequences()) {
            if (!newOrder.containsKey(readsRec.getSequenceIndex())) {
                if (ALLOW_INCOMPLETE_DICT_CONCORDANCE)
                    newOrder.put(readsRec.getSequenceIndex(), -1);
                else
                    throw new UserException("New reference sequence does not contain a matching contig for " + readsRec.getSequenceName());
            }
        }

        return newOrder;
    }

    /**
     * Helper function to print out a sequence dictionary
     */
    private void printDictionary(String name, SAMSequenceDictionary dict) {
        logger.info(name);
        for (final SAMSequenceRecord contig : dict.getSequences()) {
            logger.info("  SN=%s LN=%d%n", contig.getSequenceName(), contig.getSequenceLength());
        }
    }
}
