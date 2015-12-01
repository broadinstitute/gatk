package org.broadinstitute.hellbender.tools;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.Defaults;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.rnaseq.OverhangFixingManager;
import org.broadinstitute.hellbender.transformers.NDNCigarReadTransformer;
import org.broadinstitute.hellbender.transformers.ReadTransformer;
import org.broadinstitute.hellbender.utils.GenomeLocParser;
import org.broadinstitute.hellbender.utils.clipping.ReadClipper;
import org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.hellbender.utils.iterators.SAMRecordToReadIterator;
import org.broadinstitute.hellbender.utils.read.CigarUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.read.SAMFileGATKReadWriter;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.stream.StreamSupport;

import static org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions.*;

/**
 *
 * Splits reads that contain Ns in their cigar string (e.g. spanning splicing events).
 *
 * Identifies all N cigar elements and creates k+1 new reads (where k is the number of N cigar elements).
 * The first read includes the bases that are to the left of the first N element, while the part of the read that is to the right of the N
 * (including the Ns) is hard clipped and so on for the rest of the new reads.
 */
@CommandLineProgramProperties(
        summary = "Splits reads that contain Ns in their cigar string (e.g. spanning splicing events).",
        oneLineSummary = "Split Reads with N in Cigar",
        programGroup = ReadProgramGroup.class
)
public final class SplitNCigarReads extends CommandLineProgram {

    @Argument(fullName = INPUT_LONG_NAME, shortName= INPUT_SHORT_NAME, doc="The SAM/BAM file to read from.")
    public File INPUT;

    @Argument(fullName = OUTPUT_LONG_NAME, shortName = OUTPUT_SHORT_NAME, doc="Write output to this BAM filename instead of STDOUT")
    protected File OUTPUT;

    @Argument(fullName = StandardArgumentDefinitions.REFERENCE_LONG_NAME, shortName = StandardArgumentDefinitions.REFERENCE_SHORT_NAME, doc = "Reference sequence file.",
            common = true, optional = true)
    public File REFERENCE_SEQUENCE = Defaults.REFERENCE_FASTA;

    /**
     * This flag tells GATK to refactor cigar string with NDN elements to one element. It intended primarily for use in
     * a RNAseq pipeline since the problem might come up when using RNAseq aligner such as Tophat2 with provided transcriptoms.
     * You should only use this if you know that your reads have that problem.
     */
    @Argument(fullName = "refactor_NDN_cigar_string", shortName = "fixNDN", doc = "refactor cigar string with NDN elements to one element", optional = true)
    public boolean REFACTOR_NDN_CIGAR_READS = false;

    /**
     * For expert users only!  To minimize memory consumption you can lower this number, but then the tool may skip
     * overhang fixing in regions with too much coverage.  Just make sure to give Java enough memory!  4Gb should be
     * enough with the default value.
     */
    @Argument(fullName="maxReadsInMemory", shortName="maxInMemory", doc="max reads allowed to be kept in memory at a time by the BAM writer", optional=true)
    protected int MAX_RECORDS_IN_MEMORY = 150000;

    /**
     * If there are more than this many mismatches within the overhang regions, the whole overhang will get hard-clipped out.
     * It is still possible in some cases that the overhang could get clipped if the number of mismatches do not exceed this
     * value, e.g. if most of the overhang mismatches.
     */
    @Argument(fullName="maxMismatchesInOverhang", shortName="maxMismatches", doc="max number of mismatches allowed in the overhang", optional=true)
    protected int MAX_MISMATCHES_IN_OVERHANG = 1;

    /**
     * If there are more than this many bases in the overhang, we won't try to hard-clip them out
     */
    @Argument(fullName="maxBasesInOverhang", shortName="maxOverhang", doc="max number of bases allowed in the overhang", optional=true)
    protected int MAX_BASES_TO_CLIP = 40;

    @Argument(fullName="doNotFixOverhangs", shortName="doNotFixOverhangs", doc="do not have the walker hard-clip overhanging sections of the reads", optional=true)
    protected boolean doNotFixOverhangs = false;

    /**
     * This stores all of the already-split reads and manages any processing (e.g. clipping overhangs) that happens to them.
     * It will emit reads to the underlying writer as needed so we don't need to worry about any of that in this class.
     */
    protected OverhangFixingManager overhangManager;

    @Override
    protected Object doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        final SamReader in = SamReaderFactory.makeDefault().open(INPUT);
        final Iterable<GATKRead> readIter = new SAMRecordToReadIterator(in.iterator());
        final SAMFileGATKReadWriter outputWriter = initialize(in);

        final ReadTransformer rnaReadTransform = REFACTOR_NDN_CIGAR_READS ? new NDNCigarReadTransformer() : ReadTransformer.identity();

        StreamSupport.stream(readIter.spliterator(), false)
                .map(rnaReadTransform)
                .forEach(read -> splitNCigarRead(read, overhangManager));
        overhangManager.close();
        CloserUtil.close(in);
        CloserUtil.close(outputWriter);
        return null;
    }

    private SAMFileGATKReadWriter initialize(final SamReader in) {
        final SAMFileHeader outputHeader = ReadUtils.cloneSAMFileHeader(in.getFileHeader());
        final SAMFileGATKReadWriter outputWriter = new SAMFileGATKReadWriter(ReadUtils.createCommonSAMWriter(OUTPUT, REFERENCE_SEQUENCE, outputHeader, true, false, false));

        try {
            final IndexedFastaSequenceFile referenceReader = new CachingIndexedFastaSequenceFile(REFERENCE_SEQUENCE);
            GenomeLocParser genomeLocParser= new GenomeLocParser(referenceReader.getSequenceDictionary());
            overhangManager = new OverhangFixingManager(outputHeader, outputWriter, genomeLocParser, referenceReader, MAX_RECORDS_IN_MEMORY, MAX_MISMATCHES_IN_OVERHANG, MAX_BASES_TO_CLIP, doNotFixOverhangs);
            return outputWriter;
        } catch (FileNotFoundException ex) {
            throw new UserException.CouldNotReadInputFile(REFERENCE_SEQUENCE, ex);
        }
    }

    /**
     * Goes through the cigar string of the read and create new reads for each consecutive non-N elements (while hard clipping the rest of the read).
     * For example: for a read with cigar '1H2M2D1M2N1M2I1N1M2S' 3 new reads will be created with cigar strings: 1H2M2D1M, 1M2I and 1M2S
     *
     * @param read     the read to split
     */
    public static GATKRead splitNCigarRead(final GATKRead read, OverhangFixingManager manager) {
        final int numCigarElements = read.getCigar().numCigarElements();

        int firstCigarIndex = 0;
        for ( int i = 0; i < numCigarElements; i++ ) {
            final CigarElement cigarElement = read.getCigar().getCigarElement(i);
            if (cigarElement.getOperator() == CigarOperator.N) {
                manager.addRead(splitReadBasedOnCigar(read, firstCigarIndex, i, manager));
                firstCigarIndex = i+1;
            }
        }

        // if there are no N's in the read
        if (firstCigarIndex == 0) {
            manager.addRead(read);
        }
        //add the last section of the read: from the last N to the the end of the read
        // (it will be done for all the usual cigar string that does not end with N)
        else if (firstCigarIndex < numCigarElements) {
            manager.addRead(splitReadBasedOnCigar(read, firstCigarIndex, numCigarElements, null));
        }
        return read;
    }

    /**
     * Pull out an individual split position for a read
     *
     * @param read               the read being split
     * @param cigarStartIndex    the index of the first cigar element to keep
     * @param cigarEndIndex      the index of the last cigar element to keep
     * @param forSplitPositions  the manager for keeping track of split positions; can be null
     * @return a non-null read representing the section of the original read being split out
     */
    private static GATKRead splitReadBasedOnCigar(final GATKRead read, final int cigarStartIndex, final int cigarEndIndex, final OverhangFixingManager forSplitPositions) {
        int cigarFirstIndex = cigarStartIndex;
        int cigarSecondIndex = cigarEndIndex;

        //in case a section of the read ends or starts with D (for example the first section in 1M1D1N1M is 1M1D), we should trim this cigar element
        // it can be 'if', but it was kept as 'while' to make sure the code can work with Cigar strings that were not "cleaned"
        while(read.getCigar().getCigarElement(cigarFirstIndex).getOperator().equals(CigarOperator.D)) {
            cigarFirstIndex++;
        }
        while(read.getCigar().getCigarElement(cigarSecondIndex-1).getOperator().equals(CigarOperator.D)) {
            cigarSecondIndex--;
        }
        if(cigarFirstIndex > cigarSecondIndex) {
            throw new IllegalArgumentException("Cannot split this read (might be an empty section between Ns, for example 1N1D1N): " + read.getCigar().toString());
        }

        // we keep only the section of the read that is aligned to the reference between startRefIndex and stopRefIndex (inclusive).
        // the other sections of the read are clipped:
        final int startRefIndex = read.getUnclippedStart() + CigarUtils.countRefBasesBasedOnCigar(read, 0, cigarFirstIndex); //goes through the prefix of the cigar (up to cigarStartIndex) and move the reference index.
        final int stopRefIndex = startRefIndex + CigarUtils.countRefBasesBasedOnCigar(read,cigarFirstIndex,cigarSecondIndex)-1; //goes through a consecutive non-N section of the cigar (up to cigarEndIndex) and move the reference index.

        if ( forSplitPositions != null ) {
            final String contig = read.getContig();
            final int splitStart = startRefIndex + CigarUtils.countRefBasesBasedOnCigar(read,cigarFirstIndex,cigarEndIndex);  //we use cigarEndIndex instead of cigarSecondIndex so we won't take into account the D's at the end.
            final int splitEnd = splitStart + read.getCigar().getCigarElement(cigarEndIndex).getLength() - 1;
            forSplitPositions.addSplicePosition(contig, splitStart, splitEnd);
        }

        return ReadClipper.hardClipToRegionIncludingClippedBases(read, startRefIndex, stopRefIndex);
    }

}
