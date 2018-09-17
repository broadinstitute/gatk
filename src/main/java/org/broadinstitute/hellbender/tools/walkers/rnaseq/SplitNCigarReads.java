package org.broadinstitute.hellbender.tools.walkers.rnaseq;

import htsjdk.samtools.*;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.TwoPassReadWalker;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.transformers.MappingQualityReadTransformer;
import org.broadinstitute.hellbender.transformers.NDNCigarReadTransformer;
import org.broadinstitute.hellbender.transformers.ReadTransformer;
import org.broadinstitute.hellbender.utils.GenomeLocParser;
import org.broadinstitute.hellbender.utils.SATagBuilder;
import org.broadinstitute.hellbender.utils.clipping.ReadClipper;
import org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.hellbender.utils.read.*;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import static org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions.OUTPUT_LONG_NAME;
import static org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions.OUTPUT_SHORT_NAME;

/**
 *
 * Splits reads that contain Ns in their cigar string (e.g. spanning splicing events in RNAseq data).
 *
 * Identifies all N cigar elements and creates k+1 new reads (where k is the number of N cigar elements).
 * The first read includes the bases that are to the left of the first N element, while the part of the read that is to the right of the N
 * (including the Ns) is hard clipped and so on for the rest of the new reads. Used for post-processing RNA reads aligned against the full reference.
 *
 * <h3>Input</h3>
 *  <p>
 *	    BAM file
 *  </p>
 *
 *
 * <h3>Output</h3>
 *  <p>
 *      BAM file with reads split at N CIGAR elements and CIGAR strings updated.
 *  </p>
 *
 * <h3>Usage example</h3>
 *  <pre>
 *    gatk SplitNCigarReads \
 *      -R Homo_sapiens_assembly38.fasta \
 *      -I input.bam \
 *      -O output.bam
 *  </pre>
 */
@DocumentedFeature
@CommandLineProgramProperties(
        summary = "Splits reads that contain Ns in their cigar string (e.g. spanning splicing events).",
        oneLineSummary = "Split Reads with N in Cigar",
        programGroup = ReadDataManipulationProgramGroup.class
)
public final class SplitNCigarReads extends TwoPassReadWalker {

    // A list of tags that break upon splitting on N. These will be removed from reads in the output.
    // NOTE: for future developers who want to use these tags. For each tag you remove from this list corresponding
    // support for properly recalculating the tags must be added to repairSupplementaryTags()
    static final String[] TAGS_TO_REMOVE = {"NM","MD","NH"};
    static final String MATE_CIGAR_TAG = "MC";

    @Argument(fullName = OUTPUT_LONG_NAME, shortName = OUTPUT_SHORT_NAME, doc="Write output to this BAM filename")
    String OUTPUT;

    /**
     * This flag tells GATK to refactor cigar string with NDN elements to one element. It intended primarily for use in
     * a RNAseq pipeline since the problem might come up when using RNAseq aligner such as Tophat2 with provided transcriptomes.
     * You should only use this if you know that your reads have that problem.
     */
    @Argument(fullName = "refactor-cigar-string", shortName = "fixNDN", doc = "refactor cigar string with NDN elements to one element", optional = true)
    boolean REFACTOR_NDN_CIGAR_READS = false;

    /**
     * This flag turns off the mapping quality 255 -> 60 read transformer. The transformer is on by default to ensure that
     * uniquely mapping reads assigned STAR's default 255 MQ aren't filtered out by HaplotypeCaller.
     */
    @Argument(fullName = "skip-mapping-quality-transform", shortName = "skip-mq-transform", doc = "skip the 255 -> 60 MQ read transform", optional = true)
    boolean SKIP_MQ_TRANSFORM = false;

    /**
     * For expert users only!  To minimize memory consumption you can lower this number, but then the tool may skip
     * overhang fixing in regions with too much coverage.  Just make sure to give Java enough memory!  4Gb should be
     * enough with the default value.
     */
    @Advanced
    @Argument(fullName="max-reads-in-memory", doc="max reads allowed to be kept in memory at a time by the BAM writer", optional=true)
    int MAX_RECORDS_IN_MEMORY = 150000;

    /**
     * If there are more than this many mismatches within the overhang regions, the whole overhang will get hard-clipped out.
     * It is still possible in some cases that the overhang could get clipped if the number of mismatches do not exceed this
     * value, e.g. if most of the overhang mismatches.
     */
    @Argument(fullName="max-mismatches-in-overhang", doc="max number of mismatches allowed in the overhang", optional=true)
    int MAX_MISMATCHES_IN_OVERHANG = 1;

    /**
     * If there are more than this many bases in the overhang, we won't try to hard-clip them out
     */
    @Argument(fullName="max-bases-in-overhang", doc="max number of bases allowed in the overhang", optional=true)
    int MAX_BASES_TO_CLIP = 40;

    @Argument(fullName="do-not-fix-overhangs", doc="do not have the walker soft-clip overhanging sections of the reads", optional=true)
    boolean doNotFixOverhangs = false;

    @Argument(fullName="process-secondary-alignments", doc="have the walker split secondary alignments (will still repair MC tag without it)", optional=true)
    boolean processSecondaryAlignments = false;

    @Override
    public boolean requiresReference() {
        return true;
    }

    private SAMFileGATKReadWriter outputWriter;
    private OverhangFixingManager overhangManager;
    private IndexedFastaSequenceFile referenceReader;
    SAMFileHeader header;

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        return Collections.singletonList(ReadFilterLibrary.ALLOW_ALL_READS);
    }

    @Override
    public ReadTransformer makePreReadFilterTransformer(){
        if (REFACTOR_NDN_CIGAR_READS) {
            return new NDNCigarReadTransformer();
        }
        else {
            return ReadTransformer.identity();
        }
    }

    //FIXME: once the engine accepts read transformer arguments, remove these magic numbers?
    private static final int FROM_QUALITY = 255;
    private static final int TO_QUALITY = 60;


    @Override
    public ReadTransformer makePostReadFilterTransformer(){
        if (SKIP_MQ_TRANSFORM) {
            return ReadTransformer.identity();
        }
        else {
            return new MappingQualityReadTransformer(FROM_QUALITY, TO_QUALITY);
        }
    }

    @Override
    public void onTraversalStart() {
        header = getHeaderForSAMWriter();
        try {
            referenceReader = new CachingIndexedFastaSequenceFile(referenceArguments.getReferencePath());
            GenomeLocParser genomeLocParser = new GenomeLocParser(getBestAvailableSequenceDictionary());
            outputWriter = createSAMWriter(IOUtils.getPath(OUTPUT), false);
            overhangManager = new OverhangFixingManager(header, outputWriter, genomeLocParser, referenceReader, MAX_RECORDS_IN_MEMORY, MAX_MISMATCHES_IN_OVERHANG, MAX_BASES_TO_CLIP, doNotFixOverhangs, processSecondaryAlignments);

        } catch (FileNotFoundException ex) {
            throw new UserException.CouldNotReadInputFile(referenceArguments.getReferencePath(), ex);
        }
    }

    @Override
    protected void firstPassApply(GATKRead read, ReferenceContext bytes, FeatureContext featureContext) {
        splitNCigarRead(read,overhangManager, true, header, processSecondaryAlignments);
    }

    @Override
    protected void secondPassApply(GATKRead read, ReferenceContext bytes, FeatureContext featureContext) {
        splitNCigarRead(read,overhangManager, true, header, processSecondaryAlignments);
    }

    // Activates writing in the manager, which destinguishes each pass
    @Override
    protected void afterFirstPass() {
        overhangManager.activateWriting();
    }

    @Override
    public void closeTool() {
        if (overhangManager != null) { overhangManager.flush(); }
        if (outputWriter != null ) { outputWriter.close(); }
        try {if (referenceReader != null) { referenceReader.close(); } }
        catch (IOException ex) {
            throw new UserException.MissingReference("Could not find reference file");
        }
    }


    /**
     * Goes through the cigar string of the read and create new reads for each consecutive non-N elements (while soft clipping the rest of the read) that are supplemental to each other.
     * For example: for a read with cigar '1H2M2D1M2N1M2I1N1M2S' 3 new reads will be created with cigar strings: 1H2M2D1M6S, 1H3S1M2I2S and 1H6S1M2S
     * If the read has an MC tag it will be adjusted according to the clipping of that mate based on its cigar
     *
     * @param read     the read to split (can be null)
     * @param emitReads   a parameter used to mock behavior for repairing mate cigar string information
     */
    public static GATKRead splitNCigarRead(final GATKRead read, OverhangFixingManager manager, boolean emitReads, SAMFileHeader header, boolean secondaryAlignments) {
        final int numCigarElements = read.numCigarElements();
        List<GATKRead> splitReads = new ArrayList<>(2);

        // Run the tool on dummy mate read to determine what the mate cigar will be upon completion, if manager has a prediction then dont repair
        if (emitReads && read.hasAttribute(MATE_CIGAR_TAG)) {
            final GATKRead mateSplitting = splitNCigarRead(ArtificialReadUtils.createArtificialRead(header, TextCigarCodec.decode(read.getAttributeAsString(MATE_CIGAR_TAG))), manager, false, header, secondaryAlignments);
            read.setAttribute(MATE_CIGAR_TAG, mateSplitting.getCigar().toString());
        }
        manager.setPredictedMateInformation(read);

        // If it is a secondary alignment, repair its mate-information (assuming its mate was primary) and pass
        if ( !secondaryAlignments &&  read.isSecondaryAlignment()){
            manager.addReadGroup(Collections.singletonList(read));
            return read;
        }

        boolean sectionHasMatch = false;
        int firstCigarIndex = 0;
        for ( int i = 0; i < numCigarElements; i++ ) {
            final CigarElement cigarElement = read.getCigarElement(i);
            CigarOperator op = cigarElement.getOperator();

            // One of the "Real" operators
            if (op == CigarOperator.M || op == CigarOperator.EQ || op == CigarOperator.X ||
                    op == CigarOperator.I || op == CigarOperator.D) {
                sectionHasMatch = true;
            }

            if (op == CigarOperator.N) {
                if (sectionHasMatch) {
                    if (!emitReads) {
                        // not passing the manager ensures that no splices get added to the manager for fake reads
                        splitReads.add(splitReadBasedOnCigar(read, firstCigarIndex, i, null));
                    } else {
                        splitReads.add(splitReadBasedOnCigar(read, firstCigarIndex, i, manager));
                    }
                }
                firstCigarIndex = i+1;
                sectionHasMatch = false;
            }
        }

        // if there are no N's in the read
        if (splitReads.size() < 1) {
            if (emitReads) {
                manager.addReadGroup(Collections.singletonList(read));
            }
            return read;
        }
        //add the last section of the read: from the last N to the the end of the read
        // (it will be done for all the usual cigar string that does not end with N)
        else if ((firstCigarIndex < numCigarElements) && sectionHasMatch) {
            splitReads.add(splitReadBasedOnCigar(read, firstCigarIndex, numCigarElements, null));
        }

        if (emitReads) {
            manager.addReadGroup(splitReads);
            return read;

        // If emitReads is false then we want the mate of the read
        } else {
            return splitReads.get(0);
        }
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
        while(read.getCigarElement(cigarFirstIndex).getOperator().equals(CigarOperator.D)) {
            cigarFirstIndex++;
        }
        while(read.getCigarElement(cigarSecondIndex-1).getOperator().equals(CigarOperator.D)) {
            cigarSecondIndex--;
        }
        if(cigarFirstIndex > cigarSecondIndex) {
            throw new IllegalArgumentException("Cannot split this read (might be an empty section between Ns, for example 1N1D1N): " + read.getCigar().toString());
        }

        // we keep only the section of the read that is aligned to the reference between startRefIndex and stopRefIndex (inclusive).
        // the other sections of the read are clipped:
        final int startRefIndex = read.getUnclippedStart() + CigarUtils.countRefBasesBasedOnUnclippedAlignment(read, 0, cigarFirstIndex); //goes through the prefix of the cigar (up to cigarStartIndex) and move the reference index.
        final int stopRefIndex = startRefIndex + CigarUtils.countRefBasesBasedOnUnclippedAlignment(read,cigarFirstIndex,cigarSecondIndex)-1; //goes through a consecutive non-N section of the cigar (up to cigarEndIndex) and move the reference index.

        if ( forSplitPositions != null ) {
            final String contig = read.getContig();
            final int splitStart = startRefIndex + CigarUtils.countRefBasesBasedOnUnclippedAlignment(read,cigarFirstIndex,cigarEndIndex);  //we use cigarEndIndex instead of cigarSecondIndex so we won't take into account the D's at the end.
            final int splitEnd = splitStart + read.getCigarElement(cigarEndIndex).getLength() - 1;
            forSplitPositions.addSplicePosition(contig, splitStart, splitEnd);
        }

        return ReadClipper.softClipToRegionIncludingClippedBases(read,startRefIndex,stopRefIndex);
    }

    /**
     * A method that repairs the NH and NM tags for a group of reads
     *
     * @param readFamily         the a list of reads where the first item is the read to be marked as primary
     * @param header             the file header to associate with the given reads
     * @return a non-null read representing the section of the original read being split out
     */
    public static void repairSupplementaryTags(List<GATKRead> readFamily, SAMFileHeader header) {
        for (GATKRead read : readFamily) {
            for (String attribute : TAGS_TO_REMOVE) {
                read.clearAttribute(attribute);
            }
        }
        if (readFamily.size() > 1) {
            GATKRead primary = readFamily.remove(0);
            SATagBuilder.setReadsAsSupplemental(primary,readFamily);
        }
    }
}
