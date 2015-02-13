/*
* Copyright (c) 2012 The Broad Institute
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.hellbender.tools;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.StringUtil;
import org.apache.commons.lang3.tuple.MutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.ClippingOp;
import org.broadinstitute.hellbender.utils.ClippingRepresentation;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.clipping.ReadClipper;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Read clipping based on quality, position or sequence matching
 *
 * <p>This tool provides simple, powerful read clipping capabilities that allow you to remove low quality strings of bases, sections of reads, and reads containing user-provided sequences.</p> 
 *
 * <p>There are three arguments for clipping (quality, position and sequence), which can be used alone or in combination. In addition, you can also specify a clipping representation, which determines exactly how ClipReads applies clips to the reads (soft clips, writing Q0 base quality scores, etc.). Please note that you MUST specify at least one of the three clipping arguments, and specifying a clipping representation is not sufficient. If you do not specify a clipping argument, the program will run but it will not do anything to your reads.</p>
 *
 * <dl>
 *     <dt>Quality score based clipping</dt>
 *     <dd>
 *         Clip bases from the read in clipper from
 *         <pre>argmax_x{ \sum{i = x + 1}^l (qTrimmingThreshold - qual)</pre>
 *         to the end of the read.  This is copied from BWA.
 *
 *         Walk through the read from the end (in machine cycle order) to the beginning, calculating the
 *         running sum of qTrimmingThreshold - qual.  While we do this, we track the maximum value of this
 *         sum where the delta > 0.  After the loop, clipPoint is either -1 (don't do anything) or the
 *         clipping index in the read (from the end).
 *     </dd><br />
 *     <dt>Cycle based clipping</dt>
 *     <dd>Clips machine cycles from the read. Accepts a string of ranges of the form start1-end1,start2-end2, etc.
 *     For each start/end pair, removes bases in machine cycles from start to end, inclusive. These are 1-based values (positions).
 *     For example, 1-5,10-12 clips the first 5 bases, and then three bases at cycles 10, 11, and 12.
 *     </dd><br />
 *     <dt>Sequence matching</dt>
 *     <dd>Clips bases from that exactly match one of a number of base sequences. This employs an exact match algorithm,
 *     filtering only bases whose sequence exactly matches SEQ.</dd>
 * </dl>
 *
 *
 * <h3>Input</h3>
 * <p>
 *     Any number of BAM files.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 *     A new BAM file containing all of the reads from the input BAMs with the user-specified clipping
 *     operation applied to each read.
 * </p>
 * <p>
 *     <h4>Summary output (console)</h4>
 *     <pre>
 *     Number of examined reads              13
 *     Number of clipped reads               13
 *     Percent of clipped reads              100.00
 *     Number of examined bases              988
 *     Number of clipped bases               126
 *     Percent of clipped bases              12.75
 *     Number of quality-score clipped bases 126
 *     Number of range clipped bases         0
 *     Number of sequence clipped bases      0
 *     </pre>
 * </p>
 *
 * <h3>Example</h3>
 * <pre>
 *   java -jar GenomeAnalysisTK.jar \
 *     -T ClipReads \
 *     -R reference.fasta \
 *     -I original.bam \
 *     -o clipped.bam \
 *     -XF seqsToClip.fasta \
 *     -X CCCCC \
 *     -CT "1-5,11-15" \
 *     -QT 10
 * </pre>
 * <p>The command line shown above will apply all three arguments in combination. See the detailed examples below to see how the choice of clipping representation affects the output.</p>
 *
 *     <h4>Detailed clipping examples</h4>
 *     <p>Suppose we are given this read:</p>
 *     <pre>
 *     314KGAAXX090507:1:19:1420:1123#0        16      chrM    3116    29      76M     *       *       *
 *          TAGGACCCGGGCCCCCCTCCCCAATCCTCCAACGCATATAGCGGCCGCGCCTTCCCCCGTAAATGATATCATCTCA
 *          #################4?6/?2135;;;'1/=/<'B9;12;68?A79@,@==@9?=AAA3;A@B;A?B54;?ABA
 *     </pre>
 *
 *     <p>If we are clipping reads with -QT 10 and -CR WRITE_NS, we get:</p>
 *
 *     <pre>
 *     314KGAAXX090507:1:19:1420:1123#0        16      chrM    3116    29      76M     *       *       *
 *          NNNNNNNNNNNNNNNNNTCCCCAATCCTCCAACGCATATAGCGGCCGCGCCTTCCCCCGTAAATGATATCATCTCA
 *          #################4?6/?2135;;;'1/=/<'B9;12;68?A79@,@==@9?=AAA3;A@B;A?B54;?ABA
 *     </pre>
 *
 *     <p>Whereas with -QT 10 -CR WRITE_Q0S:</p>
 *     <pre>
 *     314KGAAXX090507:1:19:1420:1123#0        16      chrM    3116    29      76M     *       *       *
 *          TAGGACCCGGGCCCCCCTCCCCAATCCTCCAACGCATATAGCGGCCGCGCCTTCCCCCGTAAATGATATCATCTCA
 *          !!!!!!!!!!!!!!!!!4?6/?2135;;;'1/=/<'B9;12;68?A79@,@==@9?=AAA3;A@B;A?B54;?ABA
 *     </pre>
 *
 *     <p>Or -QT 10 -CR SOFTCLIP_BASES:</p>
 *     <pre>
 *     314KGAAXX090507:1:19:1420:1123#0        16      chrM    3133    29      17S59M  *       *       *
 *          TAGGACCCGGGCCCCCCTCCCCAATCCTCCAACGCATATAGCGGCCGCGCCTTCCCCCGTAAATGATATCATCTCA
 *          #################4?6/?2135;;;'1/=/<'B9;12;68?A79@,@==@9?=AAA3;A@B;A?B54;?ABA
 *     </pre>
 *
 */
@CommandLineProgramProperties(
        usage = "Read clipping based on quality, position or sequence matching.",
        usageShort = "Clip reads",
        programGroup = ReadProgramGroup.class
)
public final class ClipReads extends ReadWalker {

    private final Logger logger = LogManager.getLogger(ClipReads.class);

    /**
     * The output SAM/BAM file will be written here
     */
    @Argument(doc = "BAM output file", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    File OUTPUT;

    /**
     * If provided, ClipReads will write summary statistics about the clipping operations applied to the reads in this file.
     */
    @Argument(fullName = "outputStatistics", shortName = "os", doc = "File to output statistics", optional = true)
    File STATSOUTPUT = null;

    /**
     * If a value > 0 is provided, then the quality score based read clipper will be applied to the reads using this
     * quality score threshold.
     */
    @Argument(fullName = "qTrimmingThreshold", shortName = "QT", doc = "If provided, the Q-score clipper will be applied", optional = true)
    int qTrimmingThreshold = -1;

    /**
     * Clips machine cycles from the read. Accepts a string of ranges of the form start1-end1,start2-end2, etc.
     * For each start/end pair, removes bases in machine cycles from start to end, inclusive. These are 1-based
     * values (positions). For example, 1-5,10-12 clips the first 5 bases, and then three bases at cycles 10, 11,
     * and 12.
     */
    @Argument(fullName = "cyclesToTrim", shortName = "CT", doc = "String indicating machine cycles to clip from the reads", optional = true)
    String cyclesToClipArg = null;

    /**
     * Reads the sequences in the provided FASTA file, and clip any bases that exactly match any of the
     * sequences in the file.
     */
    @Argument(fullName = "clipSequencesFile", shortName = "XF", doc = "Remove sequences within reads matching the sequences in this FASTA file", optional = true)
    String clipSequenceFile = null;

    /**
     * Clips bases from the reads matching the provided SEQ.
     */
    @Argument(fullName = "clipSequence", shortName = "X", doc = "Remove sequences within reads matching this sequence", optional = true)
    List<String> clipSequencesArgs = null;

    /**
     * The different values for this argument determines how ClipReads applies clips to the reads.  This can range
     * from writing Ns over the clipped bases to hard clipping away the bases from the BAM.
     */
    @Argument(fullName = "clipRepresentation", shortName = "CR", doc = "How should we actually clip the bases?", optional = true)
    ClippingRepresentation clippingRepresentation = ClippingRepresentation.WRITE_NS;

    @Argument(fullName="read", doc="", optional = true)
    String onlyDoRead = null;

    /**
     * List of sequence that should be clipped from the reads
     */
    private List<SeqToClip> sequencesToClip = new ArrayList<SeqToClip>();

    /**
     * List of cycle start / stop pairs (0-based, stop is included in the cycle to remove) to clip from the reads
     */
    private List<Pair<Integer, Integer>> cyclesToClip = null;

    /**
     * Output reads is written to this BAM.
     */
    private SAMFileWriter outputBam;

    /**
     * Accumulator for the stats.
     */
    private ClippingData accumulator;

    /**
     * Output stream for the stats.
     */
    private PrintStream out;

    /**
     * The initialize function.
     */
    @Override
    public void onTraversalStart() {
        if (qTrimmingThreshold >= 0) {
            logger.info(String.format("Creating Q-score clipper with threshold %d", qTrimmingThreshold));
        }

        //
        // Initialize the sequences to clip
        //
        if (clipSequencesArgs != null) {
            int i = 0;
            for (String toClip : clipSequencesArgs) {
                i++;
                ReferenceSequence rs = new ReferenceSequence("CMDLINE-" + i, -1, StringUtil.stringToBytes(toClip));
                addSeqToClip(rs.getName(), rs.getBases());
            }
        }

        if (clipSequenceFile != null) {
            ReferenceSequenceFile rsf = ReferenceSequenceFileFactory.getReferenceSequenceFile(new File(clipSequenceFile));

            while (true) {
                ReferenceSequence rs = rsf.nextSequence();
                if (rs == null)
                    break;
                else {
                    addSeqToClip(rs.getName(), rs.getBases());
                }
            }
        }


        //
        // Initialize the cycle ranges to clip
        //
        if (cyclesToClipArg != null) {
            cyclesToClip = new ArrayList<Pair<Integer, Integer>>();
            for (String range : cyclesToClipArg.split(",")) {
                try {
                    String[] elts = range.split("-");
                    int start = Integer.parseInt(elts[0]) - 1;
                    int stop = Integer.parseInt(elts[1]) - 1;

                    if (start < 0) throw new Exception();
                    if (stop < start) throw new Exception();

                    logger.info(String.format("Creating cycle clipper %d-%d", start, stop));
                    cyclesToClip.add(new MutablePair<>(start, stop));
                } catch (Exception e) {
                    throw new RuntimeException("Badly formatted cyclesToClip argument: " + cyclesToClipArg);
                }
            }
        }

        final boolean presorted = EnumSet.of(ClippingRepresentation.WRITE_NS, ClippingRepresentation.WRITE_NS_Q0S, ClippingRepresentation.WRITE_Q0S).contains(clippingRepresentation);
        final SAMFileHeader outputHeader = getHeaderForReads().clone();
        outputBam = new SAMFileWriterFactory().makeWriter(outputHeader, presorted, OUTPUT, REFERENCE_FILE);

        accumulator = new ClippingData(sequencesToClip);
        try {
            out = new PrintStream(STATSOUTPUT);
        } catch (FileNotFoundException e) {
            throw new UserException.CouldNotCreateOutputFile(STATSOUTPUT, e);
        }
    }

    /**
     * The reads map function.
     *
     *
     * @param ref  the reference bases that correspond to our read, if a reference was provided
     * @param read the read itself, as a GATKSAMRecord
     */
    @Override
    public void apply(SAMRecord read, ReferenceContext ref) {
        if ( onlyDoRead == null || read.getReadName().equals(onlyDoRead) ) {
            if ( clippingRepresentation == ClippingRepresentation.HARDCLIP_BASES || clippingRepresentation == ClippingRepresentation.REVERT_SOFTCLIPPED_BASES )
                read = ReadClipper.revertSoftClippedBases(read);
            ReadClipperWithData clipper = new ReadClipperWithData(read, sequencesToClip);

            //
            // run all three clipping modules
            //
            clipBadQualityScores(clipper);
            clipCycles(clipper);
            clipSequences(clipper);
            accumulate(clipper);
        }
    }

    @Override
    public ClippingData onTraversalDone(){
        if ( out != null ){
           out.printf(accumulator.toString());
        }
        CloserUtil.close(outputBam);
        return accumulator;
    }
    /**
     * Helper function that adds a seq with name and bases (as bytes) to the list of sequences to be clipped
     *
     * @param name
     * @param bases
     */
    private void addSeqToClip(String name, byte[] bases) {
        SeqToClip clip = new SeqToClip(name, bases);
        sequencesToClip.add(clip);
        logger.info(String.format("Creating sequence clipper %s: %s/%s", clip.name, clip.seq, clip.revSeq));
    }


    /**
     * clip sequences from the reads that match all of the sequences in the global sequencesToClip variable.
     * Adds ClippingOps for each clip to clipper.
     *
     * @param clipper
     */
    private void clipSequences(ReadClipperWithData clipper) {
        if (sequencesToClip != null) {                // don't bother if we don't have any sequences to clip
            SAMRecord read = clipper.getRead();
            ClippingData data = clipper.getData();

            for (SeqToClip stc : sequencesToClip) {
                // we have a pattern for both the forward and the reverse strands
                Pattern pattern = read.getReadNegativeStrandFlag() ? stc.revPat : stc.fwdPat;
                String bases = read.getReadString();
                Matcher match = pattern.matcher(bases);

                // keep clipping until match.find() says it can't find anything else
                boolean found = true;   // go through at least once
                while (found) {
                    found = match.find();
                    //System.out.printf("Matching %s against %s/%s => %b%n", bases, stc.seq, stc.revSeq, found);
                    if (found) {
                        int start = match.start();
                        int stop = match.end() - 1;
                        //ClippingOp op = new ClippingOp(ClippingOp.ClippingType.MATCHES_CLIP_SEQ, start, stop, stc.seq);
                        ClippingOp op = new ClippingOp(start, stop);
                        clipper.addOp(op);
                        data.incSeqClippedBases(stc.seq, op.getLength());
                    }
                }
            }
            clipper.setData(data);
        }
    }

    /**
     * Convenence function that takes a read and the start / stop clipping positions based on the forward
     * strand, and returns start/stop values appropriate for the strand of the read.
     *
     * @param read
     * @param start
     * @param stop
     * @return
     */
    private Pair<Integer, Integer> strandAwarePositions(SAMRecord read, int start, int stop) {
        if (read.getReadNegativeStrandFlag())
            return new MutablePair<>(read.getReadLength() - stop - 1, read.getReadLength() - start - 1);
        else
            return new MutablePair<>(start, stop);
    }

    /**
     * clip bases at cycles between the ranges in cyclesToClip by adding appropriate ClippingOps to clipper.
     *
     * @param clipper
     */
    private void clipCycles(ReadClipperWithData clipper) {
        if (cyclesToClip != null) {
            SAMRecord read = clipper.getRead();
            ClippingData data = clipper.getData();

            for (Pair<Integer, Integer> p : cyclesToClip) {   // iterate over each cycle range
                int cycleStart = p.getLeft();
                int cycleStop = p.getRight();

                if (cycleStart < read.getReadLength()) {
                    // only try to clip if the cycleStart is less than the read's length
                    if (cycleStop >= read.getReadLength())
                        // we do tolerate [for convenience) clipping when the stop is beyond the end of the read
                        cycleStop = read.getReadLength() - 1;

                    Pair<Integer, Integer> startStop = strandAwarePositions(read, cycleStart, cycleStop);
                    int start = startStop.getLeft();
                    int stop = startStop.getRight();

                    //ClippingOp op = new ClippingOp(ClippingOp.ClippingType.WITHIN_CLIP_RANGE, start, stop, null);
                    ClippingOp op = new ClippingOp(start, stop);
                    clipper.addOp(op);
                    data.incNRangeClippedBases(op.getLength());
                }
            }
            clipper.setData(data);
        }
    }

    /**
     * Clip bases from the read in clipper from
     * <p/>
     * argmax_x{ \sum{i = x + 1}^l (qTrimmingThreshold - qual)
     * <p/>
     * to the end of the read.  This is blatantly stolen from BWA.
     * <p/>
     * Walk through the read from the end (in machine cycle order) to the beginning, calculating the
     * running sum of qTrimmingThreshold - qual.  While we do this, we track the maximum value of this
     * sum where the delta > 0.  After the loop, clipPoint is either -1 (don't do anything) or the
     * clipping index in the read (from the end).
     *
     * @param clipper
     */
    private void clipBadQualityScores(ReadClipperWithData clipper) {
        SAMRecord read = clipper.getRead();
        ClippingData data = clipper.getData();
        int readLen = read.getReadBases().length;
        byte[] quals = read.getBaseQualities();


        int clipSum = 0, lastMax = -1, clipPoint = -1; // -1 means no clip
        for (int i = readLen - 1; i >= 0; i--) {
            int baseIndex = read.getReadNegativeStrandFlag() ? readLen - i - 1 : i;
            byte qual = quals[baseIndex];
            clipSum += (qTrimmingThreshold - qual);
            if (clipSum >= 0 && (clipSum >= lastMax)) {
                lastMax = clipSum;
                clipPoint = baseIndex;
            }
        }

        if (clipPoint != -1) {
            int start = read.getReadNegativeStrandFlag() ? 0 : clipPoint;
            int stop = read.getReadNegativeStrandFlag() ? clipPoint : readLen - 1;
            //clipper.addOp(new ClippingOp(ClippingOp.ClippingType.LOW_Q_SCORES, start, stop, null));
            ClippingOp op = new ClippingOp(start, stop);
            clipper.addOp(op);
            data.incNQClippedBases(op.getLength());
        }
        clipper.setData(data);
    }

    private void accumulate(ReadClipperWithData clipper) {
        if ( clipper == null )
            return;

        SAMRecord clippedRead = clipper.clipRead(clippingRepresentation);
        if (outputBam != null) {
            outputBam.addAlignment(clippedRead);
        } else {
            out.println(clippedRead.getSAMString());
        }

        accumulator.nTotalReads++;
        accumulator.nTotalBases += clipper.getRead().getReadLength();
        if (clipper.wasClipped()) {
            accumulator.nClippedReads++;
            accumulator.addData(clipper.getData());
        }
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // utility classes
    //
    // --------------------------------------------------------------------------------------------------------------

    private static class SeqToClip {
        String name;
        String seq, revSeq;
        Pattern fwdPat, revPat;

        public SeqToClip(String name, byte[] bytez) {
            this.name = name;
            this.seq = new String(bytez);
            this.fwdPat = Pattern.compile(seq, Pattern.CASE_INSENSITIVE);
            this.revSeq = new String(BaseUtils.simpleReverseComplement(bytez));
            this.revPat = Pattern.compile(revSeq, Pattern.CASE_INSENSITIVE);
        }
    }

    public static class ClippingData {
        public long nTotalReads = 0;
        public long nTotalBases = 0;
        public long nClippedReads = 0;
        public long nClippedBases = 0;
        public long nQClippedBases = 0;
        public long nRangeClippedBases = 0;
        public long nSeqClippedBases = 0;

        HashMap<String, Long> seqClipCounts = new HashMap<>();

        public ClippingData(List<SeqToClip> clipSeqs) {
            for (SeqToClip clipSeq : clipSeqs) {
                seqClipCounts.put(clipSeq.seq, 0L);
            }
        }

        public void incNQClippedBases(int n) {
            nQClippedBases += n;
            nClippedBases += n;
        }

        public void incNRangeClippedBases(int n) {
            nRangeClippedBases += n;
            nClippedBases += n;
        }

        public void incSeqClippedBases(final String seq, int n) {
            nSeqClippedBases += n;
            nClippedBases += n;
            seqClipCounts.put(seq, seqClipCounts.get(seq) + n);
        }

        public void addData (ClippingData data) {
            nTotalReads += data.nTotalReads;
            nTotalBases += data.nTotalBases;
            nClippedReads += data.nClippedReads;
            nClippedBases += data.nClippedBases;
            nQClippedBases += data.nQClippedBases;
            nRangeClippedBases += data.nRangeClippedBases;
            nSeqClippedBases += data.nSeqClippedBases;

            for (String seqClip : data.seqClipCounts.keySet()) {
                Long count = data.seqClipCounts.get(seqClip);
                if (seqClipCounts.containsKey(seqClip))
                    count += seqClipCounts.get(seqClip);
                seqClipCounts.put(seqClip, count);
            }
        }

        public String toString() {
            StringBuilder s = new StringBuilder();

            s.append(Utils.dupString('-', 80) + "\n");
            s.append(String.format("Number of examined reads              %d%n", nTotalReads));
            s.append(String.format("Number of clipped reads               %d%n", nClippedReads));
            s.append(String.format("Percent of clipped reads              %.2f%n", (100.0 * nClippedReads) / nTotalReads));
            s.append(String.format("Number of examined bases              %d%n", nTotalBases));
            s.append(String.format("Number of clipped bases               %d%n", nClippedBases));
            s.append(String.format("Percent of clipped bases              %.2f%n", (100.0 * nClippedBases) / nTotalBases));
            s.append(String.format("Number of quality-score clipped bases %d%n", nQClippedBases));
            s.append(String.format("Number of range clipped bases         %d%n", nRangeClippedBases));
            s.append(String.format("Number of sequence clipped bases      %d%n", nSeqClippedBases));

            for (Map.Entry<String, Long> elt : seqClipCounts.entrySet()) {
                s.append(String.format("  %8d clip sites matching %s%n", elt.getValue(), elt.getKey()));
            }

            s.append(Utils.dupString('-', 80) + "\n");
            return s.toString();
        }
    }

    public static class ReadClipperWithData extends ReadClipper {
        private ClippingData data;

        public ReadClipperWithData(SAMRecord read, List<SeqToClip> clipSeqs) {
            super(read);
            data = new ClippingData(clipSeqs);
        }

        public ClippingData getData() {
            return data;
        }

        public void setData(ClippingData data) {
            this.data = data;
        }

        public void addData(ClippingData data) {
            this.data.addData(data);
        }
    }


}
