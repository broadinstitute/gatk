package org.broadinstitute.hellbender.tools.walkers;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.UnalignedRead;
import org.broadinstitute.hellbender.utils.read.ByteSequence;
import picard.cmdline.programgroups.VariantEvaluationProgramGroup;

import java.io.*;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;

@CommandLineProgramProperties(
        summary = "Emit FASTQ for reads that overlap a 2Kb window around each SV, fixing up " +
                "single-base insertions and deletions that appear in a single read.",
        oneLineSummary = "Perform HiFi error correction in the vicinity of an SV.",
        programGroup = VariantEvaluationProgramGroup.class
)
public final class ErrorCorrectHiFi extends VariantWalker {
    public static final int WINDOW_SIZE = 2500;
    public static final int HUGE_EVENT_SIZE = 7500;
    public static final ByteSequence BYTE_SEQUENCE_A = new ByteSequence((byte)'A');
    public static final String SCRIPT_TEXT =
"#!/bin/sh\n" +
"minimap2 -axmap-hifi %s reads.fq | samtools sort -OBAM - > align.bam &&\n" +
"    samtools index align.bam &&\n" +
"    igv -g hg38 -n corrected,%s -l %s align.bam %s\n";

    @Argument(fullName = "mmi-index", shortName = "M", doc = "A minimap2 reference index")
    public File mmiIndex;

    final File tmpDir = IOUtils.createTempDir("echf");
    final File outputFASTQ = new File(tmpDir, "reads.fq");
    final File scriptFile = new File(tmpDir, "run.sh");
    final ProcessBuilder scriptRunner = new ProcessBuilder(scriptFile.getAbsolutePath());

    @Override
    public boolean requiresReads() { return true; }

    @Override
    public void apply( final VariantContext variant,
                       final ReadsContext readsContext,
                       final ReferenceContext refContext,
                       final FeatureContext featureContext ) {
        if ( variant.hasAttribute(GATKSVVCFConstants.SVTYPE) ) {
            final String contig = variant.getContig();
            final int eventStart = variant.getStart();
            final String svType = variant.getAttributeAsString(GATKSVVCFConstants.SVTYPE, null);
            final int svLen = variant.getAttributeAsInt(GATKSVVCFConstants.SVLEN, 0);
            final String seqName = contig + "_" + eventStart + "_" + svType + "_" + Math.abs(svLen);
            final int eventEnd = variant.getEnd();
            final int windowStart = Math.max(1, eventStart - WINDOW_SIZE + 1);
            final int windowEnd = eventEnd + WINDOW_SIZE;
            if ( eventEnd - eventStart < HUGE_EVENT_SIZE ) {
                final SimpleInterval windowInterval =
                        new SimpleInterval(contig, windowStart, windowEnd);
                processWindow(readsContext.iterator(windowInterval));
                writeScript(seqName, windowInterval);
                alignAndDisplay();
            } else {
                final SimpleInterval leadInterval =
                        new SimpleInterval(contig, windowStart, eventStart + WINDOW_SIZE);
                processWindow(readsContext.iterator(leadInterval));
                writeScript(seqName + "_BND1", leadInterval);
                alignAndDisplay();

                final SimpleInterval lagInterval =
                        new SimpleInterval(contig, eventEnd - WINDOW_SIZE, windowEnd);
                processWindow(readsContext.iterator(lagInterval));
                writeScript(seqName + "_BND2", lagInterval);
                alignAndDisplay();
            }
        }
    }

    private void processWindow( final Iterator<GATKRead> readsIterator ) {
        final List<GATKRead> reads = new ArrayList<>();
        int windowStartPos = Integer.MAX_VALUE;
        int windowEndPos = 0;
        while ( readsIterator.hasNext() ) {
            final GATKRead read = readsIterator.next();
            reads.add(read);
            windowStartPos = Math.min(windowStartPos, read.getStart());
            windowEndPos = Math.max(windowEndPos, read.getEnd());
        }

        final List<CallIterator> callIteratorList =
                errorCorrectWindow(reads, windowStartPos, windowEndPos);

        try ( final BufferedWriter writer = new BufferedWriter(new FileWriter(outputFASTQ)) ) {
            for ( final CallIterator callIterator : callIteratorList ) {
                callIterator.getUnalignedRead().writeFASTQ(writer);
            }
        } catch ( final IOException ioe ) {
            throw new GATKException("Can't write output to " + outputFASTQ, ioe);
        }
    }

    public static List<CallIterator> errorCorrectWindow( final List<GATKRead> reads,
                                                         final int windowStartPos,
                                                         final int windowEndPos ) {
        final List<CallIterator> callIteratorList = new ArrayList<>(reads.size());
        for ( final GATKRead read : reads ) {
                callIteratorList.add(new CallIterator(read, windowStartPos));
        }
        for ( int refPos = windowStartPos; refPos <= windowEndPos && !callIteratorList.isEmpty(); ++refPos ) {
            for ( CallIterator callIterator : callIteratorList ) {
                if ( !callIterator.hasNext() ) {
                    continue;
                }
                final Call call = callIterator.next();
                if ( call == null ) {
                    continue;
                }
                final ByteSequence calls = call.getCalls();
                final int len = calls.length();
                // if we have a single-base deletion, or a single-base insertion
                if ( (len == 0 && call.getCigarElement().getLength() == 1) || len == 2 ) {
                    callIterator.fixup(call, BYTE_SEQUENCE_A);
                }
            }
        }
        return callIteratorList;
    }

    private void writeScript( final String seqName, final SimpleInterval interval ) {
        final String resolvedScript =
                String.format(SCRIPT_TEXT,
                        mmiIndex.getAbsolutePath(),
                        seqName,
                        interval,
                        readArguments.getReadPaths().get(0).toAbsolutePath());
        IOUtils.writeByteArrayToFile(resolvedScript.getBytes(), scriptFile);
        if ( !scriptFile.setExecutable(true) ) {
            throw new UserException("Can't make minimap2/igv script executable.");
        }
        scriptRunner.directory(tmpDir).inheritIO();
    }

    private void alignAndDisplay() {
        try {
            final Process process = scriptRunner.start();
            final int exitValue = process.waitFor();
            if ( exitValue != 0 ) {
                throw new UserException("Script that runs minimap2 and igv returned " + exitValue);
            }
        } catch ( final IOException | InterruptedException ex ) {
            throw new UserException("Script that runs minimap2 and igv failed.", ex);
        }
    }

    /**
     * What's happening at some reference position in some read.
     * Note that inserts are weird:  the calls/quals byte sequence will be one longer than the cigar
     * element's length, because it has the call at that reference position as well as the inserted
     * calls.  (It works like a vcf entry that has, for example, a ref allele of A, and an alt
     * allele of ATC.  That's a 2I, even though the alt allele is 3 base calls in length.)
     * Deletions have 0-length byte sequences regardless of the length of the cigar element.
     */
    public static final class Call {
        private final ByteSequence calls;
        private final ByteSequence quals;
        private final CigarElement cigarElement;

        public Call( final ByteSequence calls,
                     final ByteSequence quals,
                     final CigarElement cigarElement ) {
            this.calls = calls;
            this.quals = quals;
            this.cigarElement = cigarElement;
        }

        public ByteSequence getCalls() { return calls; }
        public ByteSequence getQuals() { return quals; }

        public int getMeanQual() {
            int qSum = 0;
            final int qLen = quals.length();
            for ( int idx = 0; idx != qLen; ++idx ) {
                qSum += quals.byteAt(idx);
            }
            return Math.round(1.F*qSum/qLen);
        }
        public CigarElement getCigarElement() { return cigarElement; }

        @Override
        public int hashCode() {
            if ( calls.length() != 0 ) {
                return calls.hashCode();
            }
            return 47 * cigarElement.getLength();
        }

        @Override
        public boolean equals( final Object obj ) {
            if ( !(obj instanceof final Call that) ) {
                return false;
            }
            if ( calls.length() != 0 ) {
                return calls.equals(that.calls);
            }
            return that.calls.length() == 0 && cigarElement.getLength() == that.cigarElement.getLength();
        }

        @Override
        public String toString() {
            if ( calls.length() == 1 ) return calls.toString() + getMeanQual();
            return cigarElement.toString() + getMeanQual();
        }
    }

    /**
     *  Instructions to remedy a single-base insertion (value will be 0)
     *  or a single-base deletion (value will be A, C, G, or T) in some read.
     */
    public static final class Fixup {
        private final int readIndex;
        private final byte value;

        public Fixup( final int readIndex, final byte value ) {
            this.readIndex = readIndex;
            this.value = value;
        }
        public int getReadIndex() { return readIndex; }
        public byte getValue() { return value; }
    }

    /**
     * Iterate over a read to yield a call for each reference position covered by the read.
     */
    public static final class CallIterator implements Iterator<Call> {
        private final String readName;
        private final ByteSequence calls;
        private final ByteSequence quals;
        private final List<CigarElement> cigarElements;
        private final List<Fixup> fixups;
        private final int initialReadIndex;
        private int cigarElementsIndex;
        private int cigarElementIndex;
        private int readIndex;
        private int fixupLengthDelta;

        public CallIterator( final GATKRead read, final int initialRefPos ) {
            readName = read.getName();
            calls = new ByteSequence(read.getBasesNoCopy());
            quals = new ByteSequence(read.getBaseQualitiesNoCopy());
            final List<CigarElement> elements = read.getCigarElements();
            cigarElements = new ArrayList<>(elements.size());
            fixups = new ArrayList<>();
            cigarElementsIndex = 0;
            cigarElementIndex = 0;
            readIndex = 0;
            fixupLengthDelta = 0;
            advance(elements, read.getStart(), initialRefPos);
            initialReadIndex = readIndex;
        }

        public boolean hasNext() {
            return cigarElementsIndex < cigarElements.size();
        }

        public Call next() {
            if ( !hasNext() ) {
                throw new NoSuchElementException();
            }
            final CigarElement element = cigarElements.get(cigarElementsIndex);
            final CigarOperator op = element.getOperator();
            final Call call;
            if ( op == CigarOperator.M ) {
                call = new Call(calls.subSequence(readIndex, readIndex + 1),
                                quals.subSequence(readIndex, readIndex + 1),
                                element);
                readIndex += 1;
            } else if ( op == CigarOperator.D ) {
                if ( cigarElementIndex != 0 ) {
                    call = null;
                } else {
                    call = new Call(calls.subSequence(readIndex, readIndex),
                            quals.subSequence(readIndex - 1, readIndex + 1),
                            element);
                }
            } else if ( op == CigarOperator.I ) {
                final int len = element.getLength();
                call = new Call(calls.subSequence(readIndex, readIndex + len),
                                quals.subSequence(readIndex, readIndex + len),
                                element);
                readIndex += len;
                cigarElementIndex = len - 1;
            } else if ( op == CigarOperator.SKIPPED_REGION ) {
                call = null;
            } else {
                throw new GATKException("Unexpected operator.");
            }
            if ( ++cigarElementIndex == element.getLength() ) {
                cigarElementsIndex += 1;
                cigarElementIndex = 0;
            }
            return call;
        }

        public UnalignedRead getUnalignedRead() {
            if ( fixups.isEmpty() ) {
                return new UnalignedRead(readName, calls, quals);
            }
            final int newLength = readIndex - initialReadIndex + fixupLengthDelta;
            final byte[] callBytes = new byte[newLength];
            final byte[] qualBytes = new byte[newLength];
            final Iterator<Fixup> fixupIterator = fixups.iterator();
            Fixup curFixup = fixupIterator.hasNext() ? fixupIterator.next() : null;
            int destIdx = 0;
            for ( int idx = initialReadIndex; idx != readIndex; ++idx ) {
                if ( curFixup != null && curFixup.getReadIndex() == idx ) {
                    if ( curFixup.getValue() != 0 ) {
                        callBytes[destIdx] = calls.byteAt(idx);
                        qualBytes[destIdx] = quals.byteAt(idx);
                        destIdx += 1;
                        callBytes[destIdx] = curFixup.getValue();
                        qualBytes[destIdx] = 0;
                        destIdx += 1;
                    }
                    curFixup = fixupIterator.hasNext() ? fixupIterator.next() : null;
                } else {
                    callBytes[destIdx] = calls.byteAt(idx);
                    qualBytes[destIdx] = quals.byteAt(idx);
                    destIdx += 1;
                }
            }
            if ( destIdx != callBytes.length ) {
                throw new GATKException("failed call fixup");
            }
            return new UnalignedRead(readName, new ByteSequence(callBytes), new ByteSequence(qualBytes));
        }

        public void fixup( final Call call, final ByteSequence alignedValue ) {
            final ByteSequence calls = call.getCalls();
            if ( calls.length() == 0 ) {  // i.e., we're fixing a single-base deletion
                fixups.add(new Fixup(calls.getStart(), alignedValue.byteAt(0)));
                fixupLengthDelta += 1;
            } else if ( calls.length() == 2 ) { // i.e., we're fixing a single-base insertion
                fixups.add(new Fixup(calls.getStart()+1, (byte)0));
                fixupLengthDelta -= 1;
            } else {
                throw new GATKException("unexpected fixup length");
            }
        }

        private void advance( final List<CigarElement> elements,
                              final int curRefPos,
                              final int initialRefPos ) {
            int distanceToTarget = initialRefPos - curRefPos;
            if ( distanceToTarget < 0 ) {
                cigarElements.add(new CigarElement(-distanceToTarget, CigarOperator.N));
            }
            for ( final CigarElement element : elements ) {
                final int elementLength = element.getLength();
                final CigarOperator op = element.getOperator();
                if ( distanceToTarget >= 0 ) {
                    if ( op.consumesReferenceBases() ) {
                        if ( elementLength > distanceToTarget ) {
                            cigarElementIndex = distanceToTarget;
                        }
                        distanceToTarget -= elementLength;
                    }
                    if ( op.consumesReadBases() ) {
                        if ( distanceToTarget < 0 ) {
                            readIndex += cigarElementIndex;
                        } else {
                            readIndex += elementLength;
                        }
                    }
                }
                if ( distanceToTarget < 0 ) { // not else!  previous if may have made distanceToTarget < 0.
                    switch ( op ) {
                        case I:
                            cigarElements.add(handleInsertion(element));
                            break;
                        case M: case D:
                            cigarElements.add(element);
                            break;
                        case EQ: case X:
                            cigarElements.add(new CigarElement(elementLength, CigarOperator.M));
                            break;
                        case N:
                            cigarElements.add(new CigarElement(elementLength, CigarOperator.D));
                            break;
                        case S: case H: case P:
                            break;
                    }
                }
            }
        }

        // include final aligned base from previous call with the insertion that follows
        private CigarElement handleInsertion( final CigarElement cigarElement ) {
            final int lastIdx = cigarElements.size() - 1;
            final CigarElement previousElement = cigarElements.get(lastIdx);
            final int len = previousElement.getLength();
            if ( previousElement.getOperator() == CigarOperator.M ) {
                if ( len == 1 ) {
                    cigarElements.remove(lastIdx);
                } else {
                    cigarElements.set(lastIdx, new CigarElement(len - 1, CigarOperator.M));
                    if ( lastIdx == 0 && cigarElementIndex == len - 1 ) {
                        cigarElementIndex = 0;
                        cigarElementsIndex = 1;
                    }
                }
                return new CigarElement(cigarElement.getLength() + 1, CigarOperator.I);
            }
            if ( previousElement.getOperator() == CigarOperator.D ) {
                if ( len == 1 ) {
                    cigarElements.remove(lastIdx);
                } else {
                    cigarElements.set(lastIdx, new CigarElement(len - 1, CigarOperator.D));
                }
                return cigarElement;
            }
            throw new UserException("Found I in cigar with no preceding D or M");
        }
    }
}
