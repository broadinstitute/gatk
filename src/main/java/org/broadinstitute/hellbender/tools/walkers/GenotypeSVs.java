package org.broadinstitute.hellbender.tools.walkers;

import edu.mit.broad.tedsUtils.align.Aligner;
import edu.mit.broad.tedsUtils.align.Scorer;
import edu.mit.broad.tedsUtils.align.SomewhatGlobalAligner;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.ByteSequence;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SequenceRC;
import picard.cmdline.programgroups.VariantEvaluationProgramGroup;

import java.io.*;
import java.util.*;
import java.util.regex.Pattern;

@CommandLineProgramProperties(
        summary = "For each SV, write a couple of snippets of sequence to look for.",
        oneLineSummary = "For each SV, write a couple of snippets of sequence to look for.",
        programGroup = VariantEvaluationProgramGroup.class
)
public final class GenotypeSVs extends VariantWalker {
    public static final int MINIMUM_SV_SIZE = 50;
    public static final int BIG_EVENT_LENGTH = 500;
    public static final int WINDOW_SIZE = 100;
    public static final int REF_LEADIN_LEN = 12;
    public static final int MAX_READ_LEN = REF_LEADIN_LEN + WINDOW_SIZE;
    public static final int MIN_READ_LEN = REF_LEADIN_LEN + 25;
    public static final int MIN_CALLS_TO_GENOTYPE = 7;
    public static final float MIN_FRACTION_FOR_ALLELE_CALL = .15f;

    // ]p]t -- the BND locus, p, is the end of the reference sequence that precedes t
    public static final Pattern BND_PRE_FWD = Pattern.compile("]([!-~]*):([1-9][0-9]*)]([acgtnACGTN]+)");

    // [p[t -- the BND locus, p, is the start of the reference sequence the reverse complement of which precedes t
    public static final Pattern BND_PRE_REV = Pattern.compile("\\[([!-~]*):([1-9][0-9]*)\\[([acgtnACGTN]+)");

    // t[p[ -- the BND locus, p, is the start of the reference sequence that follows t
    public static final Pattern BND_POST_FWD = Pattern.compile("([acgtnACGTN]+)\\[([!-~]*):([1-9][0-9]*)\\[");

    // t]p] -- the BND locus, p, is the end of the reference sequence the reverse complement of which follows t
    public static final Pattern BND_POST_REV = Pattern.compile("([acgtnACGTN]+)]([!-~]*):([1-9][0-9]*)]");

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="Write output to this file")
    public File outputFile;
    public BufferedWriter writer;

    @Override
    public void onTraversalStart() {
        super.onTraversalStart();
        try {
            writer = new BufferedWriter(new FileWriter(outputFile));
        } catch ( final IOException ioe ) {
            throw new UserException("Can't open output file " + outputFile, ioe);
        }
    }

    @Override
    public Object onTraversalSuccess() {
        try {
            writer.close();
        } catch ( final IOException ioe ) {
            throw new UserException("Can't close output file " + outputFile, ioe);
        }
        return super.onTraversalSuccess();
    }

    @Override
    public void apply( final VariantContext vc,
                       final ReadsContext readsContext,
                       final ReferenceContext refContext,
                       final FeatureContext featureContext) {
        final String svType;
        final int svLen;
        if ( vc.hasAttribute(GATKSVVCFConstants.SVTYPE) ) {
            svType = vc.getAttributeAsString(GATKSVVCFConstants.SVTYPE, null);
            svLen = Math.abs(vc.getAttributeAsInt(GATKSVVCFConstants.SVLEN, 0));
        } else if ( vc.isIndel() ) {
            final int refLen = vc.getReference().getBaseString().length();
            final int altLen = vc.getAlternateAllele(0).getBaseString().length();
            if ( refLen >= altLen ) {
                svType = GATKSVVCFConstants.SYMB_ALT_STRING_DEL;
                svLen = refLen - altLen;
            } else {
                svType = GATKSVVCFConstants.SYMB_ALT_STRING_INS;
                svLen = altLen - refLen;
            }
        } else {
            svType = null;
            svLen = 0;
        }
        if ( svType == null || svLen < MINIMUM_SV_SIZE ) {
            return;
        }

        switch ( svType ) {
            case GATKSVVCFConstants.SYMB_ALT_STRING_DEL -> testDeletion(vc, refContext, readsContext);
            case GATKSVVCFConstants.SYMB_ALT_STRING_INS -> testInsertion(vc, refContext, readsContext);
            case GATKSVVCFConstants.SYMB_ALT_STRING_DUP -> testDuplication(vc, refContext, readsContext);
            case GATKSVVCFConstants.SYMB_ALT_STRING_INV -> testInversion(vc, refContext, readsContext);
            case GATKSVVCFConstants.BREAKEND_STR -> testBreakend(vc, refContext, readsContext);
/*
            case GATKSVVCFConstants.SYMB_ALT_STRING_DUP -> {
                final SimpleInterval dupInterval = new SimpleInterval(contig, start + 1, end);
                final String altBases = new String(refContext.getBases(dupInterval));
                writeRef(seqName + "1", leadBases + altBases.substring(1, WINDOW_SIZE + 1));
                writeRef(seqName + "2", altBases.substring(altBases.length() - WINDOW_SIZE) + lagBases);
            }
            case GATKSVVCFConstants.SYMB_ALT_STRING_INV -> {
                final SimpleInterval invInterval = new SimpleInterval(contig, start + 1, end);
                final byte[] invBaseCalls = refContext.getBases(invInterval);
                SequenceUtil.reverseComplement(invBaseCalls);
                final String altBases = new String(invBaseCalls);
                writeRef(seqName + "1", leadBases + altBases.substring(1, WINDOW_SIZE + 1));
                writeRef(seqName + "2", altBases.substring(altBases.length() - WINDOW_SIZE) + lagBases);
            }
            case GATKSVVCFConstants.BREAKEND_STR -> {
                final String bndDescription = variant.getAlternateAllele(0).getDisplayString();
                Matcher matcher;
                if ( (matcher = BND_PRE_FWD.matcher(bndDescription)).matches() ) {
                    final String tig = matcher.group(1);
                    final int bndPosEnd = Integer.parseInt(matcher.group(2));
                    final int bndPosBeg = Math.max(1, bndPosEnd - WINDOW_SIZE + 1);
                    final SimpleInterval bndInterval = new SimpleInterval(tig, bndPosBeg, bndPosEnd);
                    final String altBases =
                            new String(refContext.getBases(bndInterval)) + matcher.group(3);
                    writeRef(seqName, altBases.substring(altBases.length() - WINDOW_SIZE) + lagBases);
                } else if ( (matcher = BND_PRE_REV.matcher(bndDescription)).matches() ) {
                    final String tig = matcher.group(1);
                    final int bndPosBeg = Integer.parseInt(matcher.group(2));
                    final int bndPosEnd = bndPosBeg + WINDOW_SIZE - 1;
                    final SimpleInterval bndInterval = new SimpleInterval(tig, bndPosBeg, bndPosEnd);
                    final byte[] bndBaseCalls = refContext.getBases(bndInterval);
                    SequenceUtil.reverseComplement(bndBaseCalls);
                    final String altBases = new String(bndBaseCalls) + matcher.group(3);
                    writeRef(seqName, altBases.substring(altBases.length() - WINDOW_SIZE) + lagBases);
                } else if ( (matcher = BND_POST_FWD.matcher(bndDescription)).matches() ) {
                    final String tig = matcher.group(2);
                    final int bndPosBeg = Integer.parseInt(matcher.group(3));
                    final int bndPosEnd = bndPosBeg + WINDOW_SIZE - 1;
                    final SimpleInterval bndInterval = new SimpleInterval(tig, bndPosBeg, bndPosEnd);
                    final String altBases =
                            matcher.group(1) + new String(refContext.getBases(bndInterval));
                    writeRef(seqName, leadBases + altBases.substring(0, WINDOW_SIZE));
                } else if ( (matcher = BND_POST_REV.matcher(bndDescription)).matches() ) {
                    final String tig = matcher.group(2);
                    final int bndPosEnd = Integer.parseInt(matcher.group(3));
                    final int bndPosBeg = Math.max(1, bndPosEnd - WINDOW_SIZE + 1);
                    final SimpleInterval bndInterval = new SimpleInterval(tig, bndPosBeg, bndPosEnd);
                    final byte[] bndBaseCalls = refContext.getBases(bndInterval);
                    SequenceUtil.reverseComplement(bndBaseCalls);
                    final String altBases = matcher.group(1) + new String(bndBaseCalls);
                    writeRef(seqName, leadBases + altBases.substring(0, WINDOW_SIZE));
                } else {
                    logger.warn("Can't interpret the BND description: " + bndDescription);
                }
            }
*/
            default -> logger.warn("Don't know how to handle SVTYPE=" + svType +
                                    " for variant " + vc.getID() +
                                    " at " + vc.getContig() + ":" + vc.getStart());
        }
    }

    private void testDeletion( final VariantContext vc,
                               final ReferenceContext refContext,
                               final ReadsContext readsContext ) {
        final int start = Math.max(1, vc.getStart() + 1 - REF_LEADIN_LEN);
        final String contig = vc.getContig();
        final ByteSequence refCalls = getRefCalls(contig, start, MAX_READ_LEN, refContext);

        final int altStart = vc.getEnd() + 1;
        final ByteSequence altWindow = getRefCalls(contig, altStart, WINDOW_SIZE, refContext);
        final ByteSequence altCalls = refCalls.subSequence(0, REF_LEADIN_LEN).append(altWindow);

        final SimpleInterval alignInterval = new SimpleInterval(contig, start, altStart);
        final ScoreSummary scoreSummary = new ScoreSummary();
        alignReads(scoreSummary, refCalls, altCalls, alignInterval, readsContext, false);
        try {
            scoreSummary.writeSummary(String.format("%s:%09d", contig, vc.getStart()), writer);
        } catch ( final IOException ioe ) {
            throw new UserException("Can't write data for " + vc.getID() +
                        " at " + vc.getContig() + ":" + vc.getStart() + " to " + outputFile, ioe);
        }
    }

    private void testInsertion( final VariantContext vc,
                                final ReferenceContext refContext,
                                final ReadsContext readsContext ) {
        final int start = Math.max(1, vc.getStart() + 1 - REF_LEADIN_LEN);
        final String contig = vc.getContig();
        final SimpleInterval alignInterval = new SimpleInterval(contig, start, start);
        final ScoreSummary scoreSummary = new ScoreSummary();

        final ByteSequence refCalls = getRefCalls(contig, start, MAX_READ_LEN, refContext);

        final ByteSequence altAlleleCalls = new ByteSequence(vc.getAlternateAllele(0).getBases());
        final int altLength = altAlleleCalls.length() - 1;
        if ( altLength >= WINDOW_SIZE ) {
            final ByteSequence altCalls =
                    refCalls.subSequence(0, REF_LEADIN_LEN)
                            .append(altAlleleCalls.subSequence(1, WINDOW_SIZE + 1));
            alignReads(scoreSummary, refCalls, altCalls, alignInterval, readsContext, false);
            final ByteSequence altCalls2 =
                    altAlleleCalls.subSequence(altLength - WINDOW_SIZE + 1, altLength + 1)
                            .append(refCalls.subSequence(REF_LEADIN_LEN, REF_LEADIN_LEN + REF_LEADIN_LEN));
            alignReads(scoreSummary, refCalls, altCalls2, alignInterval, readsContext, true);
        } else {
            final ByteSequence altCalls =
                    refCalls.subSequence(0, REF_LEADIN_LEN)
                            .append(altAlleleCalls.subSequence(1))
                            .append(refCalls.subSequence(REF_LEADIN_LEN, REF_LEADIN_LEN + WINDOW_SIZE - altLength));
            alignReads(scoreSummary, refCalls, altCalls, alignInterval, readsContext, false);
        }

        try {
            scoreSummary.writeSummary(String.format("%s:%09d", contig, vc.getStart()), writer);
        } catch ( final IOException ioe ) {
            throw new UserException("Can't write data for " + vc.getID() +
                    " at " + vc.getContig() + ":" + vc.getStart() + " to " + outputFile, ioe);
        }
    }

    private void testDuplication( final VariantContext vc,
                                  final ReferenceContext refContext,
                                  final ReadsContext readsContext ) {
    }

    private void testInversion( final VariantContext vc,
                                final ReferenceContext refContext,
                                final ReadsContext readsContext ) {
    }

    private void testBreakend( final VariantContext vc,
                               final ReferenceContext refContext,
                               final ReadsContext readsContext ) {
    }

    public static ByteSequence getRefCalls( final String contig, final int start, final int length,
                                      final ReferenceContext refContext ) {
        final SimpleInterval refInterval = new SimpleInterval(contig, start, start + length - 1);
        return new ByteSequence(refContext.getBases(refInterval));
    }

    private static int getMinReadLen( final ByteSequence refCalls, final ByteSequence altCalls ) {
        int minReadLen = 0;
        while ( minReadLen < MAX_READ_LEN ) {
            if ( refCalls.charAt(minReadLen) != altCalls.charAt(minReadLen) ) {
                break;
            }
            minReadLen += 1;
        }
        minReadLen += 1; // bump index of 1st disagreement to make a length out of it
        if ( minReadLen < MIN_READ_LEN ) {
            minReadLen = MIN_READ_LEN;
        }
        return minReadLen;
    }

    private static void alignReads( final ScoreSummary scoreSummary,
                                    final ByteSequence refCalls,
                                    final ByteSequence altCalls,
                                    final SimpleInterval alignInterval,
                                    final ReadsContext readsContext,
                                    final boolean alignRC ) {
        final int start = alignInterval.getStart();
        final int end = alignInterval.getEnd();
        final List<SimpleInterval> intervalsToQuery;
        if ( end - start < BIG_EVENT_LENGTH ) {
            intervalsToQuery = Collections.singletonList(alignInterval);
        } else {
            intervalsToQuery = new ArrayList<>(2);
            final String contig = alignInterval.getContig();
            intervalsToQuery.add(new SimpleInterval(contig, start, start));
            intervalsToQuery.add(new SimpleInterval(contig, end, end));
        }
        final int minReadLen = getMinReadLen(refCalls, altCalls);
        for ( final SimpleInterval interval : intervalsToQuery ) {
            final Iterator<GATKRead> readItr = readsContext.iterator(interval);
            while ( readItr.hasNext() ) {
                final GATKRead read = readItr.next();
                if ( minReadLen > MAX_READ_LEN ) {
                    scoreSummary.score(0, 0);
                    continue;
                }
                int readStartPos = getReadPositionForRefPosition(read, start);
                // TODO: evaluate validity of this condition for all event types
                if ( readStartPos == -1 && start != end ) {
                    readStartPos = getReadPositionForRefPosition(read, end) - REF_LEADIN_LEN;
                }
                if ( readStartPos < 0 ) {
                    continue;
                }
                final ByteSequence allReadCalls = new ByteSequence(read.getBasesNoCopy());
                final int readLen = Math.min(allReadCalls.length() - readStartPos, MAX_READ_LEN);
                if ( readLen >= minReadLen ) {
                    final ByteSequence readCalls = allReadCalls.subSequence(readStartPos, readStartPos + readLen);
                    final Scorer scorer = new Scorer();
                    if ( alignRC ) {
                        final CharSequence readCallsRC = new SequenceRC(readCalls);
                        scoreSummary.score(
                                Math.round(new SomewhatGlobalAligner(readCallsRC, new SequenceRC(refCalls), scorer, minReadLen).getScore().getScore()),
                                Math.round(new SomewhatGlobalAligner(readCallsRC, new SequenceRC(altCalls), scorer, minReadLen).getScore().getScore()));
                    } else {
                        scoreSummary.score(
                                Math.round(new SomewhatGlobalAligner(readCalls, refCalls, scorer, minReadLen).getScore().getScore()),
                                Math.round(new SomewhatGlobalAligner(readCalls, altCalls, scorer, minReadLen).getScore().getScore()));
                    }
                }
            }
        }
    }

    public static int getReadPositionForRefPosition( final GATKRead read, final int refPosition ) {
        int curRefPosition = read.getStart();
        final List<CigarElement> cigarElements = read.getCigarElements();
        if ( curRefPosition > refPosition ) {
            final CigarElement firstElement = cigarElements.get(0);
            if ( firstElement.getOperator() == CigarOperator.S ) {
                final int virtualStart = firstElement.getLength() - (curRefPosition - refPosition);
                if ( virtualStart >= 0 ) {
                    return virtualStart;
                }
            }
            return -1;
        }
        int curReadPosition = 0;
        for ( final CigarElement cigarElement : cigarElements ) {
            final int opLen = cigarElement.getLength();
            final CigarOperator op = cigarElement.getOperator();
            if ( op.consumesReferenceBases() ) {
                final int nextRefPosition = curRefPosition + opLen;
                if ( nextRefPosition > refPosition ) {
                    if ( op.consumesReadBases() ) {
                        curReadPosition += opLen - (nextRefPosition - refPosition);
                    }
                    return curReadPosition;
                }
                curRefPosition = nextRefPosition;
            }
            if ( op.consumesReadBases() ) {
                curReadPosition += opLen;
            }
        }
        return -1;
    }

//    public record ScorePair( int refScore, int altScore ) {}

    public final static class ScoreSummary {
//        final Map<ScorePair, Integer> scores =
//                new TreeMap<>(Comparator.comparing(ScorePair::refScore).thenComparing(ScorePair::altScore));
        int refBetter = 0;
        int same = 0;
        int altBetter = 0;

        public void score( final int refScore, final int altScore ) {
//            scores.merge(new ScorePair(refScore, altScore), 1, Integer::sum);
            if ( refScore > altScore ) {
                refBetter += 1;
            } else if ( altScore > refScore ) {
                altBetter += 1;
            } else {
                same += 1;
            }
        }

        public void writeSummary( final String eventID,
                                  final BufferedWriter writer ) throws IOException {
//            writer.write('#');
            writer.write(eventID);
            writer.write('\t');
            final int totalCounts = refBetter + same + altBetter;
            if ( totalCounts < MIN_CALLS_TO_GENOTYPE ) {
                writer.write("./.");
            } else if ( 1.0f * refBetter / totalCounts > MIN_FRACTION_FOR_ALLELE_CALL ) {
                if ( 1.0f * altBetter / totalCounts > MIN_FRACTION_FOR_ALLELE_CALL ) {
                    writer.write("0/1");
                } else {
                    writer.write("0/0");
                }
            } else if ( 1.0f * altBetter / totalCounts > MIN_FRACTION_FOR_ALLELE_CALL ) {
                writer.write("1/1");
            } else {
                writer.write("NA");
            }
            writer.newLine();
/*
            for ( final Map.Entry<ScorePair, Integer> entry : scores.entrySet() ) {
                writer.write(Integer.toString(entry.getKey().refScore()));
                writer.write('\t');
                writer.write(Integer.toString(entry.getKey().altScore()));
                writer.write('\t');
                writer.write(entry.getValue().toString());
                writer.newLine();
            }
 */
        }
    }
}
