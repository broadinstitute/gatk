package org.broadinstitute.hellbender.tools.walkers;

import edu.mit.broad.tedsUtils.align.EndsFreeAligner;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.util.SequenceUtil;
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
import org.broadinstitute.hellbender.utils.read.GATKRead;
import picard.cmdline.programgroups.VariantEvaluationProgramGroup;

import java.io.*;
import java.util.Iterator;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

@CommandLineProgramProperties(
        summary = "For each SV, write a couple of snippets of sequence to look for.",
        oneLineSummary = "For each SV, write a couple of snippets of sequence to look for.",
        programGroup = VariantEvaluationProgramGroup.class
)
public final class GenotypeSVs extends VariantWalker {
    public static final int MINIMUM_SV_SIZE = 50;
    public static final int WINDOW_SIZE = 50;

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
    public void apply( final VariantContext variant,
                       final ReadsContext readsContext,
                       final ReferenceContext referenceContext,
                       final FeatureContext featureContext) {
        final String svType;
        final int svLen;
        if ( variant.hasAttribute(GATKSVVCFConstants.SVTYPE) ) {
            svType = variant.getAttributeAsString(GATKSVVCFConstants.SVTYPE, null);
            svLen = Math.abs(variant.getAttributeAsInt(GATKSVVCFConstants.SVLEN, 0));
        } else if ( variant.isIndel() ) {
            final int refLen = variant.getReference().getBaseString().length();
            final int altLen = variant.getAlternateAllele(0).getBaseString().length();
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

        final String contig = variant.getContig();
        final int start = variant.getStart();
        final int end = variant.getEnd();

        switch ( svType ) {
            case GATKSVVCFConstants.SYMB_ALT_STRING_DEL ->
                    testDeletion(contig, start, end, referenceContext, readsContext);
/*
            case GATKSVVCFConstants.SYMB_ALT_STRING_INS -> {
                final String altBases = variant.getAlternateAllele(0).getBaseString();
                writeRef(seqName + "1", leadBases + altBases.substring(1, WINDOW_SIZE + 1));
                writeRef(seqName + "2", altBases.substring(altBases.length() - WINDOW_SIZE) + lagBases);
            }
            case GATKSVVCFConstants.SYMB_ALT_STRING_DUP -> {
                final SimpleInterval dupInterval = new SimpleInterval(contig, start + 1, end);
                final String altBases = new String(referenceContext.getBases(dupInterval));
                writeRef(seqName + "1", leadBases + altBases.substring(1, WINDOW_SIZE + 1));
                writeRef(seqName + "2", altBases.substring(altBases.length() - WINDOW_SIZE) + lagBases);
            }
            case GATKSVVCFConstants.SYMB_ALT_STRING_INV -> {
                final SimpleInterval invInterval = new SimpleInterval(contig, start + 1, end);
                final byte[] invBaseCalls = referenceContext.getBases(invInterval);
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
                            new String(referenceContext.getBases(bndInterval)) + matcher.group(3);
                    writeRef(seqName, altBases.substring(altBases.length() - WINDOW_SIZE) + lagBases);
                } else if ( (matcher = BND_PRE_REV.matcher(bndDescription)).matches() ) {
                    final String tig = matcher.group(1);
                    final int bndPosBeg = Integer.parseInt(matcher.group(2));
                    final int bndPosEnd = bndPosBeg + WINDOW_SIZE - 1;
                    final SimpleInterval bndInterval = new SimpleInterval(tig, bndPosBeg, bndPosEnd);
                    final byte[] bndBaseCalls = referenceContext.getBases(bndInterval);
                    SequenceUtil.reverseComplement(bndBaseCalls);
                    final String altBases = new String(bndBaseCalls) + matcher.group(3);
                    writeRef(seqName, altBases.substring(altBases.length() - WINDOW_SIZE) + lagBases);
                } else if ( (matcher = BND_POST_FWD.matcher(bndDescription)).matches() ) {
                    final String tig = matcher.group(2);
                    final int bndPosBeg = Integer.parseInt(matcher.group(3));
                    final int bndPosEnd = bndPosBeg + WINDOW_SIZE - 1;
                    final SimpleInterval bndInterval = new SimpleInterval(tig, bndPosBeg, bndPosEnd);
                    final String altBases =
                            matcher.group(1) + new String(referenceContext.getBases(bndInterval));
                    writeRef(seqName, leadBases + altBases.substring(0, WINDOW_SIZE));
                } else if ( (matcher = BND_POST_REV.matcher(bndDescription)).matches() ) {
                    final String tig = matcher.group(2);
                    final int bndPosEnd = Integer.parseInt(matcher.group(3));
                    final int bndPosBeg = Math.max(1, bndPosEnd - WINDOW_SIZE + 1);
                    final SimpleInterval bndInterval = new SimpleInterval(tig, bndPosBeg, bndPosEnd);
                    final byte[] bndBaseCalls = referenceContext.getBases(bndInterval);
                    SequenceUtil.reverseComplement(bndBaseCalls);
                    final String altBases = matcher.group(1) + new String(bndBaseCalls);
                    writeRef(seqName, leadBases + altBases.substring(0, WINDOW_SIZE));
                } else {
                    logger.warn("Can't interpret the BND description: " + bndDescription);
                }
            }
*/
            default -> logger.warn("Don't know how to handle SVTYPE=" + svType + " at " + contig + ":" + start);
        }
    }

    private void testDeletion( final String contig, final int start, final int end,
                               final ReferenceContext referenceContext,
                               final ReadsContext readsContext ) {
        final int alignStart = Math.max(1, start - 5);
        final int alignEnd = end + WINDOW_SIZE;
        final String refCalls =
            new String(referenceContext.getBases(new SimpleInterval(contig, alignStart, start))) +
              new String(referenceContext.getBases(new SimpleInterval(contig, end + 1, alignEnd)));
        final Iterator<GATKRead> readItr =
                readsContext.iterator(new SimpleInterval(contig, alignStart, alignEnd));
        while ( readItr.hasNext() ) {
            final GATKRead read = readItr.next();
            if ( read.getStart() > alignStart ) {
                break;
            }
            final int readStartPos = getReadPositionForRefPosition(read, alignStart);
            if ( readStartPos == -1 ) {
                continue;
            }
            final String allReadCalls = read.getBasesString();
            final int readEndPos = Math.min(readStartPos + refCalls.length(), allReadCalls.length());
            final String readCalls = allReadCalls.substring(readStartPos, readEndPos);
            final EndsFreeAligner aligner = new EndsFreeAligner(readCalls, refCalls);
            System.out.println(aligner.getAlignment());
        }
    }

    private int getReadPositionForRefPosition( final GATKRead read, final int refPosition ) {
        int refPos = read.getStart();
        int readPos = 0;
        for ( final CigarElement cigarElement : read.getCigarElements() ) {
            final int opLen = cigarElement.getLength();
            final CigarOperator op = cigarElement.getOperator();
            if ( op.consumesReferenceBases() ) {
                final int nextRefPos = refPos + opLen;
                if ( nextRefPos > refPosition ) {
                    if ( op.consumesReadBases() ) {
                        readPos += opLen + nextRefPos - refPosition;
                    }
                    return readPos;
                }
                refPos += opLen;
            }
            if ( op.consumesReadBases() ) {
                readPos += opLen;
            }
        }
        return -1;
    }
}
