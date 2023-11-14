package org.broadinstitute.hellbender.tools.walkers.featuremapping;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import org.apache.commons.math3.util.Precision;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.tools.FlowBasedArgumentCollection;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.haplotype.FlowBasedHaplotype;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.FlowBasedRead;
import org.broadinstitute.hellbender.utils.read.FlowBasedReadUtils;

import java.util.Arrays;
import java.util.List;

public class FlowFeatureMapperUtils {

    // some local (mostly debug) configuration
    static class Args {
        SAMFileHeader header;
        FlowBasedArgumentCollection fbArgs;
        boolean debugNegatives;
        List<String> debugReadName;
        double limitScore;
        boolean keepNegatives;
        double negativeScoreOverride = 0;
    };

    private static final Double     LOWEST_PROB = 0.0001;

    private static final Logger logger = LogManager.getLogger(FlowFeatureMapper.class);

    public static double scoreFeature(final MappedFeature fr, Args args) {
        return scoreFeature(fr, (byte)0, args);
    }

    public static double scoreFeature(final MappedFeature fr, byte altBase, Args args) {

        // build haplotypes
        final FlowBasedReadUtils.ReadGroupInfo rgInfo = FlowBasedReadUtils.getReadGroupInfo(args.header, fr.read);
        final FlowBasedHaplotype[]    haplotypes = buildHaplotypes(fr, rgInfo.flowOrder, altBase);

        // create flow read
        final FlowBasedRead flowRead;
        if ( fr.flowRead == null ) {
            flowRead = new FlowBasedRead(fr.read, rgInfo.flowOrder, rgInfo.maxClass, args.fbArgs);
        } else {
            flowRead = new FlowBasedRead(fr.flowRead, true);
        }

        final int diffLeft = haplotypes[0].getStart() - flowRead.getStart() + fr.offsetDelta;
        final int diffRight = flowRead.getEnd() - haplotypes[0].getEnd();
        flowRead.applyBaseClipping(Math.max(0, diffLeft), Math.max(diffRight, 0), false);

        if ( !flowRead.isValid() ) {
            return -1;
        }

        // compute alternative score
        final int         hapKeyLength = Math.min(haplotypes[0].getKeyLength(), haplotypes[1].getKeyLength());
        final double      readScore = computeLikelihoodLocal(flowRead, haplotypes[0], hapKeyLength, false);
        final double      refScore = computeLikelihoodLocal(flowRead, haplotypes[1], hapKeyLength, false);
        double            score = readScore - refScore;
        if ( !Double.isNaN(args.limitScore) ) {
            score = Math.min(score, args.limitScore);
        }

        if ( ((Double.isNaN(score) || (score < 0)) && args.debugNegatives)
                || (args.debugReadName != null && args.debugReadName.contains(fr.read.getName())) ) {
            logger.info("**** debug read: " + fr.read);
            logger.info("readBases: " + fr.read.getBasesString());
            logger.info("flowRead: " + flowRead);
            logger.info("flowBases: " + flowRead.getBasesString());
            logger.info("flowOrder: " + flowRead.getFlowOrder());
            logger.info("flowKey: " + flowRead.getKeyLength() + " " + Arrays.toString(flowRead.getKey()));
            logger.info("readHaplotype: " + haplotypes[0]);
            logger.info("readHapKey: " + haplotypes[0].getKeyLength() + " " + Arrays.toString(haplotypes[0].getKey()));
            computeLikelihoodLocal(flowRead, haplotypes[0], hapKeyLength, true);
            logger.info("refrHaplotype: " + haplotypes[1]);
            logger.info("refrHapKey: " + haplotypes[1].getKeyLength() + " " + Arrays.toString(haplotypes[1].getKey()));
            computeLikelihoodLocal(flowRead, haplotypes[1], hapKeyLength, true);
            logger.info("score: " + score);

            // analyze read
            final FlowBasedRead flowRead2;
            flowRead2 = new FlowBasedRead(fr.read, rgInfo.flowOrder, rgInfo.maxClass, args.fbArgs);
            final int[]        key2 = flowRead2.getKey();
            for ( int i = 0 ; i < key2.length ; i++ ) {
                final double      p1 = flowRead2.getProb(i, key2[i]);
                for ( int j = 0 ; j < rgInfo.maxClass ; j++ ) {
                    final double      p2 = flowRead2.getProb(i, j);
                    if ( p2 > p1 )
                        logger.info(String.format("prob at %s key[%d]=%d, %f is lower than at %d which is %f",
                                flowRead2.getName(), i, key2[i], p1, j, p2));
                }
            }
        }

        if ( score < 0 && !args.keepNegatives && score != -1.0 ) {
            score = args.negativeScoreOverride;
        }

        return score;
    }

    public static double computeLikelihoodLocal(final FlowBasedRead read, final FlowBasedHaplotype haplotype, final int hapKeyLength, final boolean debug) {

        final byte[] flowOrder = haplotype.getFlowOrderArray();
        final byte   readFlowOrder0 = read.getFlowOrderArray()[0];
        int startingPoint = 0;
        for (int i = 0; i < flowOrder.length; i++) {
            if (flowOrder[i] == readFlowOrder0) {
                startingPoint = i;
                break;
            }
        }
        final int[]         key = haplotype.getKey();

        // debug support
        StringBuffer        debugMessage = null;
        if ( debug )
            debugMessage = new StringBuffer(Integer.toString(startingPoint) + " hmer prob |");
        double              result = 0 ;
        for (int i = 0; i < read.getKeyLength(); i++) {
            int     index = i + startingPoint;
            double  prob = 0;
            int     locationToFetch = 0;
            if ( index < hapKeyLength ) {
                locationToFetch = Math.min(key[index] & 0xff, read.getMaxHmer() + 1);
                prob = read.getProb(i, locationToFetch);
            } else {
                if ( debug ) {
                    debugMessage.append(" clip");
                }
                break;
            }
            if ( Precision.equals(prob, 0.0) ) {
                prob = LOWEST_PROB;
            }
            result += Math.log10(prob);

            if ( debug ) {
                debugMessage.append(String.format(" %d %.4f", locationToFetch, prob));
            }
        }

        if ( debug ) {
            debugMessage.append(" | " + result);
            logger.info("debugMessage: " + debugMessage);
        }

        return result;
    }

    static FlowBasedHaplotype[] buildHaplotypes(final MappedFeature fr, final String flowOrder, byte altBase) {

        // build bases for flow haplotypes
        // NOTE!!!: this code assumes length of feature on read and reference is the same
        // this is true for SNP but not for INDELs - it will have to be re-written!
        // TODO: write for INDEL
        byte[] bases = fr.read.getBasesNoCopy();
        int         offset = fr.readBasesOffset;
        int         refStart = fr.start;
        int         refModOfs = 0;

        // install alt base?
        byte orgBase = 0;
        if ( altBase != 0 ) {
            orgBase = fr.refBases[0];
            fr.refBases[0] = altBase;
        }

        if ( offset > 0 ) {
            // reach into hmer before
            offset--;
            refModOfs++;
            refStart--;

            // extend until start of hmer
            final byte        hmerBase = bases[offset];
            while ( offset > 0 && bases[offset-1] == hmerBase ) {
                offset--;
                refModOfs++;
                refStart--;
            }
        }
        final byte[]      sAltBases = Arrays.copyOfRange(bases, offset, bases.length);
        final byte[]      sRefBases = Arrays.copyOf(sAltBases, sAltBases.length);
        System.arraycopy(fr.refBases, 0, sRefBases, refModOfs, fr.refBases.length);

        // construct haplotypes
        final SimpleInterval genomeLoc = new SimpleInterval(fr.read.getContig(), refStart, refStart + sAltBases.length - 1);
        final Cigar cigar = new Cigar();
        cigar.add(new CigarElement(sAltBases.length, CigarOperator.M));
        final Haplotype altHaplotype = new Haplotype(sAltBases, false);
        final Haplotype      refHaplotype = new Haplotype(sRefBases, true);
        altHaplotype.setGenomeLocation(genomeLoc);
        refHaplotype.setGenomeLocation(genomeLoc);
        altHaplotype.setCigar(cigar);
        refHaplotype.setCigar(cigar);

        // prepare flow based haplotypes
        final FlowBasedHaplotype[] result = {
                new FlowBasedHaplotype(altHaplotype, flowOrder),
                new FlowBasedHaplotype(refHaplotype, flowOrder)
        };

        // restore changes
        if ( altBase != 0 ) {
            fr.refBases[0] = orgBase;
        }

        // return
        return result;
    }

}
