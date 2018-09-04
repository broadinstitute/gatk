package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalTree;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.*;
import java.util.function.Function;

/**
 * Figures out what kind of BreakpointEvidence, if any, a read represents.
 */
public class ReadClassifier implements Function<GATKRead, Iterator<BreakpointEvidence>> {
    @VisibleForTesting static final int MIN_SOFT_CLIP_LEN = 30; // minimum length of an interesting soft clip
    @VisibleForTesting static final int MIN_INDEL_LEN = 40; // minimum length of an interesting indel
    private static final byte MIN_QUALITY = 15; // minimum acceptable quality in a soft-clip window
    private static final int MAX_LOW_QUALITY_SCORES = 3; // maximum # of low quality base calls in soft-clip window
    private static final float MAX_NON_OUTLIER_ZISH_SCORE = 6.f; // maximum fragment-length "z" score for a normal fragment
    private final ReadMetadata readMetadata;
    private final GATKRead sentinel;
    private final int allowedShortFragmentOverhang;
    private final SVReadFilter filter;
    private final KSWindowFinder smallIndelFinder;
    private final SVIntervalTree<SVInterval> regionsToIgnore;

    public ReadClassifier(final ReadMetadata readMetadata,
                          GATKRead sentinel,
                          final int allowedShortFragmentOverhang,
                          SVReadFilter filter,
                          final SVIntervalTree<SVInterval> regionsToIgnore) {
        this.readMetadata = readMetadata;
        this.sentinel = sentinel;
        this.allowedShortFragmentOverhang = allowedShortFragmentOverhang;
        this.filter = filter;
        this.regionsToIgnore = regionsToIgnore;
        smallIndelFinder = new KSWindowFinder(readMetadata, filter);
    }

    @Override
    public Iterator<BreakpointEvidence> apply( final GATKRead read ) {
        if ( read == sentinel ) {
            final List<BreakpointEvidence> evidenceList = new ArrayList<>();
            smallIndelFinder.checkHistograms(evidenceList);
            return evidenceList.iterator();
        }

        if ( !filter.isMappedToPrimaryContig(read, readMetadata)) return Collections.emptyIterator();
        if ( !filter.isEvidence(read) ) return Collections.emptyIterator();
        if (regionsToIgnore != null) {
            final int readContigId = readMetadata.getContigID(read.getContig());
            final SVInterval clippedReadInterval = new SVInterval(readContigId, read.getStart(), read.getEnd());
            if (filter.containedInRegionToIgnore(clippedReadInterval, regionsToIgnore)) {
                return Collections.emptyIterator();
            }
        }

        final List<BreakpointEvidence> evidenceList = new ArrayList<>();
        checkForSplitRead(read, evidenceList);

        checkDiscordantPair(read, evidenceList);
        smallIndelFinder.testReadAndGatherEvidence(read, evidenceList);

        return evidenceList.iterator();
    }

    private void checkForSplitRead( final GATKRead read, final List<BreakpointEvidence> evidenceList ) {
        final List<CigarElement> cigarElements = read.getCigar().getCigarElements();
        if ( hasInitialSoftClip(cigarElements, read) ) {
            evidenceList.add(new BreakpointEvidence.SplitRead(read, readMetadata, true));
        }
        if ( hasFinalSoftClip(cigarElements, read) ) {
            evidenceList.add(new BreakpointEvidence.SplitRead(read, readMetadata, false));
        }
        checkBigIndel(cigarElements, read, evidenceList);
    }

    private static boolean hasInitialSoftClip( final List<CigarElement> cigarElements, final GATKRead read ) {
        final ListIterator<CigarElement> itr = cigarElements.listIterator();
        if ( !itr.hasNext() ) return false;

        CigarElement firstEle = itr.next();
        if ( firstEle.getOperator() == CigarOperator.HARD_CLIP && itr.hasNext() ) {
            firstEle = itr.next();
        }
        final int clipStart = firstEle.getLength() - MIN_SOFT_CLIP_LEN;
        return firstEle.getOperator() == CigarOperator.SOFT_CLIP &&
                clipStart >= 0 &&
                isHighQualityRegion(read.getBaseQualities(), clipStart);
    }

    private static boolean hasFinalSoftClip( final List<CigarElement> cigarElements, final GATKRead read ) {
        final ListIterator<CigarElement> itr = cigarElements.listIterator(cigarElements.size());
        if ( !itr.hasPrevious() ) return false;

        CigarElement lastEle = itr.previous();
        if ( lastEle.getOperator() == CigarOperator.HARD_CLIP && itr.hasPrevious() ) {
            lastEle = itr.previous();
        }
        return lastEle.getOperator() == CigarOperator.SOFT_CLIP &&
                lastEle.getLength() >= MIN_SOFT_CLIP_LEN &&
                isHighQualityRegion(read.getBaseQualities(), read.getLength() - lastEle.getLength());
    }

    private static boolean isHighQualityRegion( final byte[] quals, int idx ) {
        int lowQuals = 0;
        for ( final int end = idx+MIN_SOFT_CLIP_LEN; idx != end; ++idx ) {
            if ( quals[idx] < MIN_QUALITY ) {
                lowQuals += 1;
                if ( lowQuals > MAX_LOW_QUALITY_SCORES ) return false;
            }
   
        }
        return true;
    }

    private void checkBigIndel( final List<CigarElement> cigarElements,
                               final GATKRead read,
                               final List<BreakpointEvidence> evidenceList ) {
        int locus = read.getStart();
        for ( final CigarElement ele : cigarElements ) {
            final CigarOperator op = ele.getOperator();
            if ( ele.getLength() >= MIN_INDEL_LEN ) {
                if ( op == CigarOperator.INSERTION ) {
                    evidenceList.add(new BreakpointEvidence.LargeIndel(read, readMetadata, locus));
                } else if ( op == CigarOperator.DELETION ) {
                    evidenceList.add(new BreakpointEvidence.LargeIndel(read, readMetadata, locus+ele.getLength()/2));
                }
            }
            if ( op.consumesReferenceBases() ) {
                locus += ele.getLength();
            }
        }
    }

    private void checkDiscordantPair( final GATKRead read, final List<BreakpointEvidence> evidenceList ) {
        if (! filter.isPrimaryLine(read)) return;

        if ( read.mateIsUnmapped() ) {
            evidenceList.add(new BreakpointEvidence.MateUnmapped(read, readMetadata));
        } else {

            final int contigID1 = readMetadata.getContigID(read.getContig());
            final int contigID2 = readMetadata.getContigID(read.getMateContig());

            if (contigID1 != contigID2) {
                if ( !readMetadata.ignoreCrossContigID(contigID1) && !readMetadata.ignoreCrossContigID(contigID2) ) {
                    evidenceList.add(new BreakpointEvidence.InterContigPair(read, readMetadata));
                }
            } else {
                if ( read.isReverseStrand() == read.mateIsReverseStrand() ) {
                    evidenceList.add(new BreakpointEvidence.SameStrandPair(read, readMetadata));
                } else if (read.isReverseStrand() ?
                        read.getStart() + allowedShortFragmentOverhang < read.getMateStart() :
                        read.getStart() - allowedShortFragmentOverhang > read.getMateStart()) {
                    evidenceList.add(new BreakpointEvidence.OutiesPair(read, readMetadata));
                } else {
                    final float zIshScore = readMetadata.getZishScore(read.getReadGroup(), Math.abs(read.getFragmentLength()));
                    if ( zIshScore > MAX_NON_OUTLIER_ZISH_SCORE) {
                        evidenceList.add(new BreakpointEvidence.WeirdTemplateSize(read, readMetadata));
                    }
                }
            }
        }
    }
}
