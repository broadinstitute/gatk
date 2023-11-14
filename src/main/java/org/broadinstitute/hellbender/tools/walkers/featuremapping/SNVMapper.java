package org.broadinstitute.hellbender.tools.walkers.featuremapping;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.text.similarity.LevenshteinDistance;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.FlowBasedRead;

import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.function.Consumer;

/**
 * An implementation of a feature mapper that finds SNPs (SVN)
 *
 * This class only finds SNP that are surrounded by a specific number of bases identical to the reference.
 */

public class SNVMapper implements FeatureMapper {

    static final int SUROUND_DEFUALT = 5;

    final int         surroundBefore;
    final int         surroundAfter;
    final LevenshteinDistance levDistance = new LevenshteinDistance();
    final Integer     smqSize;
    final Integer     smqSizeMean;
    final boolean     reportAllAlts;
    final boolean     tagBasesWithAdjacentRefDiff;
    final boolean     generateNonRefFeatures;

    boolean     ignoreSurround;
    int         spanBefore;
    int         spanAfter;
    int         minCigarElementLength;

    public SNVMapper(FlowFeatureMapperArgumentCollection fmArgs) {
        surroundBefore = fmArgs.snvIdenticalBases;
        surroundAfter = (fmArgs.snvIdenticalBasesAfter != 0) ?  fmArgs.snvIdenticalBasesAfter : surroundBefore;
        smqSize = fmArgs.surroundingMediaQualitySize;
        smqSizeMean = fmArgs.surroundingMeanQualitySize;
        reportAllAlts = fmArgs.reportAllAlts;
        tagBasesWithAdjacentRefDiff = fmArgs.tagBasesWithAdjacentRefDiff;
        generateNonRefFeatures = false;

        init();
    }

    public SNVMapper(ApplySNVQRArgumentCollection aqArgs) {
        surroundBefore = SUROUND_DEFUALT;
        surroundAfter = SUROUND_DEFUALT;
        smqSize = null;
        smqSizeMean = null;
        reportAllAlts = true;
        tagBasesWithAdjacentRefDiff = true;
        generateNonRefFeatures = true;

        init();
    }

    private void init() {

        // ignore surround
        ignoreSurround = reportAllAlts || tagBasesWithAdjacentRefDiff;
        spanBefore = ignoreSurround ? 0 : surroundBefore;
        spanAfter = ignoreSurround ? 0 : surroundAfter;
        minCigarElementLength = spanBefore + 1 + spanAfter;

        // adjust minimal read length
        FlowBasedRead.setMinimalReadLength(1 + 1 + spanAfter);
    }

    @Override
    public void forEachOnRead(GATKRead read, ReferenceContext referenceContext, Consumer<? super MappedFeature> action) {

        // prepare list
        List<MappedFeature>     features = new LinkedList<>();

        // access bases
        final byte[]      bases = read.getBasesNoCopy();
        final byte[]      ref = referenceContext.getBases();

        // calculate edit distance
        int               startSoftClip = read.getStart() - read.getSoftStart();
        int               endSoftClip = read.getSoftEnd() - read.getEnd();
        String            basesString;
        if ( startSoftClip == 0 && endSoftClip == 0 ) {
            basesString = new String(bases);
        } else {
            basesString = new String(Arrays.copyOfRange(bases, startSoftClip, bases.length - endSoftClip));
        }
        int               refEditDistance = calcEditDistance(basesString, new String(ref));

        // count bases delta on M cigar elements
        int         nonIdentMBases = 0;
        int         readOfs = 0;
        int         refOfs = 0;
        int         cigarElementCount = read.getCigar().numCigarElements();
        for ( int cigarElementIndex = 0 ; cigarElementIndex < cigarElementCount ; cigarElementIndex++ ) {
            final CigarElement cigarElement = read.getCigar().getCigarElement(cigarElementIndex);
            final int     length = cigarElement.getLength();
            if ( cigarElement.getOperator().consumesReadBases() && cigarElement.getOperator().consumesReferenceBases() ) {
                for ( int ofs = 0 ; ofs < length ; ofs++ ) {
                    if ( ref[refOfs+ofs] != 'N' && bases[readOfs+ofs] != ref[refOfs+ofs] ) {
                        nonIdentMBases++;
                    }
                }
            }
            if (cigarElement.getOperator().consumesReadBases()) {
                readOfs += length;
            }
            if (cigarElement.getOperator().consumesReferenceBases()) {
                refOfs += length;
            }
        }
        int         hardLength = read.getUnclippedEnd() - read.getUnclippedStart() + 1;

        // walk the cigar (again) looking for features
        readOfs = 0;
        refOfs = 0;
        int numCigarElements = read.numCigarElements();
        for ( int cigarElementIndex = 0 ; cigarElementIndex < numCigarElements ; cigarElementIndex++ ) {
            final CigarElement cigarElement = read.getCigarElement(cigarElementIndex);
            final int     length = cigarElement.getLength();

            // worth looking into?
            if ( length >= minCigarElementLength &&
                    cigarElement.getOperator().consumesReadBases() &&
                    (generateNonRefFeatures || cigarElement.getOperator().consumesReferenceBases()) ) {
                readOfs += spanBefore;
                refOfs += spanBefore;
                final boolean hasRef = cigarElement.getOperator().consumesReferenceBases();
                for ( int ofs = spanBefore ; ofs < length - spanAfter ; ofs++, readOfs++, refOfs += (hasRef) ? 1 : 0 ) {

                    final byte refBase = hasRef ? ref[refOfs] : bases[readOfs];
                    if ( refBase != 'N' && (reportAllAlts || (bases[readOfs] != refBase)) ) {

                        // check that this is really a SNV (must be surrounded by identical ref)
                        boolean     surrounded = true;
                        for ( int i = 0 ; i < surroundBefore && surrounded ; i++ ) {
                            final int bIndex = readOfs-1-i;
                            final int rIndex = refOfs-1-i;
                            if ( bIndex < 0 || bIndex >= bases.length || rIndex < 0 || rIndex >= ref.length ) {
                                surrounded = false;
                                continue;
                            }
                            if ( !hasRef || (bases[bIndex] != ref[rIndex]) ) {
                                surrounded = false;
                            }
                        }
                        for (int i = 0; i < surroundAfter && surrounded ; i++ ) {
                            final int bIndex = readOfs+1+i;
                            final int rIndex = refOfs+1+i;
                            if ( bIndex < 0 || bIndex >= bases.length || rIndex < 0 || rIndex >= ref.length ) {
                                surrounded = false;
                                continue;
                            }
                            if ( !hasRef || (bases[bIndex] != ref[rIndex]) ) {
                                surrounded = false;
                            }
                        }
                        if ( (!reportAllAlts && !tagBasesWithAdjacentRefDiff) && !surrounded ) {
                            continue;
                        }

                        // add this feature
                        MappedFeature feature = MappedFeature.makeSNV(read, readOfs, refBase, referenceContext.getStart() + refOfs, readOfs - refOfs);
                        if ( (reportAllAlts || tagBasesWithAdjacentRefDiff) && !surrounded )
                            feature.adjacentRefDiff = true;
                        feature.nonIdentMBasesOnRead = nonIdentMBases;
                        feature.refEditDistance = refEditDistance;
                        if ( !read.isReverseStrand() )
                            feature.index = readOfs;
                        else
                            feature.index = hardLength - readOfs;

                        // reverse quals?
                        if ( smqSize != null || smqSizeMean != null  ) {

                            // prepare qualities
                            final byte[] quals;
                            if (!read.isReverseStrand()) {
                                quals = read.getBaseQualitiesNoCopy();
                            } else {
                                quals = read.getBaseQualities();
                                ArrayUtils.reverse(quals);
                            }

                            // surrounding median quality?
                            if ( smqSize != null ) {
                                feature.smqLeft = calcSmq(quals, feature.index - 1 - smqSize, feature.index - 1, true);
                                feature.smqRight = calcSmq(quals, feature.index + 1, feature.index + 1 + smqSize, true);
                                if ( read.isReverseStrand() ) {

                                    // left and right are reversed
                                    int tmp = feature.smqLeft;
                                    feature.smqLeft = feature.smqRight;
                                    feature.smqRight = tmp;
                                }
                            }
                            if ( smqSizeMean != null ) {
                                feature.smqLeftMean = calcSmq(quals, feature.index - 1 - smqSizeMean, feature.index - 1, false);
                                feature.smqRightMean = calcSmq(quals, feature.index + 1, feature.index + 1 + smqSizeMean, false);
                                if ( read.isReverseStrand() ) {

                                    // left and right are reversed
                                    int tmp = feature.smqLeftMean;
                                    feature.smqLeftMean = feature.smqRightMean;
                                    feature.smqRightMean = tmp;
                                }
                            }

                        }

                        features.add(feature);
                    }
                }
                readOfs += spanAfter;
                refOfs += spanAfter;

            } else {

                // manual advance
                if (cigarElement.getOperator().consumesReadBases()) {
                    readOfs += length;
                }
                if (cigarElement.getOperator().consumesReferenceBases()) {
                    refOfs += length;
                }
            }
        };

        // report features
        for ( MappedFeature feature : features ) {
            feature.featuresOnRead = features.size();
            action.accept(feature);
        }
    }

    private int calcEditDistance(String s1, String s2) {

        // identical
        if ( s1.length() == s2.length() && s1.equals(s2) ) {
            return 0;
        }

        // find first difference & trim
        int     minLength = Math.min(s1.length(), s2.length());
        int     diffIndex = 0;
        for ( ; diffIndex < minLength ; diffIndex++ ) {
            if (s1.charAt(diffIndex) != s2.charAt(diffIndex)) {
                break;
            }
        }
        s1 = s1.substring(diffIndex);
        s2 = s2.substring(diffIndex);
        minLength = Math.min(s1.length(), s2.length());

        // find last difference & trim
        int diffLength1 = s1.length();
        int diffLength2 = s2.length();
        for ( ; diffLength1 > 0 && diffLength2 > 0 ; diffLength1--, diffLength2-- ) {
            if (s1.charAt(diffLength1 - 1) != s2.charAt(diffLength2 - 1)) {
                break;
            }
        }
        s1 = s1.substring(0, diffLength1);
        s2 = s2.substring(0, diffLength2);

        // find edit distance on the shorter strings
        return levDistance.apply(s1, s2);
    }

    private int calcSmq(final byte[] quals, int from, int to, boolean median) {

        // limit from/to
        from = Math.max(0, Math.min(quals.length, from));
        to = Math.max(0, Math.min(quals.length, to - 1));
        if ( from > to ) {
            throw new GATKException("invalid qualities range: from > to");
        }

        // calc median
        byte[] range = Arrays.copyOfRange(quals, from, to + 1);
        if ( range.length == 0 ) {
            throw new GATKException("invalid qualities range: can't be empty");
        }

        if ( median ) {
            Arrays.sort(range);
            int midIndex = range.length / 2;
            if ((range.length % 2) == 1) {
                // odd
                return range[midIndex];
            } else {
                // even
                return (range[midIndex - 1] + range[midIndex]) / 2;
            }
        } else {
            int sum = 0;
            for ( int i = 0 ; i < range.length ; i++ ) {
                sum += range[i];
            }
            return sum / range.length;
        }
    }

    public FilterStatus noFeatureButFilterAt(GATKRead read, ReferenceContext referenceContext, int start) {

        // access bases
        final byte[]      bases = read.getBasesNoCopy();
        final byte[]      ref = referenceContext.getBases();

        // walk the cigar
        int         readOfs = 0;
        int         refOfs = 0;
        int numCigarElements = read.numCigarElements();
        for ( int cigarElementIndex = 0 ; cigarElementIndex < numCigarElements ; cigarElementIndex++ ) {
            final CigarElement cigarElement = read.getCigarElement(cigarElementIndex);

            final int     length = cigarElement.getLength();

            // worth looking into?
            boolean     includes = (start >= referenceContext.getStart() + refOfs) &&
                    (start < referenceContext.getStart() + refOfs + length);
            if ( includes && length >= minCigarElementLength &&
                    cigarElement.getOperator().consumesReadBases() &&
                    cigarElement.getOperator().consumesReferenceBases() ) {

                // break out if not enough clearing
                if ( (start < referenceContext.getStart() + refOfs + spanBefore) ||
                        (start >= referenceContext.getStart() + refOfs + length - spanAfter) )
                    return FilterStatus.Filtered;

                int         delta = start - (referenceContext.getStart() + refOfs);
                readOfs += delta;
                refOfs += delta;

                if ( bases[readOfs] == ref[refOfs] ) {

                    // check that this is really a SNV (must be surrounded by identical ref)
                    boolean     surrounded = true;
                    for ( int i = 0 ; i < surroundBefore && surrounded ; i++ ) {
                        final int bIndex = readOfs-1-i;
                        final int rIndex = refOfs-1-i;
                        if ( bIndex < 0 || bIndex >= bases.length || rIndex < 0 || rIndex >= ref.length ) {
                            surrounded = false;
                            continue;
                        }
                        if ( bases[bIndex] != ref[rIndex] ) {
                            surrounded = false;
                        }
                    }
                    for (int i = 0; i < surroundAfter && surrounded ; i++ ) {
                        final int bIndex = readOfs+1+i;
                        final int rIndex = refOfs+1+i;
                        if ( bIndex < 0 || bIndex >= bases.length || rIndex < 0 || rIndex >= ref.length ) {
                            surrounded = false;
                            continue;
                        }
                        if ( bases[bIndex] != ref[rIndex] ) {
                            surrounded = false;
                        }
                    }
                    if ( !surrounded ) {
                        continue;
                    }

                    // this is it! no feature but filtered in
                    return FilterStatus.NoFeatureAndFiltered;
                } else
                    return FilterStatus.Filtered;

            } else {

                // manual advance
                if (cigarElement.getOperator().consumesReadBases()) {
                    readOfs += length;
                }
                if (cigarElement.getOperator().consumesReferenceBases()) {
                    refOfs += length;
                }
            }
        };

        // if here, false
        return FilterStatus.None;
    }

}
