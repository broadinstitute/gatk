package org.broadinstitute.hellbender.tools.walkers.featuremapping;

import htsjdk.samtools.CigarElement;
import org.apache.commons.text.similarity.LevenshteinDistance;
import org.broadinstitute.hellbender.engine.ReferenceContext;
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

    final int         identBefore;
    final int         identAfter;
    final int         minCigarElementLength;
    final LevenshteinDistance levDistance = new LevenshteinDistance();

    public SNVMapper(FlowFeatureMapperArgumentCollection fmArgs) {
        identBefore = fmArgs.snvIdenticalBases;
        identAfter = (fmArgs.snvIdenticalBasesAfter != 0) ?  fmArgs.snvIdenticalBasesAfter : identBefore;
        minCigarElementLength = identBefore + 1 + identAfter;

        // adjust minimal read length
        FlowBasedRead.setMinimalReadLength(1 + 1 + identAfter);
    }

    @Override
    public void forEachOnRead(GATKRead read, ReferenceContext referenceContext, Consumer<? super FlowFeatureMapper.MappedFeature> action) {

        // prepare list
        List<FlowFeatureMapper.MappedFeature>     features = new LinkedList<>();

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
        int               refEditDistance = levDistance.apply(basesString, new String(ref));

        // count bases delta on M cigar elements
        int         nonIdentMBases = 0;
        int         readOfs = 0;
        int         refOfs = 0;
        for ( final CigarElement cigarElement : read.getCigarElements() ) {
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
        for ( final CigarElement cigarElement : read.getCigarElements() ) {

            final int     length = cigarElement.getLength();

            // worth looking into?
            if ( length >= minCigarElementLength &&
                    cigarElement.getOperator().consumesReadBases() &&
                    cigarElement.getOperator().consumesReferenceBases() ) {
                readOfs += identBefore;
                refOfs += identBefore;
                for ( int ofs = identBefore ; ofs < length - identAfter ; ofs++, readOfs++, refOfs++ ) {

                    if ( ref[refOfs] != 'N'  && bases[readOfs] != ref[refOfs] ) {

                        // check that this is really a SNV (must be surrounded by identical ref)
                        boolean     surrounded = true;
                        for ( int i = 0 ; i < identBefore && surrounded ; i++ ) {
                            if ( bases[readOfs-1-i] != ref[refOfs-1-i] ) {
                                surrounded = false;
                            }
                        }
                        for ( int i = 0 ; i < identAfter && surrounded ; i++ ) {
                            if ( bases[readOfs+1+i] != ref[refOfs+1+i] ) {
                                surrounded = false;
                            }
                        }
                        if ( !surrounded ) {
                            continue;
                        }

                        // add this feature
                        FlowFeatureMapper.MappedFeature feature = FlowFeatureMapper.MappedFeature.makeSNV(read, readOfs, ref[refOfs], referenceContext.getStart() + refOfs, readOfs - refOfs);
                        feature.nonIdentMBasesOnRead = nonIdentMBases;
                        feature.refEditDistance = refEditDistance;
                        if ( !read.isReverseStrand() )
                            feature.index = readOfs;
                        else
                            feature.index = hardLength - readOfs;
                        features.add(feature);
                    }
                }
                readOfs += identAfter;
                refOfs += identAfter;

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
        for ( FlowFeatureMapper.MappedFeature feature : features ) {
            feature.featuresOnRead = features.size();
            action.accept(feature);
        }
    }

    public boolean noFeatureButFilterAt(GATKRead read, ReferenceContext referenceContext, int start) {

        // access bases
        final byte[]      bases = read.getBasesNoCopy();
        final byte[]      ref = referenceContext.getBases();

        // walk the cigar
        int         readOfs = 0;
        int         refOfs = 0;
        for ( final CigarElement cigarElement : read.getCigarElements() ) {

            final int     length = cigarElement.getLength();

            // worth looking into?
            boolean     includes = (start >= referenceContext.getStart() + refOfs) &&
                    (start < referenceContext.getStart() + refOfs + length);
            if ( includes && length >= minCigarElementLength &&
                    cigarElement.getOperator().consumesReadBases() &&
                    cigarElement.getOperator().consumesReferenceBases() ) {

                // break out if not enough clearing
                if ( (start < referenceContext.getStart() + refOfs + identBefore) ||
                        (start >= referenceContext.getStart() + refOfs + length - identAfter) )
                    return false;

                int         delta = start - (referenceContext.getStart() + refOfs);
                readOfs += delta;
                refOfs += delta;

                if ( bases[readOfs] == ref[refOfs] ) {

                    // check that this is really a SNV (must be surrounded by identical ref)
                    boolean     surrounded = true;
                    for ( int i = 0 ; i < identBefore && surrounded ; i++ ) {
                        if ( bases[readOfs-1-i] != ref[refOfs-1-i] ) {
                            surrounded = false;
                        }
                    }
                    for ( int i = 0 ; i < identAfter && surrounded ; i++ ) {
                        if ( bases[readOfs+1+i] != ref[refOfs+1+i] ) {
                            surrounded = false;
                        }
                    }
                    if ( !surrounded ) {
                        continue;
                    }

                    // this is it! no feature but filtred in
                    return true;
                } else
                    return false;

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
        return false;
    }

}
