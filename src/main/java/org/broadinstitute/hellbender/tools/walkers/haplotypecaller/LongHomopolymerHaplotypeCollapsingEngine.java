package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.SequenceUtil;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWOverhangStrategy;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWParameters;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAligner;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAlignment;

import java.util.*;

/**
 * Utility class, useful for flow based applications, implementing a workaround for long homopolymers handling.
 *
 * <h3>Why is LongHomopolymerHaplotypeCollapsingEngine useful</h3>
 *
 * Flow based applications represent sequences of bases as a series of homopolymers encoded as a flow key.
 * The flow key  consists of the lengths for each homopolymers, with spacers (0 length) according to a running flow order.
 * (a full definition of flow space encoding is beyond the scope of this class).
 *
 * Lengths expressed in the key are limited to small number, usually between 10 and 12. This corrosponds to
 * limitations present in the source of the flow reads.
 *
 * References and haplotypes derived from them do not possess this limitation. Therefore, they may contain homopolymers of lengths
 * exceeding the ones represented by flow based reads. This presents a challenge when performing operations combining
 * flow based reads and haplotypes originating from a reference. One such application is HaplotypeCalling.
 *
 * To overcome this, it is suggested that such applications may benefit from using a synthetic view of the reference which
 * collapses homopolymers longer than the ones that can be represented by flow based reads. By doing so, the artifacts of
 * this mismatch in spaces (base, flow) may be drastically reduced.
 *
 * Still, one is mostly interested in results (Haplotypes) in the original reference space, where no limit is imposed
 * on the length of a homopolymer. This gives rise to the common operation provided by this class: the 'uncollapsing' of
 * haplotypes so that long homopolymers are restored according to the reference they are aligned to.
 *
 * To perform this 'uncollapsing' the class employs several heuristics.
 *
 * Note that although the class uses "Collapsing" in its name, the main operation performed is actually
 * uncollapsing  (collapsing is a contracting operation, where uncollapsing is an expansion operation).
 * Reads are already collapsed as a result of the limit flow-based format puts on hmer lengths. It is the function
 * of this class to uncollapse such data back to a hmer length compatible with the reference.
 * 
 * To determine if a sequence (of bases) can benefit from (un)collapsing, use {@link #needsCollapsing(byte[], int, Logger)}.
 *
 * To uncollapse a set of haplotypes, use {@link #uncollapseHmersInHaplotypes(Collection, boolean, byte[])}.
 *
 */
public class LongHomopolymerHaplotypeCollapsingEngine {

    final private int             hmerSizeThreshold;
    final private boolean         partialMode;

    final private byte[]          fullRef;
    final private Locatable       refLoc;
    final private Logger          logger;
    final private boolean         debug;
    final private SmithWatermanAligner aligner;
    final private SWParameters alignmentParameters;

    // uncollapseByRef result tuple
    private static class UncollapseResult {
        byte[]      bases;
        int         offset;
        boolean     collapsed;
    }
    public LongHomopolymerHaplotypeCollapsingEngine(final int hmerSizeThreshold, final boolean partialMode, final byte[] fullRef, final Locatable refLoc, final Logger logger,
                                                    final boolean debug, final SmithWatermanAligner aligner, final SWParameters swParameters) {

        this.hmerSizeThreshold = hmerSizeThreshold;
        this.partialMode = partialMode;
        this.fullRef = fullRef;
        this.refLoc = refLoc;
        this.logger = logger;
        this.debug = debug;
        this.aligner = aligner;
        this.alignmentParameters = swParameters;
        if ( debug ) {
            logger.info("LongHomopolymerHaplotypeCollapsingEngine: >" + hmerSizeThreshold + "hmer, refLoc: " + refLoc + " fullRef:");
            logger.info(new String(fullRef));
        }
    }

    /**
     * Determine if given bases (normally a reference region) will require collapsing support for haplotypes created over it.
     *
     * Collapsing is applied when the reference contains hmers longer than a given threshold. Given flow representation,
     * such hmer will not be correctly represented in haplotypes.
     *
     * This method is called by tools, such as HaplotypeCaller, to determine if this collapsing/uncollapsing engine
     * should be deployed for a specific region.
     *
     * @param bases - bases, of a reference region, to be examined
     * @param hmerSizeThreshold - the hmer length threshold, above which, the existance of such an hmer will mark the region as needing collapsing
     * @param logger - logger to use for debug message
     * @return - boolean indicating of region related to set of the bases will benefit from collapsing
     */
    public static boolean needsCollapsing(final byte[] bases, final int hmerSizeThreshold, final Logger logger) {

        byte    lastBase = 0;
        int     baseSameCount = 0;

        if ( logger.isDebugEnabled() ) {
            logger.debug("checking for >" + hmerSizeThreshold + "hmer in:");
            logger.debug(new String(bases));
        }

        // check if has at least one sequence of stable bases larger than threshold
        for  ( byte base : bases ) {
            if ( base == lastBase ) {
                if ( ++baseSameCount >= hmerSizeThreshold ) {
                    if ( logger.isDebugEnabled() ) {
                        logger.debug("will collapse. found a stable sequence of at least " + (baseSameCount + 1) + " of " + Character.toString((char) lastBase));
                    }
                    return true;
                }
            } else {
                lastBase = base;
                baseSameCount = 0;
            }
        }

        if ( logger.isDebugEnabled() ) {
            logger.debug("will not collapse");
        }
        return false;
    }

    /**
     * Given that collapsing and uncollapsing of a set of different (by sequence) haplotypes may result in
     * duplicate haplotype, this method is a helper designed to re-index a result set of haplotypes and
     * arrange them into group of distinct haplotypes.
     *
     * @param haplotypes - input haplotypes
     * @return - a map consisting of entries in which the keys are again distinct (different from each other). The
     * Associated values are all the haplotypes of the same sequence as the key (include itself).
     */
    public static Map<Haplotype, List<Haplotype>> identicalBySequence(final List<Haplotype> haplotypes) {

        // create a map where each node's value contains a list of haplotypes with an identical sequence
        final Map<String, List<Haplotype>> sequenceMap = new LinkedHashMap<>();
        haplotypes.forEach(h -> sequenceMap.computeIfAbsent(h.getBaseString(), k -> new LinkedList<>()).add(h));

        // create a map where each node's key is a haplotype and the value is the list of all haplotypes (incl self) with same sequence
        final LinkedHashMap<Haplotype, List<Haplotype>> result = new LinkedHashMap<>();
        sequenceMap.values().forEach(h -> result.put(h.get(0), h));

        // find reference haplotype, there should always be one
        Haplotype refHaplotype = AlleleFiltering.findReferenceHaplotype(haplotypes);
        if ( refHaplotype == null )
            throw new IllegalArgumentException("Reference haplotype missing from the list of alleles");

        // make sure that result map contains the reference haplotype as one of the keys. rearrange if not
        if ( !result.containsKey(refHaplotype) ) {
            for ( Haplotype k : result.keySet()) {
                //  is the reference haplotype on this noode?
                if ( result.get(k).contains(refHaplotype) ) {

                    // switch node to be key'ed by the reference haplotype
                    result.put(refHaplotype, result.get(k));
                    result.remove(k);
                    break;
                }
            }
        }

        return result;
    }

    /**
     * This method is the main functional part of the uncollapsing engine. It takes a list of haplotypes and
     * generates a corrosponding list of new haplotypes which have their bases "uncollapse" according to the
     * given reference. For a definition of uncollapsing, please see the class definition.
     *
     * @param haplotypes                       - input haplotypes
     * @param limitToHmerSizeThreshold         - limit resulting long hmers to MmerSizeThreshold
     * @param refBasesArg                      - related reference sequence
     * @return                                 - list of uncollapsed haplotypes
     */
    public List<Haplotype> uncollapseHmersInHaplotypes(final Collection<Haplotype> haplotypes, final boolean limitToHmerSizeThreshold, final byte[] refBasesArg) {

        final List<Haplotype>       result = new LinkedList<>();
        final Map<Locatable, byte[]> refMap = new LinkedHashMap<>();
        int                         alignmentStartHapwrtRef = 0;
        byte[]                      refBases = refBasesArg;

        // locate reference haplotype, if needed, also collect refMap
        for ( Haplotype h : haplotypes ) {
            if ( h.isReference() ) {
                if ( refBases == null ) {
                    refBases = h.getBases();
                    alignmentStartHapwrtRef = h.getAlignmentStartHapwrtRef();
                }
            } else {
                byte[] ref = refMap.get(h.getGenomeLocation());
                if ( ref == null ) {
                    ref = uncollapsedPartialRef(h.getGenomeLocation());
                    refMap.put(h.getGenomeLocation(), ref);
                }
            }
        }

        // uncollapse haplotypes
        for ( Haplotype h : haplotypes ) {
            final Haplotype alignedHaplotype = uncollapseSingleHaplotype(h, limitToHmerSizeThreshold, refMap);
            alignedHaplotype.setUniquenessValue(result.size());
            result.add(alignedHaplotype);
        }

        // if we had a reference, generate cigar against it
        if ( refBases != null ) {
            for ( Haplotype h : result ) {
                if ( !h.isReference() ) {
                    final SmithWatermanAlignment alignment = aligner.align(refBases, h.getBases(),
                            alignmentParameters, SWOverhangStrategy.INDEL);
                    h.setCigar(alignment.getCigar());
                    h.setAlignmentStartHapwrtRef(alignment.getAlignmentOffset() + alignmentStartHapwrtRef);
                }
            }
        }

        return result;
    }

    private Haplotype uncollapseSingleHaplotype(final Haplotype h, final boolean limitToHmerSizeThreshold, final Map<Locatable, byte[]> refMap) {

        if ( !h.isReference() ) {

            // find ref for this location
            byte[] ref = refMap.get(h.getGenomeLocation());
            if ( ref == null ) {
                ref = uncollapsedPartialRef(h.getGenomeLocation());
                refMap.put(h.getGenomeLocation(), ref);
            }

            // uncollapse haplotype, select direction producing a stronger/longer result
            final UncollapseResult      resultFwd = uncollapseByRef(h.getBases(), ref, false);
            final UncollapseResult      resultRev = uncollapseByRef(h.getBases(), ref, true);
            final UncollapseResult      result = (resultRev.bases.length > resultFwd.bases.length) ? resultRev : resultFwd;

            // limit?
            if ( limitToHmerSizeThreshold ) {
                result.bases = collapseBases(result.bases);
            }

            // build haplotype
            Haplotype alignedHaplotype = new Haplotype(result.bases, h.isReference());
            alignedHaplotype.setScore(h.getScore());
            alignedHaplotype.setGenomeLocation(h.getGenomeLocation());
            alignedHaplotype.setEventMap(h.getEventMap());
            alignedHaplotype.setAlignmentStartHapwrtRef(result.offset);
            alignedHaplotype.setCollapsed(result.collapsed);

            return alignedHaplotype;
        } else {
            return h;
        }
    }

    private byte[] collapseBases(final byte[] fullBases) {

        // collapsed sequence would not be longer than full sequence
        final byte[]          collapsedBases = new byte[fullBases.length];

        // loop while trimming
        byte    lastBase = 0;
        int     baseSameCount = 0;
        int     dstOfs = 0;
        boolean firstHomopolymer = true;
        for  ( byte base : fullBases ) {
            if ( base == lastBase ) {
                baseSameCount++;
                if ( !firstHomopolymer && (baseSameCount >= hmerSizeThreshold) ) {
                    // collapsing, do not store
                } else {
                    // stable but under threshold, store
                    collapsedBases[dstOfs++] = base;
                }
            } else {
                // unstable, simply store
                if ( lastBase != 0 ) {
                    firstHomopolymer = false;
                }
                lastBase = base;
                baseSameCount = 0;

                collapsedBases[dstOfs++] = base;
            }
        }

        // adjust size of collapsedBases
        // not very efficient as it allocates copies the data.
        // do we really need the array to be the right size?
        return Arrays.copyOf(collapsedBases, dstOfs);
    }

    private byte[] uncollapsedPartialRef(final Locatable ucLoc) {

        final int         ucOfs = ucLoc.getStart() - refLoc.getStart();
        final int         size = ucLoc.getLengthOnReference();
        final byte[]      bases =  Arrays.copyOfRange(fullRef, ucOfs, ucOfs + size);

        if ( debug ) {
            logger.info("uncollapsedPartialRef: cOfs: " + ucOfs + ", size: " + size + ", bases:");
            logger.info(new String(bases));
        }

        return bases;
    }

    private UncollapseResult uncollapseByRef(final byte[] basesArg, final byte[] refArg, final boolean rev) {

        // prepeare result and working vars
        final UncollapseResult  uncollapseResult = new UncollapseResult();
        byte[]                  bases = basesArg;
        byte[]                  ref = refArg;

        if ( debug ) {
            logger.info("bases, ref, finalResult:");
            logger.info(new String(bases));
            logger.info(new String(ref));
        }

        if ( rev ) {
            final byte[]      basesRev = Arrays.copyOf(bases, bases.length);
            final byte[]      refRev = Arrays.copyOf(ref, ref.length);
            SequenceUtil.reverseComplement(basesRev);
            SequenceUtil.reverseComplement(refRev);
            bases = basesRev;
            ref = refRev;
        }

        // use aligner to get CIGAR
        final SmithWatermanAlignment alignment = aligner.align(ref, bases,
                alignmentParameters, SWOverhangStrategy.INDEL);
        if ( debug )
            logger.info("alignment.offset: " + alignment.getAlignmentOffset() + ", cigar: " + alignment.getCigar());

        // collect max length by walking the cigar and adding up delete operators (some of which may actually be replaced)
        int     resultLength = bases.length;
        for ( CigarElement c : alignment.getCigar() ) {
            if (c.getOperator() == CigarOperator.D)
                resultLength += c.getLength();
        }
        final byte[] result = new byte[resultLength];

        // prepare offsets
        int         basesOfs = alignment.getAlignmentOffset();
        int         refOfs = 0;
        int         resultOfs = 0;

        // start walking on cigars and make adjustments
        for ( CigarElement c : alignment.getCigar() ) {
            if (c.getOperator() != CigarOperator.D) {
                // not D - simple case
                if (c.getOperator().consumesReadBases()) {
                    System.arraycopy(bases, basesOfs, result, resultOfs, c.getLength());
                    basesOfs += c.getLength();
                    resultOfs += c.getLength();
                }
            } else {

                // check if the incoming bases contain a homopolymer
                final byte[] fwdSlice = Arrays.copyOfRange(bases, basesOfs, Math.min(basesOfs + hmerSizeThreshold, bases.length));
                final byte[] bckSlice = Arrays.copyOfRange(bases, Math.max(0, basesOfs - hmerSizeThreshold), basesOfs);
                if ( needsCollapsing(fwdSlice, hmerSizeThreshold - 1, logger) ||
                        needsCollapsing(bckSlice, hmerSizeThreshold - 1, logger) ) {


                    // check for a delete at the end of an hmer or at the beginning
                    if (onHomoPolymer(ref, refOfs - hmerSizeThreshold, ref[refOfs], hmerSizeThreshold)) {
                        // fill with base until end of homopolymer on the ref
                        final byte            base = ref[refOfs];
                        for (int size = 0; (size < c.getLength()) ; size++) {
                            if ( partialMode && ref[refOfs + size] != base )
                                break;
                            result[resultOfs++] = base;
                        }
                        uncollapseResult.collapsed = true;

                    } else if (onHomoPolymer(ref, refOfs + c.getLength(), ref[refOfs + c.getLength() - 1], hmerSizeThreshold)) {
                        // fill with base until start of homopolymer on the ref
                        final byte            base = ref[refOfs + c.getLength() - 1];
                        for (int size = 0; (size < c.getLength()) ; size++) {
                            if ( partialMode && ref[refOfs + c.getLength() - 1 - size] != base )
                                break;
                            result[resultOfs++] = base;
                        }
                        uncollapseResult.collapsed = true;
                    }
                }
            }
            if ( c.getOperator().consumesReferenceBases() )
                refOfs += c.getLength();
        }

        // return adjusted result
        final byte[] finalResult = (result.length == resultOfs) ? result : Arrays.copyOf(result, resultOfs);

        if ( rev )
            SequenceUtil.reverseComplement(finalResult);

        if ( debug ) {
            logger.info(new String(finalResult));
        }

        // return offset
        uncollapseResult.offset = alignment.getAlignmentOffset();

        // return bases
        uncollapseResult.bases = finalResult;
        return uncollapseResult;
    }

    private boolean onHomoPolymer(final byte[] bases, final int ofs, final byte base, final int length)
    {
        for ( int tick = 0 ; tick < hmerSizeThreshold ; tick++ ) {
            if ( sameBase(bases, ofs + tick, base, length) )
                return true;
        }
        return false;
    }

    private boolean sameBase(final byte[] bases, final int ofsArg, final byte base, final int lengthArg)
    {
        int     ofs = ofsArg;
        int     length = lengthArg;
        try {
            // has enough bases?
            if ( (ofs + length) > bases.length )
                return false;
            while ( length-- != 0 ) {
                if ( bases[ofs++] != base ) {
                    return false;
                }
            }
            return true;
        } catch (ArrayIndexOutOfBoundsException e) {
            return false;
        }
    }

    /**
     * Replace the haplotypes contained in the assembly result set with a fresh set.
     *
     * This method is used to replace the haplotypes in the aseembly result with uncollapsed ones.
     *
     * @param assemblyResultSet - assembly result to modify
     * @param list - fresh set of haplotypes
     */
    public void replaceAllHaplotypes(final AssemblyResultSet assemblyResultSet, final Set<Haplotype> list) {
        assemblyResultSet.clearHaplotypes();
        for ( Haplotype h : list )
            assemblyResultSet.add(h);
    }
}
