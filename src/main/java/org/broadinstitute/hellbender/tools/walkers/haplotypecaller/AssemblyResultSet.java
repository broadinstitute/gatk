
package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.commons.lang3.mutable.MutableInt;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWOverhangStrategy;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWParameters;
import org.broadinstitute.hellbender.engine.AssemblyRegion;
import org.broadinstitute.hellbender.tools.walkers.ReferenceConfidenceVariantContextMerger;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.AbstractReadThreadingGraph;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.haplotype.Event;
import org.broadinstitute.hellbender.utils.haplotype.EventMap;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.logging.OneShotLogger;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.read.CigarUtils;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAligner;

import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Collection of read assembly using several kmerSizes.
 *
 * <p>
 *     There could be a different assembly per each kmerSize. In turn, haplotypes are result of one of those
 *     assemblies.
 * </p>
 *
 * <p>
 *     Where there is more than one possible kmerSize that generates a haplotype we consider the smaller one.
 * </p>
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.com&gt;
 */
public final class AssemblyResultSet {
    public static final Comparator<Event> HAPLOTYPE_EVENT_COMPARATOR = Comparator.comparingInt(Event::getStart)
            // Decide arbitrarily so as not to accidentally throw away overlapping variants
            .thenComparingInt(e -> e.refAllele().length())
            .thenComparing(e -> e.altAllele());
    private final Map<Integer,AssemblyResult> assemblyResultByKmerSize;
    private final Set<Haplotype> haplotypes;
    private Set<Haplotype> originalAssemblyHaps = new LinkedHashSet<>();
    private final Map<Haplotype,AssemblyResult> assemblyResultByHaplotype;
    private AssemblyRegion regionForGenotyping;
    private byte[] fullReferenceWithPadding;
    private SimpleInterval paddedReferenceLoc;
    private boolean variationPresent;
    private Haplotype refHaplotype;
    private final SortedSet<Integer>  kmerSizes;
    private SortedSet<Event> variationEvents;
    private OptionalInt lastMaxMnpDistanceUsed = OptionalInt.empty();
    private boolean debug;
    private static final Logger logger = LogManager.getLogger(AssemblyResultSet.class);
    public static final OneShotLogger haplotypeDeletionWarningLogger = new OneShotLogger(AssemblyBasedCallerUtils.class);
    private LongHomopolymerHaplotypeCollapsingEngine haplotypeCollapsingEngine; // this is nullable - indicating no collapsing engine (flow-based specific)
    private boolean isPartiallyDeterminedList = false;

    /**
     * Constructs a new empty assembly result set.
     */
    public AssemblyResultSet() {
        assemblyResultByKmerSize = new LinkedHashMap<>(4);
        haplotypes = new LinkedHashSet<>(10);
        assemblyResultByHaplotype = new LinkedHashMap<>(10);
        kmerSizes = new TreeSet<>();
    }

    /**
     * Trims an assembly result set down based on a new set of trimmed haplotypes.
     *
     * @param trimmedAssemblyRegion the trimmed down active region.
     *
     * @throws NullPointerException if any argument in {@code null} or
     *      if there are {@code null} entries in {@code originalByTrimmedHaplotypes} for trimmed haplotype keys.
     * @throws IllegalArgumentException if there is no reference haplotype amongst the trimmed ones.
     *
     * @return never {@code null}, a new trimmed assembly result set.
     */
    public AssemblyResultSet trimTo(final AssemblyRegion trimmedAssemblyRegion) {

        final Map<Haplotype,Haplotype> originalByTrimmedHaplotypes = calculateOriginalByTrimmedHaplotypes(trimmedAssemblyRegion.getPaddedSpan());
        if (refHaplotype == null) {
            throw new IllegalStateException("refHaplotype is null");
        }
        Utils.nonNull(trimmedAssemblyRegion);
        final AssemblyResultSet result = new AssemblyResultSet();

        for (final Haplotype trimmed : originalByTrimmedHaplotypes.keySet()) {
            final Haplotype original = originalByTrimmedHaplotypes.get(trimmed);
            if (original == null) {
                throw new IllegalStateException("all trimmed haplotypes must have an original one");
            }
            final AssemblyResult as = assemblyResultByHaplotype.get(original);
            if (as == null) {
                result.add(trimmed);
            } else {
                result.add(trimmed, as);
            }
        }

        result.setRegionForGenotyping(trimmedAssemblyRegion);
        result.setFullReferenceWithPadding(fullReferenceWithPadding);
        result.setPaddedReferenceLoc(paddedReferenceLoc);
        result.variationPresent = haplotypes.stream().anyMatch(Haplotype::isNonReference);
        if (result.refHaplotype == null) {
            throw new IllegalStateException("missing reference haplotype in the trimmed set");
        }
        return result;
    }

    private Map<Haplotype, Haplotype> calculateOriginalByTrimmedHaplotypes(final Locatable span) {


        if ( debug ) {
            logger.info("Trimming active region " + getRegionForGenotyping() + " with " + getHaplotypeCount() + " haplotypes");
        }

        final List<Haplotype> haplotypeList = getHaplotypeList();

        // trim down the haplotypes
        final Map<Haplotype, Haplotype> originalByTrimmedHaplotypes = trimDownHaplotypes(span, haplotypeList);

        // create the final list of trimmed haplotypes
        final List<Haplotype> trimmedHaplotypes = new ArrayList<>(originalByTrimmedHaplotypes.keySet());

        // resort the trimmed haplotypes.
        Collections.sort(trimmedHaplotypes, Haplotype.SIZE_AND_BASE_ORDER);
        final Map<Haplotype, Haplotype> sortedOriginalByTrimmedHaplotypes = mapOriginalToTrimmed(originalByTrimmedHaplotypes, trimmedHaplotypes);

        if ( debug ) {
            logger.info("Trimmed region to " + span + " and reduced number of haplotypes from " +
                    haplotypeList.size() + " to only " + trimmedHaplotypes.size());

            for (final Haplotype remaining : trimmedHaplotypes) {
                logger.info("Remains: " + remaining + " cigar " + remaining.getCigar());
            }
        }
        return sortedOriginalByTrimmedHaplotypes;
    }

    // trim haplotypes to span merging haplotypes that are equal in bases
    private Map<Haplotype, Haplotype> trimDownHaplotypes(final Locatable span, final List<Haplotype> haplotypeList) {

        final Map<Haplotype,Haplotype> originalByTrimmedHaplotypes = new HashMap<>();

        // Note that two haplotypes one of which is marked as reference are different,
        // so trivially after trimming duplicate haplotypes could occur. We remove
        // these duplcations. If one of the haplotypes before trimming was marked as
        // reference - we mark the trimmed haplotype reference

        for ( final Haplotype h : haplotypeList ) {
            final Haplotype trimmed = h.trim(span, true);

            if ( trimmed != null ) {
                if (originalByTrimmedHaplotypes.containsKey(trimmed)) {
                    if (h.isReference()) {
                        originalByTrimmedHaplotypes.remove(trimmed);
                        originalByTrimmedHaplotypes.put(trimmed, h);
                    }
                } else {
                    originalByTrimmedHaplotypes.put(trimmed, h);
                }
            } else if (h.isReference()) {
                throw new IllegalStateException("trimming eliminates the reference haplotype");
            } else if ( debug ) {
                logger.info("Throwing out haplotype " + h + " with cigar " + h.getCigar() +
                        " because it starts with or ends with an insertion or deletion when trimmed to " + span);
            }
        }

        // Now set reference status originalByTrimmedHaplotypes
        final Map<Haplotype, Haplotype> fixedOriginalByTrimmedHaplotypes = new HashMap<>();
        for (Haplotype h : originalByTrimmedHaplotypes.keySet()){
            if (originalByTrimmedHaplotypes.get(h).isReference()){
                Haplotype fixedHap = new Haplotype(h.getBases(), true);
                fixedHap.setCigar(h.getCigar());
                fixedHap.setGenomeLocation(h.getGenomeLocation());
                fixedHap.setScore(h.getScore());
                fixedHap.setAlignmentStartHapwrtRef(h.getAlignmentStartHapwrtRef());
                fixedOriginalByTrimmedHaplotypes.put(fixedHap, originalByTrimmedHaplotypes.get(h));
            } else {
                fixedOriginalByTrimmedHaplotypes.put(h, originalByTrimmedHaplotypes.get(h));
            }
        }
        return fixedOriginalByTrimmedHaplotypes;
    }

    private static Map<Haplotype, Haplotype> mapOriginalToTrimmed(final Map<Haplotype, Haplotype> originalByTrimmedHaplotypes, final List<Haplotype> trimmedHaplotypes) {
        final Map<Haplotype,Haplotype> sortedOriginalByTrimmedHaplotypes = new LinkedHashMap<>(trimmedHaplotypes.size());
        for (final Haplotype trimmed : trimmedHaplotypes) {
            sortedOriginalByTrimmedHaplotypes.put(trimmed, originalByTrimmedHaplotypes.get(trimmed));
        }
        return sortedOriginalByTrimmedHaplotypes;
    }

    /**
     * Query the reference haplotype in the result set.
     * @return {@code null} if none wasn't yet added, otherwise a reference haplotype.
     */
    public Haplotype getReferenceHaplotype() {
        return refHaplotype;
    }

    /**
     * Checks whether there is any variation present in the assembly result set.
     *
     * <p>
     *     This is equivalent to whether there is more than one haplotype.
     * </p>
     *
     * @return {@code true} if there is variation present, {@code false} otherwise.
     */
    public boolean isVariationPresent() {
        return variationPresent && haplotypes.size() > 1;
    }

    /**
     * Dumps debugging information into a print-writer.
     *
     * @param pw where to dump the information.
     *
     * @throws NullPointerException if {@code pw} is {@code null}.
     */
    private void debugDump(final PrintWriter pw) {
        if (getHaplotypeList().isEmpty()) {
            return;
        }
        pw.println("Active Region " + regionForGenotyping.getSpan());
        pw.println("Extended Act Region " + getRegionForGenotyping().getPaddedSpan());
        pw.println("Ref haplotype coords " + getHaplotypeList().get(0).getGenomeLocation());
        pw.println("Haplotype count " + haplotypes.size());
        final Map<Integer,Integer> kmerSizeToCount = new HashMap<>();

        for (final Map.Entry<Haplotype,AssemblyResult> e : assemblyResultByHaplotype.entrySet()) {
            final AssemblyResult as = e.getValue();
            final int kmerSize = as.getKmerSize();
            if (kmerSizeToCount.containsKey(kmerSize)) {
                kmerSizeToCount.put(kmerSize,kmerSizeToCount.get(kmerSize) + 1);
            } else {
                kmerSizeToCount.put(kmerSize,1);
            }
        }
        pw.println("Kmer sizes count " + kmerSizeToCount.entrySet().size());
        final Integer[] kmerSizes = kmerSizeToCount.keySet().toArray(new Integer[kmerSizeToCount.size()]);
        Arrays.sort(kmerSizes);
        pw.println("Kmer sizes values " + Arrays.toString(kmerSizes));
        for (final int size : kmerSizes) {
            pw.println("Kmer size " + size + " count " + kmerSizeToCount.get(size));
        }
    }

    /**
     * Adds a haplotype to the result set without indicating a generating assembly result.
     *
     * <p>
     *     It is possible to call this method with the same haplotype several times. In that the second and further
     *     calls won't have any effect (thus returning {@code false}).
     * </p>
     *
     * @param h the haplotype to add to the assembly result set.
     *
     * @throws NullPointerException if {@code h} is {@code null}
     * @throws IllegalArgumentException if {@code h} does not have a genome location.
     *
     * @return {@code true} if the assembly result set has been modified as a result of this call.
     */
    public boolean add(final Haplotype h) {

        Utils.nonNull(h, "input haplotype cannot be null");
        Utils.nonNull(h.getGenomeLocation(), "haplotype genomeLocation cannot be null");
        if (haplotypes.contains(h)) {
            return false;
        }
        haplotypes.add(h);
        updateReferenceHaplotype(h);
        return true;
    }

    /**
     * Adds simultaneously a haplotype and the generating assembly-result.
     *
     * <p>
     *     Haplotypes and their assembly-result can be added multiple times although just the first call will have
     *     any effect (return value is {@code true}).
     * </p>
     *
     *
     * @param h haplotype to add.
     * @param ar assembly-result that is assumed to have given rise to that haplotype.
     *
     * @throws NullPointerException if {@code h} or {@code ar} is {@code null}.
     * @throws IllegalArgumentException if {@code h} has not defined genome location.
     *
     * @return {@code true} iff this called changes the assembly result set.
     */
    public boolean add(final Haplotype h, final AssemblyResult ar) {
        Utils.nonNull(h, "input haplotype cannot be null");
        Utils.nonNull(ar, "input assembly-result cannot be null");
        Utils.nonNull(h.getGenomeLocation(), "the haplotype provided must have a genomic location");

        final boolean assemblyResultAdditionReturn =  add(ar);

        if (haplotypes.contains(h)) {
            final AssemblyResult previousAr = assemblyResultByHaplotype.get(h);
            if (previousAr == null) {
                assemblyResultByHaplotype.put(h, ar);
                return true;
            } else if (!previousAr.equals(ar)) {
                throw new IllegalStateException("there is already a different assembly result for the input haplotype");
            } else {
                return assemblyResultAdditionReturn;
            }
        } else {
            haplotypes.add(h);
            assemblyResultByHaplotype.put(h,ar);
            updateReferenceHaplotype(h);
            if (h.isNonReference()) {
                variationPresent = true;
            }
            return true;
        }
    }

    /**
     * Add a assembly-result object.
     *
     * @param ar the assembly result to add.
     *
     * @throws NullPointerException if {@code ar} is {@code null}.
     * @throws IllegalStateException if there is an assembly result with the same kmerSize.
     * @return {@code true} iff this addition changed the assembly result set.
     */
    private boolean add(final AssemblyResult ar) {
        Utils.nonNull(ar);
        final int kmerSize = ar.getKmerSize();
        if (assemblyResultByKmerSize.containsKey(kmerSize)) {
            if (!assemblyResultByKmerSize.get(kmerSize).equals(ar)) {
                throw new IllegalStateException("a different assembly result with the same kmerSize was already added");
            }
            return false;
        } else {
            assemblyResultByKmerSize.put(kmerSize, ar);
            kmerSizes.add(kmerSize);
            return true;
        }
    }

    /**
     * Returns the current region for genotyping.
     *
     * @return might be {@code null}.
     */
    public AssemblyRegion getRegionForGenotyping() {
        return regionForGenotyping;
    }

    /**
     * Sets the region for genotyping.
     *
     * @param regionForGenotyping the new value.
     */
    public void setRegionForGenotyping(final AssemblyRegion regionForGenotyping) {
        this.regionForGenotyping = regionForGenotyping;
    }

    /**
     * Returns the current full reference with padding.
     *
     * @return might be {@code null}. The result must not be modified by the caller.
     */
    public byte[] getFullReferenceWithPadding() {
        return fullReferenceWithPadding;
    }

    /**
     * Sets the full reference with padding base sequence.
     *
     * @param fullReferenceWithPadding the new value. The array must not be modified by the caller.
     */
    public void setFullReferenceWithPadding(final byte[] fullReferenceWithPadding) {
        this.fullReferenceWithPadding = fullReferenceWithPadding;
    }

    /**
     * Returns the padded reference location.
     *
     * @return might be {@code null}
     */
    public SimpleInterval getPaddedReferenceLoc() {
        return paddedReferenceLoc;
    }

    /**
     * Changes the padded reference location.
     * @param paddedReferenceLoc the new value.
     */
    public void setPaddedReferenceLoc(final SimpleInterval paddedReferenceLoc) {
        this.paddedReferenceLoc = paddedReferenceLoc;
    }

    /**
     * Returns the number of haplotypes in the assembly result set.
     * @return {@code 0} or greater.
     */
    public int getHaplotypeCount() {
        return haplotypes.size();
    }

    /**
     * Returns the haplotypes as a list.
     *
     * <p>
     *     The result is unmodifiable.
     * </p>
     *
     * @return never {@code null}, but perhaps a empty list if no haplotype was generated during assembly.
     */
    public List<Haplotype> getHaplotypeList() {
        return Arrays.asList(haplotypes.toArray(new Haplotype[haplotypes.size()]));
    }

    /**
     * Returns the maximum kmerSize available.
     *
     * @throws IllegalStateException if no assembly-result was added to the set, thus there is no kmerSize.
     *
     * @return greater than 0.
     */
    public int getMaximumKmerSize() {
        if (kmerSizes.isEmpty()) {
            throw new IllegalStateException("there is yet no kmerSize in this assembly result set");
        }
        return kmerSizes.last();
    }

    /**
     * Indicates whether there are more than one kmerSize in the set.
     *
     * @return {@code true} iff there is more than one kmerSize assembly in the set.
     */
    public boolean hasMultipleKmerSizes() {
        return kmerSizes.size() > 1;
    }

    /**
     * Returns the minimum kmerSize available.
     *
     * @throws IllegalStateException if no assembly-result was added to the set, thus there is no kmerSize.
     *
     * @return greater than 0.
     */
    public int getMinimumKmerSize() {
        if (kmerSizes.isEmpty()) {
            throw new IllegalStateException("there is yet no kmerSize in this assembly result set");
        }
        return kmerSizes.first();
    }

    /**
     * Returns a read-threading graph in the assembly set that has a particular kmerSize.
     *
     * @param kmerSize the requested kmerSize.
     *
     * @return {@code null} if there is no read-threading-graph amongst assembly results with that kmerSize.
     */
    public AbstractReadThreadingGraph getUniqueReadThreadingGraph(final int kmerSize) {
        final AssemblyResult assemblyResult = assemblyResultByKmerSize.get(kmerSize);
        if (assemblyResult == null) {
            return null;
        }
        return assemblyResult.getThreadingGraph();
    }

    /**
     * Dumps debugging information into a logger.
     *
     * @param logger where to dump the information.
     *
     * @throws NullPointerException if {@code logger} is {@code null}.
     */
    public void debugDump(final Logger logger) {
        final StringWriter sw = new StringWriter();
        final PrintWriter pw = new PrintWriter(sw);
        debugDump(pw);
        final String str = sw.toString();
        final String[] lines = str.split("\n");
        for (final String line : lines) {
            if (line.isEmpty()) {
                continue;
            }
            logger.debug(line);
        }
    }

    /**
     * Given whether a new haplotype that has been already added to {@link #haplotypes} collection is the
     * reference haplotype and updates {@link #refHaplotype} accordingly.
     *
     * <p>
     *     This method assumes that the colling code has verified that the haplotype was not already in {@link #haplotypes}
     *     I.e. that it is really a new one. Otherwise it will result in an exception if it happen to be a reference
     *     haplotype and this has already be set. This is the case even if the new haplotypes and the current reference
     *     are equal.
     * </p>
     *
     * @param newHaplotype the new haplotype.
     * @throws NullPointerException if {@code newHaplotype} is {@code null}.
     * @throws IllegalStateException if there is already a reference haplotype.
     */
    private void updateReferenceHaplotype(final Haplotype newHaplotype) {
        if (!newHaplotype.isReference()) {
            return;
        }
        if (refHaplotype == null) {
            refHaplotype = newHaplotype;
        } else {// assumes that we have checked wether the haplotype is already in the collection and so is no need to check equality.
            throw new IllegalStateException("the assembly-result-set already have a reference haplotype that is different");
        }
    }

    /**
     * Returns a sorted set of variant events that best explain the haplotypes found by the assembly
     * across kmerSizes.
     *
     * <p/>
     * The result is sorted incrementally by location.
     * @param maxMnpDistance Phased substitutions separated by this distance or less are merged into MNPs.  More than
     *                       two substitutions occuring in the same alignment block (ie the same M/X/EQ CIGAR element)
     *                       are merged until a substitution is separated from the previous one by a greater distance.
     *                       That is, if maxMnpDistance = 1, substitutions at 10,11,12,14,15,17 are partitioned into a MNP
     *                       at 10-12, a MNP at 14-15, and a SNP at 17.  May not be negative.
     * @return never {@code null}, but perhaps an empty collection.
     */
    public SortedSet<Event> getVariationEvents(final int maxMnpDistance) {
        ParamUtils.isPositiveOrZero(maxMnpDistance, "maxMnpDistance may not be negative.");

        final boolean sameMnpDistance = lastMaxMnpDistanceUsed.isPresent() && maxMnpDistance == lastMaxMnpDistanceUsed.getAsInt();
        lastMaxMnpDistanceUsed = OptionalInt.of(maxMnpDistance);

        if (variationEvents == null || !sameMnpDistance || haplotypes.stream().anyMatch(hap -> hap.isNonReference() && hap.getEventMap() == null)) {
            regenerateVariationEvents(maxMnpDistance);
        }
        return variationEvents;
    }

    public void regenerateVariationEvents(int maxMnpDistance) {
        final List<Haplotype> haplotypeList = getHaplotypeList();
        EventMap.buildEventMapsForHaplotypes(haplotypeList, fullReferenceWithPadding, paddedReferenceLoc, debug, maxMnpDistance);
        variationEvents = getAllVariantContexts(haplotypeList);
        lastMaxMnpDistanceUsed = OptionalInt.of(maxMnpDistance);
        variationPresent = haplotypeList.stream().anyMatch(Haplotype::isNonReference);
    }

    /**
     * Get all of the VariantContexts in the event maps for all haplotypes, sorted by their start position and then arbitrarily by indel length followed by bases
     * @param haplotypes the set of haplotypes to grab the VCs from
     * @return a sorted set of variant contexts
     */
    private static SortedSet<Event> getAllVariantContexts(final List<Haplotype> haplotypes ) {
        return haplotypes.stream().flatMap(h -> h.getEventMap().getEvents().stream())
                .collect(Collectors.toCollection(() -> new TreeSet<>(HAPLOTYPE_EVENT_COMPARATOR)));
    }

    public void setDebug(boolean debug) {
        this.debug = debug;
    }

    public LongHomopolymerHaplotypeCollapsingEngine getHaplotypeCollapsingEngine() {
        return haplotypeCollapsingEngine;
    }

    public void setHaplotypeCollapsingEngine(LongHomopolymerHaplotypeCollapsingEngine haplotypeCollapsingEngine) {
        this.haplotypeCollapsingEngine = haplotypeCollapsingEngine;
    }

    public void clearHaplotypes() {
        haplotypes.clear();
        refHaplotype = null;
    }
    public void replaceAllHaplotypes(Collection<Haplotype> list) {
        haplotypes.clear();
        refHaplotype = null;
        for ( Haplotype h : list ) {
            add(h);
            if (h.isNonReference()) {
                variationPresent = true;
            }
        }
    }

    // For PDHMM use: Remove a haplotype from this AssemblyResultSet and update all of the various references in the
    //                object to that haplotype to be current.
    // WARNING: Deleting haplotypes in this manner is highly dangerous and will likely lead to lost variants
    private void removeHaplotype(final Haplotype hap) {
        haplotypes.remove(hap);
        assemblyResultByHaplotype.remove(hap);
        for (Integer kmerSize : assemblyResultByKmerSize.keySet()) {
            Set<Haplotype> discovered = assemblyResultByKmerSize.get(kmerSize).getDiscoveredHaplotypes();
            discovered.remove(hap);
            assemblyResultByKmerSize.get(kmerSize).setDiscoveredHaplotypes(discovered);
        }
    }

    public void setPartiallyDeterminedMode(final boolean isPartiallyDetermined) {
        this.isPartiallyDeterminedList = isPartiallyDetermined;
    }

    public boolean isPartiallyDeterminedList() {
        return isPartiallyDeterminedList;
    }

    // For PDHMM use: store original assembly haplotypes if we are injecting artificial haplotypes later
    // WARNING: This is likely not set in every case, use with caution
    public void storeAssemblyHaplotypes() {
        originalAssemblyHaps = new LinkedHashSet<>(haplotypes);
    }

    public boolean hasOverwrittenHaps() {
        return !originalAssemblyHaps.isEmpty();
    }

    /**
     * Remove haplotypes with alleles we wish to filter
     * // TODO this is a bad algorithm for bad people -- it might eliminate good alleles on the same haplotypes
     */
    public void removeHaplotypesWithBadAlleles(final AssemblyBasedCallerArgumentCollection argumentCollection,
                                               final Collection<Event> badPileupEvents) {
        if (badPileupEvents.isEmpty()) {
            return; // nothing to do
        }
        List<Haplotype> haplotypesWithFilterAlleles = new ArrayList<>();
        // If we are not using partially determined haplotypes, we discard every haplotype containing a bad pileup allele.


        for(Event badEvent : badPileupEvents) {
            for (Haplotype hap : getHaplotypeList()) {
                // NOTE: The event map may be null due to edge cases in the assembly engine + SW realignment that can cause
                //       haplotypes to re-merge into the reference and end up with an empty event map. (Also the event map
                //       code explicitly expects to fail in some instances)
                if (hap.getEventMap() != null && hap.getEventMap().getEvents().stream().anyMatch(badEvent::equals)) {
                    if (argumentCollection.pileupDetectionArgs.debugPileupStdout) System.err.println("Flagging hap " + hap + " for containing variant " + badEvent);
                    haplotypesWithFilterAlleles.add(hap);
                }
            }
        }

        if (!haplotypesWithFilterAlleles.isEmpty()) {
            if (argumentCollection.pileupDetectionArgs.debugPileupStdout) System.out.println("Found Assembly Haps with filtered Variants:\n"+haplotypesWithFilterAlleles.stream().map(Haplotype::toString).collect(Collectors.joining("\n")));
            haplotypeDeletionWarningLogger.warn(() -> "Haplotypes from Assembly are being filtered by heuristics from the PileupCaller. This can lead to missing variants. See --"+PileupDetectionArgumentCollection.PILEUP_DETECTION_FILTER_ASSEMBLY_HAPS_THRESHOLD+" for more details");
            haplotypesWithFilterAlleles.forEach(this::removeHaplotype);
        }
    }

    /**
     * Helper method that handles the actual "GGA-like" Merging of haplotype alleles into an assembly result set.
     *
     * First this method will filter out haplotypes that contain alleles that have failed the pileup calling filtering steps,
     * Then the list will attempt to poke into the haplotype list artificial haplotypes that have the found alleles present.
     */
    @SuppressWarnings({"deprecation"})
    public void injectPileupEvents(final AssemblyRegion region, final AssemblyBasedCallerArgumentCollection argumentCollection,
                                   final SmithWatermanAligner aligner, final Collection<Event> goodPileupEvents) {
        if (goodPileupEvents.isEmpty()) {
            return; // nothing to do
        }

        if (debug) {
            logger.info("Number of haplotypes pre-pileup injection: " + this.getHaplotypeCount());
        }

        final Haplotype refHaplotype = getReferenceHaplotype();
        final Map<Integer, List<Event>> assembledEventByStart = getVariationEvents(argumentCollection.maxMnpDistance).stream()
                .collect(Collectors.groupingBy(Event::getStart));
        final Collection<Event> assembledIndels = getVariationEvents(argumentCollection.maxMnpDistance).stream().
                filter(Event::isIndel).toList();

        Set<Haplotype> baseHaplotypes = new TreeSet<>();
        baseHaplotypes.addAll(getHaplotypeList().stream()
                .sorted(Comparator.comparingInt((Haplotype hap) -> hap.isReference() ? 1 : 0).thenComparingDouble(Haplotype::getScore).reversed())
                .limit(AssemblyBasedCallerUtils.NUM_HAPLOTYPES_TO_INJECT_FORCE_CALLING_ALLELES_INTO)
                .toList());

        //TODO its unclear whether the correct answer here is to use the hardclipped pileup reads (which we used in generating the pileup alleles for specificty reasons)
        //TODO or if it would be more accurate to use the softclipped bases here in filtering down the haplotypes. I suspect the latter but I will evaluate later.
        Map<Kmer, MutableInt> kmerReadCounts = AssemblyBasedCallerUtils.getKmerReadCounts(region.getHardClippedPileupReads(), argumentCollection.pileupDetectionArgs.filteringKmerSize);

        for (final Event event : goodPileupEvents.stream().sorted(Comparator.comparingInt(Event::getStart)).toList()) {

            if (argumentCollection.pileupDetectionArgs.debugPileupStdout) System.out.println("Processing new Haplotypes for Pileup Allele that was not in the assembly: " + event);

            // skip SNPs that are too close to assembled indels.
            if (!event.isIndel() && assembledIndels.stream().anyMatch(indel -> event.withinDistanceOf(indel, argumentCollection.pileupDetectionArgs.snpAdjacentToAssemblyIndel))) {
                continue;
            }
            final List<Event> assembledEvents = assembledEventByStart.getOrDefault(event.getStart(), Collections.emptyList());

            if (isEventPresentInAssembly(event, assembledEvents) || isSymbolic(event.altAllele())) {
                continue;
            }

            final Set<Haplotype> newPileupHaplotypes = new LinkedHashSet<>();
            for (final Haplotype baseHaplotype : baseHaplotypes) {
                final Haplotype insertedHaplotype = makeHaplotypeWithInsertedEvent(baseHaplotype, refHaplotype, event, aligner, argumentCollection.getHaplotypeToReferenceSWParameters());
                if (insertedHaplotype != null) {
                    newPileupHaplotypes.add(insertedHaplotype);
                }
            }

            final Set<Haplotype> refactoredHaps = AssemblyBasedCallerUtils.filterPileupHaplotypes(newPileupHaplotypes, kmerReadCounts,
                    argumentCollection.pileupDetectionArgs.numHaplotypesToIterate, argumentCollection.pileupDetectionArgs.filteringKmerSize);
            if (argumentCollection.pileupDetectionArgs.debugPileupStdout) {
                System.out.println("Constructed the following new Pileup Haplotypes after filtering:\n"+
                        refactoredHaps.stream().map(Haplotype::toString).collect(Collectors.joining("\n")));
            }
            
            baseHaplotypes.addAll(refactoredHaps);
        }


        baseHaplotypes.forEach(this::add);
        regenerateVariationEvents(argumentCollection.maxMnpDistance);
    }

    @VisibleForTesting
    public void addGivenAlleles(final List<Event> givenAlleles, final int maxMnpDistance, final SmithWatermanAligner aligner,
                                final SWParameters haplotypeToReferenceSWParameters) {
        if (givenAlleles.isEmpty()) {
            return;
        }
        final Haplotype refHaplotype = getReferenceHaplotype();
        final Map<Integer, List<Event>> assembledEventsByStart = getVariationEvents(maxMnpDistance).stream()
                .collect(Collectors.groupingBy(Event::getStart));

        // choose the highest-scoring haplotypes along with the reference for building force-calling haplotypes
        final List<Haplotype> baseHaplotypes = getHaplotypeList().stream()
                .sorted(Comparator.comparing(Haplotype::isReference).thenComparingDouble(hap -> hap.getScore()).reversed())
                .limit(AssemblyBasedCallerUtils.NUM_HAPLOTYPES_TO_INJECT_FORCE_CALLING_ALLELES_INTO)
                .collect(Collectors.toList());

        for (final Event given : givenAlleles) {
            final List<Event> assembledEvents = assembledEventsByStart.getOrDefault(given.getStart(), Collections.emptyList());

            if (isEventPresentInAssembly(given, assembledEvents) || isSymbolic(given.altAllele())) {
                continue;
            }

            for (final Haplotype baseHaplotype : baseHaplotypes) {
                final Haplotype insertedHaplotype = makeHaplotypeWithInsertedEvent(baseHaplotype, refHaplotype, given, aligner, haplotypeToReferenceSWParameters);
                if (insertedHaplotype != null) {
                    add(insertedHaplotype);
                }
            }
        }
        regenerateVariationEvents(maxMnpDistance);
    }

    private static boolean isEventPresentInAssembly(final Event event, final List<Event> assembledEvents) {
        return assembledEvents.stream().anyMatch(event::equals);    // note that Events are forced to have a minimal representation
    }

    private static boolean isSymbolic(final Allele a) {
        return a.equals(Allele.NO_CALL) || a.getDisplayString().equals(String.valueOf(VCFConstants.NULL_ALLELE)) || a.equals(Allele.SPAN_DEL) || a.isSymbolic();
    }

    // returns a new haplotype with the event inserted or null if that is not possible
    private static Haplotype makeHaplotypeWithInsertedEvent(final Haplotype baseHaplotype, final Haplotype refHaplotype, final Event event,
                                                            final SmithWatermanAligner aligner, final SWParameters haplotypeToReferenceSWParameters) {
        if (baseHaplotype.getEventMap() != null && baseHaplotype.getEventMap().getEvents().stream().anyMatch(event::overlaps)) {
            return null; // Can't insert because the base haplotype already has an event there.
        }

        final Haplotype insertedHaplotype = baseHaplotype.insertAllele(event.refAllele(), event.altAllele(), event.getStart());
        if (insertedHaplotype != null) { // can be null if the requested allele can't be inserted into the haplotype
            final Cigar cigar = CigarUtils.calculateCigar(refHaplotype.getBases(), insertedHaplotype.getBases(), aligner, haplotypeToReferenceSWParameters, SWOverhangStrategy.INDEL);
            insertedHaplotype.setCigar(cigar);
            insertedHaplotype.setGenomeLocation(refHaplotype.getGenomeLocation());
            insertedHaplotype.setAlignmentStartHapwrtRef(refHaplotype.getAlignmentStartHapwrtRef());

        }
        return insertedHaplotype;
    }
}
