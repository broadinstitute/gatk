
package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.AssemblyRegion;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.ReadThreadingGraph;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.collections.CountSet;
import org.broadinstitute.hellbender.utils.haplotype.EventMap;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;

import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.*;

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

    private final Map<Integer,AssemblyResult> assemblyResultByKmerSize;
    private final Set<Haplotype> haplotypes;
    private final Map<Haplotype,AssemblyResult> assemblyResultByHaplotype;
    private AssemblyRegion regionForGenotyping;
    private byte[] fullReferenceWithPadding;
    private SimpleInterval paddedReferenceLoc;
    private boolean variationPresent;
    private Haplotype refHaplotype;
    private boolean wasTrimmed = false;
    private final CountSet kmerSizes;
    private SortedSet<VariantContext> variationEvents;
    private boolean debug;
    private static final Logger logger = LogManager.getLogger(AssemblyResultSet.class);

    /**
     * Constructs a new empty assembly result set.
     */
    public AssemblyResultSet() {
        assemblyResultByKmerSize = new LinkedHashMap<>(4);
        haplotypes = new LinkedHashSet<>(10);
        assemblyResultByHaplotype = new LinkedHashMap<>(10);
        kmerSizes = new CountSet(4);
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

        final Map<Haplotype,Haplotype> originalByTrimmedHaplotypes = calculateOriginalByTrimmedHaplotypes(trimmedAssemblyRegion);
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
        if (result.refHaplotype == null) {
            throw new IllegalStateException("missing reference haplotype in the trimmed set");
        }
        result.wasTrimmed = true;
        return result;
    }

    private Map<Haplotype, Haplotype> calculateOriginalByTrimmedHaplotypes(final AssemblyRegion trimmedAssemblyRegion) {
        if ( debug ) {
            logger.info("Trimming active region " + getRegionForGenotyping() + " with " + getHaplotypeCount() + " haplotypes");
        }

        final List<Haplotype> haplotypeList = getHaplotypeList();

        // trim down the haplotypes
        final Map<Haplotype, Haplotype> originalByTrimmedHaplotypes = trimDownHaplotypes(trimmedAssemblyRegion, haplotypeList);

        // create the final list of trimmed haplotypes
        final List<Haplotype> trimmedHaplotypes = new ArrayList<>(originalByTrimmedHaplotypes.keySet());

        // resort the trimmed haplotypes.
        Collections.sort(trimmedHaplotypes, Haplotype.SIZE_AND_BASE_ORDER);
        final Map<Haplotype, Haplotype> sortedOriginalByTrimmedHaplotypes = mapOriginalToTrimmed(originalByTrimmedHaplotypes, trimmedHaplotypes);

        if ( debug ) {
            logger.info("Trimmed region to " + trimmedAssemblyRegion.getSpan() + " size " +
                    trimmedAssemblyRegion.getSpan().size() + " reduced number of haplotypes from " +
                    haplotypeList.size() + " to only " + trimmedHaplotypes.size());

            for (final Haplotype remaining : trimmedHaplotypes) {
                logger.info("Remains: " + remaining + " cigar " + remaining.getCigar());
            }
        }
        return sortedOriginalByTrimmedHaplotypes;
    }

    private Map<Haplotype, Haplotype> trimDownHaplotypes(final AssemblyRegion trimmedAssemblyRegion, final List<Haplotype> haplotypeList) {
        final Map<Haplotype,Haplotype> originalByTrimmedHaplotypes = new HashMap<>();

        for ( final Haplotype h : haplotypeList ) {
            final Haplotype trimmed = h.trim(trimmedAssemblyRegion.getExtendedSpan());

            if ( trimmed != null ) {
                if (originalByTrimmedHaplotypes.containsKey(trimmed)) {
                    if (trimmed.isReference()) {
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
                        " because it starts with or ends with an insertion or deletion when trimmed to " +
                        trimmedAssemblyRegion.getExtendedSpan());
            }
        }
        return originalByTrimmedHaplotypes;
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
        pw.println("Extended Act Region " + getRegionForGenotyping().getExtendedSpan());
        pw.println("Ref haplotype coords " + getHaplotypeList().get(0).getGenomeLocation());
        pw.println("Haplotype count " + haplotypes.size());
        final Map<Integer,Integer> kmerSizeToCount = new HashMap<>();

        for (final Map.Entry<Haplotype,AssemblyResult> e : assemblyResultByHaplotype.entrySet()) {
            final AssemblyResult as = e.getValue();
            final int kmerSize = as.getGraph().getKmerSize();
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
        return kmerSizes.max();
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
        return kmerSizes.min();
    }

    /**
     * Returns a read-threading graph in the assembly set that has a particular kmerSize.
     *
     * @param kmerSize the requested kmerSize.
     *
     * @return {@code null} if there is no read-threading-graph amongst assembly results with that kmerSize.
     */
    public ReadThreadingGraph getUniqueReadThreadingGraph(final int kmerSize) {
        final AssemblyResult assemblyResult = assemblyResultByKmerSize.get(kmerSize);
        if (assemblyResult == null) {
            return null;
        }
        return assemblyResult.getThreadingGraph();
    }

    /**
     * Checks whether this assembly result set was trimmed.
     *
     * @return {@code true} iff this assembly result set was trimmed.
     */
    public boolean wasTrimmed() {
        return wasTrimmed;
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
     *
     * @return never {@code null}, but perhaps an empty collection.
     */
    public SortedSet<VariantContext> getVariationEvents() {
        if (variationEvents == null) {
            final List<Haplotype> haplotypeList = getHaplotypeList();
            EventMap.buildEventMapsForHaplotypes(haplotypeList, fullReferenceWithPadding, paddedReferenceLoc, debug);
            variationEvents = EventMap.getAllVariantContexts(haplotypeList);
        }
        return variationEvents;
    }
}
