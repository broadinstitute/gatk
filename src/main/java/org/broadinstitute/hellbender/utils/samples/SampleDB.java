package org.broadinstitute.hellbender.utils.samples;

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.variant.variantcontext.Genotype;

import org.broadinstitute.hellbender.exceptions.GATKException;

import java.util.*;

/**
 * Simple database for managing samples
 */
public class SampleDB {
    /**
     * This is where Sample objects are stored. Samples are usually accessed by their ID, which is unique, so
     * this is stored as a HashMap.
     */
    private final HashMap<String, Sample> samples = new HashMap<>();

    /**
     * Package private for use by SampleDBBuilder .
     */
    SampleDB() {}

    /**
     * Protected function to add a single sample to the database
     *
     * @param newSample to be added
     */
    protected SampleDB addSample(final Sample newSample) {
        Sample updatedSample = newSample;

        Sample prevSample = samples.get(newSample.getID());
        if (prevSample != null) {
            updatedSample = prevSample.mergeSamples(newSample);
        }
        samples.put(newSample.getID(), updatedSample);
        return this;
    }

    // --------------------------------------------------------------------------------
    //
    // Functions for getting a sample from the DB
    //
    // --------------------------------------------------------------------------------

    /**
     * Get a sample by its ID
     * If an alias is passed in, return the main sample object
     * @param id
     * @return sample Object with this ID, or null if this does not exist
     */
    public Sample getSample(String id) {
        return samples.get(id);
    }

    // --------------------------------------------------------------------------------
    //
    // Functions for accessing samples in the DB
    //
    // --------------------------------------------------------------------------------

    /**
     * Get number of sample objects
     * @return size of samples map
     */
    public int sampleCount() {
        return samples.size();
    }

    public Set<Sample> getSamples() {
        return new HashSet<>(samples.values());
    }

    public Collection<String> getSampleNames() {
        return Collections.unmodifiableCollection(samples.keySet());
    }


    /**
     * Takes a collection of sample names and returns their corresponding sample objects
     * Note that, since a set is returned, if you pass in a list with duplicates names there will not be any duplicates in the returned set
     * @param sampleNameList Set of sample names
     * @return Corresponding set of samples
     */
    public Set<Sample> getSamples(Collection<String> sampleNameList) {
        HashSet<Sample> samples = new HashSet<>();
        for (String name : sampleNameList) {
            try {
                samples.add(getSample(name));
            }
            catch (Exception e) {
                throw new GATKException("Could not get sample with the following ID: " + name, e);
            }
        }
        return samples;
    }

    // --------------------------------------------------------------------------------
    //
    // Higher level pedigree functions
    //
    // --------------------------------------------------------------------------------

    /**
     * Returns a sorted set of the family IDs in all samples
     * @return Sorted set of the family IDs in all samples (excluding null ids)
     */
    public final Set<String> getFamilyIDs() {
        return getFamilies().keySet();
    }

    /**
     * Returns a map from family ID -> set of family members.
     * @return Map from family ID -> set of family members for all samples with non-null family ids
     */
    public final Map<String, Set<Sample>> getFamilies() {
        return getFamilies(null);
    }

    /**
     * @param sampleIds - all samples to include. If null is passed then all samples are returned.
     * @return Map from family ID -> set of family members for all samples in sampleIds with
     * non-null family ids
     */
    public final Map<String, Set<Sample>> getFamilies(Collection<String> sampleIds) {
        final Map<String, Set<Sample>> families = new TreeMap<>();

        for (final Sample sample : samples.values()) {
            if (sampleIds == null || sampleIds.contains(sample.getID())) {
                final String famID = sample.getFamilyID();
                if (famID != null) {
                    if (!families.containsKey(famID)) {
                        families.put(famID, new TreeSet<>());
                    }
                    families.get(famID).add(sample);
                }
            }
        }
        return families;
    }

    /**
     * Return all samples with a given family ID
     * @param familyId
     * @return Set of all samples with the given family id.
     */
    public Set<Sample> getFamily(String familyId) {
        return getFamilies().get(familyId);
    }

    /**
     * Returns children of a given sample.
     * @param sample
     * @return All children of a given sample. See note on the efficiency of getFamily() - since this depends on
     * getFamily() it's also not efficient
     */
    public Set<Sample> getChildren(Sample sample) {
        final HashSet<Sample> children = new HashSet<>();
        for (final Sample familyMember : getFamily(sample.getFamilyID())) {
            if (getMother(familyMember) == sample || getFather(familyMember) == sample) {
                children.add(familyMember);
            }
        }
        return children;
    }

    public Set<String> getFounderIds(){
        Set<String> founders = new HashSet<>();
        for (Sample sample : getSamples()) {
            if (getParents(sample).size() < 1) {
                founders.add(sample.getID());
            }
        }
        return founders;
    }

    /**
     * Get the sample's mother
     * @param offSpring child of mother to return
     * @return sample object with relationship mother, if exists, or null
     */
    public Sample getMother(Sample offSpring) {
        String maternalID = offSpring.getMaternalID();
        return null == maternalID ? null : samples.get(maternalID);
    }

    /**
     * Get the sample's father
     * @param offSpring child of father to return
     * @return sample object with relationship father, if exists, or null
     */
    public Sample getFather(Sample offSpring) {
        String paternalID = offSpring.getPaternalID();
        return null == paternalID ? null : samples.get(paternalID);
    }

    /**
     * Get the sample's father and mother
     * @param offSpring child of parents to return
     * @return sample objects with relationship parents, if exists, or null
     */
    public ArrayList<Sample> getParents(Sample offSpring) {
        ArrayList<Sample> parents = new ArrayList<>(2);
        Sample parent = getMother(offSpring);
        if (parent != null)
            parents.add(parent);
        parent = getFather(offSpring);
        if(parent != null)
            parents.add(parent);
        return parents;
    }

}
