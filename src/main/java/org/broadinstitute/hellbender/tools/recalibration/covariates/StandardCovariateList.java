package org.broadinstitute.hellbender.tools.recalibration.covariates;

import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.tools.recalibration.ReadCovariates;
import org.broadinstitute.hellbender.tools.recalibration.RecalibrationArgumentCollection;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

/**
 * We always us the same covariates - this is the list.
 */
public final class StandardCovariateList implements Iterable<Covariate>{

    // enforce the order with RG first and QS next.
    private final List<Covariate> theList;

    private final Covariate readGroupCovariate;
    private final Covariate qualityScoreCovariate;
    private final Covariate contextCovariate;
    private final Covariate cycleCovariate;

    public StandardCovariateList(){
        this.readGroupCovariate = new ReadGroupCovariate();
        this.qualityScoreCovariate = new QualityScoreCovariate();
        this.contextCovariate = new ContextCovariate();
        this.cycleCovariate = new CycleCovariate();
       this.theList = Arrays.asList(readGroupCovariate, qualityScoreCovariate, contextCovariate, cycleCovariate);
    }

    public final int size(){
        return theList.size();
    }

    @Override
    public Iterator<Covariate> iterator() {
        return theList.iterator();
    }


    public Covariate getReadGroupCovariate() {
        return readGroupCovariate;
    }

    public Covariate getQualityScoreCovariate() {
        return qualityScoreCovariate;
    }

    public Covariate getContextCovariate() {
        return contextCovariate;
    }

    public Covariate getCycleCovariate() {
        return cycleCovariate;
    }

    /**
     * Return a human-readable string representing the used covariates
     *
     * @return a non-null comma-separated string
     */
    public String covariateNames() {
        final List<String> names = new ArrayList<>(this.size());
        for ( final Covariate cov : theList ) {
            names.add(cov.getClass().getSimpleName());
        }
        return String.join(",", names);
    }

    private void recordValueInStorage(Covariate cov, SAMRecord read, ReadCovariates resultsStorage) {
        int index = theList.indexOf(cov);
        resultsStorage.setCovariateIndex(index);
        cov.recordValues(read, resultsStorage);
    }

    public void recordAllValuesInStorage(SAMRecord read, ReadCovariates resultsStorage) {
        // Loop through the list of requested covariates and compute the values of each covariate for all positions in this read
        for(Covariate cov : this) {
            recordValueInStorage(cov, read, resultsStorage);
        }
    }

    /**
     * Get the covariate by the index.
     * XXX Exposing the index this is nasty, we need to eliminate this.
     */
    public Covariate get(int covIndex) {
        return theList.get(covIndex);
    }

    /**
     * Return the index at which the optional covariates start.
     * XXX This is a totally bogus. The list nature of this datastructure should not leak to clients but a lot of code depends on it.
     */
    public int getOptionalCovariatesStartIndex() {
        return Math.min(theList.indexOf(contextCovariate), theList.indexOf(cycleCovariate));
    }

    public StandardCovariateList initializeAll(RecalibrationArgumentCollection rac) {
        for (Covariate c : theList){
            c.initialize(rac);
        }
        return this;
    }
}
