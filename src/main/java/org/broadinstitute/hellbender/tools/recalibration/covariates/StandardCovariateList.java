package org.broadinstitute.hellbender.tools.recalibration.covariates;

import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.tools.recalibration.ReadCovariates;
import org.broadinstitute.hellbender.tools.recalibration.RecalibrationArgumentCollection;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

/**
 * We always use the same covariates - this is the list.
 */
public final class StandardCovariateList implements Iterable<Covariate>{

    private final Covariate readGroupCovariate;
    private final Covariate qualityScoreCovariate;
    private final List<Covariate> additionalCovariates;
    private final List<Covariate> allCovariates;

    public StandardCovariateList(){
        this.readGroupCovariate = new ReadGroupCovariate();
        this.qualityScoreCovariate = new QualityScoreCovariate();
        ContextCovariate contextCovariate = new ContextCovariate();
        CycleCovariate cycleCovariate = new CycleCovariate();

        this.additionalCovariates = Arrays.asList(contextCovariate, cycleCovariate);
        this.allCovariates = Arrays.asList(readGroupCovariate, qualityScoreCovariate, contextCovariate, cycleCovariate);
    }

    public List<String> getStandardCovariateClassNames() {
        final List<String> names = new ArrayList<>(allCovariates.size());
        for ( final Covariate cov : allCovariates) {
            names.add(cov.getClass().getSimpleName());
        }
        return names;
    }

    public final int size(){
        return allCovariates.size();
    }

    @Override
    public Iterator<Covariate> iterator() {
        return allCovariates.iterator();
    }


    public Covariate getReadGroupCovariate() {
        return readGroupCovariate;
    }

    public Covariate getQualityScoreCovariate() {
        return qualityScoreCovariate;
    }

    public Iterable<Covariate> getAdditionalCovariates() {
        return additionalCovariates;
    }

    /**
     * Return a human-readable string representing the used covariates
     *
     * @return a non-null comma-separated string
     */
    public String covariateNames() {
        return String.join(",", getStandardCovariateClassNames());
    }

    /**
     * Get the covariate by the index.
     */
    public Covariate get(int covIndex) {
        return allCovariates.get(covIndex);
    }

    /**
     * Returns the index of the covariate by class name or -1 if not found.
     */
    public int indexByClass(Class<? extends Covariate> clazz){
        for(int i = 0; i < allCovariates.size(); i++){
            Covariate cov = allCovariates.get(i);
            if (cov.getClass().equals(clazz))  {
                return i;
            }
        }
        return -1;
    }

    private void recordValueInStorage(Covariate cov, SAMRecord read, ReadCovariates resultsStorage) {
        int index = indexByClass(cov.getClass());
        resultsStorage.setCovariateIndex(index);
        cov.recordValues(read, resultsStorage);
    }

    public void recordAllValuesInStorage(SAMRecord read, ReadCovariates resultsStorage) {
        // Loop through the list of requested covariates and compute the values of each covariate for all positions in this read
        for(Covariate cov : this) {
            recordValueInStorage(cov, read, resultsStorage);
        }
    }

    public StandardCovariateList initializeAll(RecalibrationArgumentCollection rac) {
        for (Covariate cov : allCovariates){
            cov.initialize(rac);
        }
        return this;
    }

    public Covariate getCovariateByParsedName(String covName) {
        for (Covariate cov : allCovariates){
            if (cov.parseNameForReport().equals(covName)) {
                return cov;
            }
        }
        throw new IllegalStateException("unknown covariate " + covName);
    }

}
