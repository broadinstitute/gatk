package org.broadinstitute.hellbender.tools.walkers.annotator;

import com.google.common.annotations.VisibleForTesting;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.samples.PedigreeValidationType;
import org.broadinstitute.hellbender.utils.samples.SampleDBBuilder;
import org.broadinstitute.hellbender.utils.samples.Trio;

import java.io.File;
import java.util.*;

/**
 * A common interface for handling annotations that require pedigree file information either in the form of explicitly
 * selected founderIDs or in the form of an imported pedigreeFile.
 *
 * In order to use the the behavior, simply extend Pedigree annotation and access its constructors, then call
 * getFounderGenotypes() to extract only the genotypes corresponding to requested founder samples or that appear as founders
 * in the provided pedigree file. If no founderIDs or pedigreeFiles are present, then it defaults to returning all genotypes.
 *
 * Alternatively, if a pedigree file has been supplied (not founderIDs) then extending classes can call getTrios() to
 * return a set of Trio objects corresponding to a parsing of pedigree file.
 */
public abstract class PedigreeAnnotation extends InfoFieldAnnotation {
    protected transient final Logger logger = LogManager.getLogger(this.getClass());

    protected File pedigreeFile;
    private final Set<Trio> trios = new HashSet<>();

    /**
     * No-arg constructor is required for command line plugins.
     */
    protected PedigreeAnnotation(){}

    public PedigreeAnnotation(final File pedigreeFile){
        setPedigreeFile(pedigreeFile);
    }

    @VisibleForTesting
    public PedigreeAnnotation(final Set<Trio> trios) {
        setTrios(trios);
    }

    public void setTrios(final Set<Trio> trios) {
        if (!getTrios().isEmpty()) {
            logger.warn("Replacing existing trios");
            this.trios.clear();
        }
        this.trios.addAll(trios);
    }

    /**
     * Gets the trios from the provided pedigree file or trio args
     */
    protected Set<Trio> getTrios() {
        return trios;
    }

    /**
     * Setter for pedigree file and founderIDs to be used by the GATKAnnotationPluginDescriptor to handle duplicated annotation
     * arguments between InbreedingCoeff and ExcessHet
     */
    public void setPedigreeFile(final File pedigreeFile) {
        this.pedigreeFile = pedigreeFile;
        final SampleDBBuilder sampleDBBuilder = new SampleDBBuilder(PedigreeValidationType.STRICT);
        sampleDBBuilder.addSamplesFromPedigreeFiles(Collections.singletonList(pedigreeFile));
        setTrios(sampleDBBuilder.getFinalSampleDB().getTrios());
    }

}
