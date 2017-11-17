package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.utils.samples.PedigreeValidationType;
import org.broadinstitute.hellbender.utils.samples.SampleDBBuilder;

import java.io.File;
import java.util.*;

/**
 * A common interface for handling annotations that require pedigree file information either in the form of explicitly
 * selected founderIDs or in the form of an imported pedigreeFile.
 *
 * In order to use the the behavior, simply extend Pedigree annotation and access its constructors, then call
 * getFounderGenotypes() to extract only the genotypes corresponding to requested founder samples or that appear as founders
 * in the provided pedigree file. If no founderIDs or pedigreeFiles are present, then it defaults to returning all genotypes.
 */
public abstract class PedigreeAnnotation extends InfoFieldAnnotation {
    private Collection<String> founderIds;
    private File pedigreeFile = null;
    private boolean hasAddedPedigreeFounders = false;

    protected GenotypesContext getFounderGenotypes(VariantContext vc) {
        if ((pedigreeFile!= null) && (!hasAddedPedigreeFounders)) {
            founderIds.addAll(initializeSampleDB(pedigreeFile));
            hasAddedPedigreeFounders=true;
        }
        return (founderIds == null || founderIds.isEmpty()) ? vc.getGenotypes() : vc.getGenotypes(new HashSet<>(founderIds));
    }

    public PedigreeAnnotation(final Set<String> founderIds){
        //If available, get the founder IDs and cache them. the IC will only be computed on founders then.
        this.founderIds = founderIds == null? new ArrayList<>() : new ArrayList<>(founderIds);
    }

    public PedigreeAnnotation(final File pedigreeFile){
        //If available, get the founder IDs and cache them. the IC will only be computed on founders then.
        this.pedigreeFile = pedigreeFile;
        founderIds = initializeSampleDB(pedigreeFile);
        hasAddedPedigreeFounders = true;
    }

    /**
     * Entry-point function to initialize the samples database from input data
     */
    private Set<String> initializeSampleDB(File pedigreeFile) {
        final SampleDBBuilder sampleDBBuilder = new SampleDBBuilder(PedigreeValidationType.STRICT);
        sampleDBBuilder.addSamplesFromPedigreeFiles(Collections.singletonList(pedigreeFile));
        return sampleDBBuilder.getFinalSampleDB().getFounderIds();
    }

    /**
     * Setter for pedigree file and founderIDs to be used by the GATKAnnotationPluginDescriptor to handle duplicated annotaiton
     * arguments between InbreedingCoeff and ExcessHet
     */
    public void setPedigreeFile(File pedigreeFile) {
        this.pedigreeFile = pedigreeFile;
        hasAddedPedigreeFounders = false;
    }
    public void setFounderIds(List<String> founderIds) {
        this.founderIds = founderIds;
        hasAddedPedigreeFounders = false;
    }
}
