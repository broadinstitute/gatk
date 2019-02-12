package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
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
    private Collection<String> founderIds;
    private File pedigreeFile = null;
    private boolean hasAddedPedigreeFounders = false;

    protected GenotypesContext getFounderGenotypes(VariantContext vc) {
        if ((pedigreeFile!= null) && (!hasAddedPedigreeFounders)) {
            initializeSampleDBAndSetFounders(pedigreeFile);
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
        initializeSampleDBAndSetFounders(pedigreeFile);
    }

    /**
     * Entry-point function to initialize the founders database from input data
     */
    private void initializeSampleDBAndSetFounders(File pedigreeFile) {
        final SampleDBBuilder sampleDBBuilder = new SampleDBBuilder(PedigreeValidationType.STRICT);
        sampleDBBuilder.addSamplesFromPedigreeFiles(Collections.singletonList(pedigreeFile));

        Set<String> founderIdsToAdd = sampleDBBuilder.getFinalSampleDB().getFounderIds();
        if (this.founderIds == null || this.founderIds.isEmpty()) {
            this.founderIds = founderIdsToAdd;
        } else {
            this.founderIds.addAll(founderIdsToAdd);
        }
        hasAddedPedigreeFounders = true;
    }

    /**
     * Computes the trios from the provided pedigree file
     */
    protected Set<Trio> getTrios() {
        if (pedigreeFile!= null) {
            final SampleDBBuilder sampleDBBuilder = new SampleDBBuilder(PedigreeValidationType.STRICT);
            sampleDBBuilder.addSamplesFromPedigreeFiles(Collections.singletonList(pedigreeFile));
            return sampleDBBuilder.getFinalSampleDB().getTrios();
        }
        return Collections.emptySet();
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
