package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.utils.samples.PedigreeValidationType;
import org.broadinstitute.hellbender.utils.samples.SampleDBBuilder;

import java.io.File;
import java.util.*;

public abstract class FounderAnnotation extends PedigreeAnnotation {
    final private Set<String> founderIds = new HashSet<>();

    /**
     * No-arg constructor is required for subclasses to be command line plugins.
     */
    protected FounderAnnotation() {}

    public FounderAnnotation(final File pedigreeFile){
        super.setPedigreeFile(pedigreeFile);
        addFoundersFromPedigreeFile(pedigreeFile);
    }

    public FounderAnnotation(final Set<String> founderIds) {
        super();
        setFounderIds(founderIds);
    }

    public void setFounderIds(Set<String> founderIds) {
        if (!founderIds.isEmpty()) {
            logger.warn("Replacing existing founder IDs");
            founderIds.clear();
        }
        founderIds.addAll(founderIds);
    }

    public Collection<String> getFounderIds() {
        return founderIds;
    }

    /**
     * Setter for pedigree file and founderIDs to be used by the GATKAnnotationPluginDescriptor to handle duplicated annotation
     * arguments between InbreedingCoeff and ExcessHet
     */
    public void setPedigreeFile(final File pedigreeFile) {
        super.setPedigreeFile(pedigreeFile);
        addFoundersFromPedigreeFile(pedigreeFile);
    }

    private void addFoundersFromPedigreeFile(final File pedigreeFile) {
        final SampleDBBuilder sampleDBBuilder = new SampleDBBuilder(PedigreeValidationType.STRICT);
        sampleDBBuilder.addSamplesFromPedigreeFiles(Collections.singletonList(pedigreeFile));
        final Set<String> founderIDsFromFile = sampleDBBuilder.getFinalSampleDB().getFounderIds();
        setFounderIds(founderIDsFromFile);
    }

    protected GenotypesContext getFounderGenotypes(VariantContext vc) {
        return (founderIds.isEmpty() ?
                vc.getGenotypes() :
                vc.getGenotypes(founderIds));
    }

}
