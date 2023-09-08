package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.utils.samples.PedigreeValidationType;
import org.broadinstitute.hellbender.utils.samples.SampleDBBuilder;
import org.broadinstitute.hellbender.utils.samples.Trio;

import javax.annotation.Nullable;
import java.util.*;

/**
 * A common interface for handling annotations that require pedigree file information either in the form of explicitly
 * selected founderIDs or in the form of an imported pedigreeFile.
 *
 * In order to use the behavior, simply extend Pedigree annotation and access its constructors, then call
 * getFounderGenotypes() to extract only the genotypes corresponding to requested founder samples or that appear as founders
 * in the provided pedigree file. If no founderIDs or pedigreeFiles are present, then it defaults to returning all genotypes.
 *
 * Alternatively, if a pedigree file has been supplied (not founderIDs) then extending classes can call getTrios() to
 * return a set of Trio objects corresponding to a parsing of pedigree file.
 */
public abstract class PedigreeAnnotation implements VariantAnnotation {
    private Collection<String> founderIds;
    private GATKPath pedigreeFile = null;
    private boolean hasAddedPedigreeFounders = false;
    protected transient final Logger logger = LogManager.getLogger(this.getClass());

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

    public PedigreeAnnotation(final GATKPath pedigreeFile){
        //If available, get the founder IDs and cache them. the IC will only be computed on founders then.
        this.pedigreeFile = pedigreeFile;
        initializeSampleDBAndSetFounders(pedigreeFile);
    }

    /**
     * Entry-point function to initialize the founders database from input data
     */
    private void initializeSampleDBAndSetFounders(GATKPath pedigreeFile) {
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
     * Setter for pedigree file and founderIDs to be used by the GATKAnnotationPluginDescriptor to handle duplicated annotation
     * arguments between InbreedingCoeff and ExcessHet
     */
    public void setPedigreeFile(GATKPath pedigreeFile) {
        this.pedigreeFile = pedigreeFile;
        hasAddedPedigreeFounders = false;
    }
    public void setFounderIds(List<String> founderIds) {
        this.founderIds = founderIds;
        hasAddedPedigreeFounders = false;
    }

    /**
     * Provide input arguments so as to warn the user if they are providing an incorrect subset of pedigree inputs.
     *
     * This is expected to be called immediately after calling setPedigreeFile() and setFounderIDs() during argument
     * propagation in the {@link org.broadinstitute.hellbender.cmdline.GATKPlugin.GATKAnnotationPluginDescriptor}.
     */
    public void validateArguments() {
        validateArguments(founderIds, pedigreeFile);
    }
    void validateArguments(Collection<String> founderIds, GATKPath pedigreeFile) {
        if ((founderIds == null || founderIds.isEmpty()) && pedigreeFile == null) {
            logger.warn(this.getClass().getSimpleName() + " annotation will not be calculated, no 'founder-id' or 'pedigree' arguments provided");
        }
    }

    /**
     * Warning generator for when a pedigree file is required and founderIDs cannot be used for this annotation.
     * @param founderIds
     * @param pedigreeFile
     */
    protected String validateArgumentsWhenPedigreeRequired(Collection<String> founderIds, GATKPath pedigreeFile) {
        if (pedigreeFile == null) {
            if ((founderIds != null && !founderIds.isEmpty())) {
                return "PossibleDenovo annotation will not be calculated, must provide a valid PED file (-ped). Founder-id arguments cannot be used for this annotation";
            } else {
                return "PossibleDenovo Annotation will not be calculated, must provide a valid PED file (-ped) from the command line.";
            }
        } else {
            if ((founderIds != null && !founderIds.isEmpty())) {
                return "PossibleDenovo annotation does not take founder-id arguments, trio information will be extracted only from the provided PED file";
            }
        }
        return null;
    }

    /**
     * This is a getter for the pedigree file, which could be null.
     *
     * @return The pedigree file
     */
    public @Nullable GATKPath getPedigreeFile() {
        return pedigreeFile;
    }

    /**
     * Helper function to check if the variant context has GQs for the trio
     * @param vc variant context
     * @param trio trio to check for GQs in the variant context
     */
    protected static boolean contextHasTrioGQs(final VariantContext vc, final Trio trio) {
        final String mom = trio.getMaternalID();
        final String dad = trio.getPaternalID();
        final String kid = trio.getChildID();

        return   (!mom.isEmpty() && vc.hasGenotype(mom) && vc.getGenotype(mom).hasGQ())
                && (!dad.isEmpty() && vc.hasGenotype(dad) && vc.getGenotype(dad).hasGQ())
                && (!kid.isEmpty() && vc.hasGenotype(kid) && vc.getGenotype(kid).hasGQ());
    }
}
