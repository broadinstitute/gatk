package org.broadinstitute.hellbender.tools.walkers.readorientation;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Nucleotide;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

/**
 * This enum encapsulates the domain of the discrete latent random variable z
 */
public enum ArtifactState {
    /* Artifact States - provide the reverse complement of each state in the constructor
     * The constructor must take string value of the rc states because it doesn't allow forward referencing
     * of an ArtifactState
     */
    F1R2_A("F2R1_T", Nucleotide.A), F1R2_C("F2R1_G", Nucleotide.C), F1R2_G("F2R1_C", Nucleotide.G), F1R2_T("F2R1_A", Nucleotide.T),
    F2R1_A("F1R2_T", Nucleotide.A), F2R1_C("F1R2_G", Nucleotide.C), F2R1_G("F1R2_C", Nucleotide.G), F2R1_T("F1R2_A", Nucleotide.T),

    // Non Artifact States
    HOM_REF("HOM_REF", Nucleotide.INVALID), GERMLINE_HET("GERMLINE_HET", Nucleotide.INVALID),
    SOMATIC_HET("SOMATIC_HET", Nucleotide.INVALID), HOM_VAR("HOM_VAR", Nucleotide.INVALID);

    final private String reverseComplementState;
    final private Nucleotide altAllele;

    ArtifactState(final String reverseComplementState, final Nucleotide altAllele){
        this.reverseComplementState = reverseComplementState;
        this.altAllele = altAllele;
    }

    public static List<ArtifactState> getStates(){
        return Arrays.stream(ArtifactState.values()).collect(Collectors.toList());
    }

    static List<ArtifactState> getF1R2ArtifactStates(){
        return Arrays.asList(F1R2_A, F1R2_C, F1R2_G, F1R2_T);
    }

    static List<ArtifactState> getF2R1ArtifactStates(){
        return Arrays.asList(F2R1_A, F2R1_C, F2R1_G, F2R1_T);
    }

    public static List<ArtifactState> getNonArtifactStates(){
        return Arrays.asList(HOM_REF, GERMLINE_HET, SOMATIC_HET, HOM_VAR);
    }

    // For a given reference base, return the ref to ref artifact states (e.g. AGT -> G), which we want to ignore
    public static List<ArtifactState> getRefToRefArtifacts(final Nucleotide refAllele){
        switch (refAllele){
            case A : return Arrays.asList( F1R2_A, F2R1_A );
            case C : return Arrays.asList( F1R2_C, F2R1_C );
            case G : return Arrays.asList( F1R2_G, F2R1_G );
            case T : return Arrays.asList( F1R2_T, F2R1_T );
            default: throw new UserException(String.format("Invalid nucleotide given: %s", refAllele));
        }
    }

    public Nucleotide getAltAlleleOfArtifact(){
        return altAllele;
    }

    public ArtifactState getRevCompState(){
        return ArtifactState.valueOf(reverseComplementState);
    }

    static List<ArtifactState> artifactStates = Arrays.asList(F1R2_A, F1R2_C, F1R2_G, F1R2_T, F2R1_A, F2R1_C, F2R1_G, F2R1_T);

    public static ArtifactState getF1R2StateForAlt(final Nucleotide altAllele) {
        switch (altAllele) {
            case A:
                return ArtifactState.F1R2_A;
            case C:
                return ArtifactState.F1R2_C;
            case G:
                return ArtifactState.F1R2_G;
            case T:
                return ArtifactState.F1R2_T;
            default:
                throw new UserException(String.format("Alt allele must be in {A, C, G, T} but got %s", altAllele));
        }
    }

    public static ArtifactState getF2R1StateForAlt(final Nucleotide altAllele) {
        switch (altAllele) {
            case A:
                return ArtifactState.F2R1_A;
            case C:
                return ArtifactState.F2R1_C;
            case G:
                return ArtifactState.F2R1_G;
            case T:
                return ArtifactState.F2R1_T;
            default:
                throw new UserException(String.format("Alt allele must be in {A, C, G, T} but got %s", altAllele));
        }
    }

    /**
     *
     * @param probability
     * @param state
     */
    public static void setStatePrior(final double[] prior, final double probability, final ArtifactState state){
        prior[state.ordinal()] = probability;
    }

}
