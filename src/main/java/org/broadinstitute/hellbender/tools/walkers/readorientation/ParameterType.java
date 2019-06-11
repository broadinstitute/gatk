package org.broadinstitute.hellbender.tools.walkers.readorientation;

public enum ParameterType {
    ARTIFACT_PRIOR("state_prior"), MEAN_ALT_F1R2_FRACTION("mean_f1r2");

    final String label;
    ParameterType(final String label){
        this.label = label;
    }

    public String getLabel(){
        return this.label;
    }
}
