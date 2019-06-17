package org.broadinstitute.hellbender.tools.walkers.readorientation;

import org.broadinstitute.hellbender.exceptions.UserException;

import java.util.Arrays;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

public enum ParameterType {
    ARTIFACT_PRIOR("state_prior"), MEAN_ALT_F1R2_FRACTION("mean_f1r2");

    final String label;
    ParameterType(final String label){
        this.label = label;
    }

    public String getLabel(){
        return this.label;
    }

    public static List<String> getParameterTypes(){
        return Arrays.stream(ParameterType.values()).map(t -> t.getLabel()).collect(Collectors.toList());
    }

    public static ParameterType getParameterTypeFromLabel(final String query){
        final Optional<ParameterType> type = Arrays.stream(ParameterType.values()).filter(pt -> pt.getLabel().equals(query)).findAny();
        if (! type.isPresent()){
            throw new UserException("The label " + query + " is not a valid label for the enum ParameterType");
        }

        return type.get();
    }
}
