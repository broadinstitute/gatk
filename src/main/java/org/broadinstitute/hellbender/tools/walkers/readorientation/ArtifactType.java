package org.broadinstitute.hellbender.tools.walkers.readorientation;

/**
 * Created by tsato on 3/28/18.
 */

// Artifact type is a misnomer: use read orientation
public enum ArtifactType {
    F1R2,
    F2R1;

    public static ArtifactType getOtherType(ArtifactType artifact){
        return artifact == F1R2 ? F2R1 : F1R2;
    }
}

