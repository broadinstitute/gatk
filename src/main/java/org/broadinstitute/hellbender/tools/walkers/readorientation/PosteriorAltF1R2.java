package org.broadinstitute.hellbender.tools.walkers.readorientation;

import java.util.HashMap;
import java.util.Map;

public class PosteriorAltF1R2 {
    String refContext;
    Map<ArtifactState, BetaDistributionShape> state2Posterior;

    public PosteriorAltF1R2(final String refContext){
        this.refContext = refContext;
        this.state2Posterior = new HashMap<>();
    }

    public void add(final ArtifactState artifactState, final BetaDistributionShape posterior){
        state2Posterior.put(artifactState, posterior);
    }

    public BetaDistributionShape get(final ArtifactState artifactState){
        return state2Posterior.get(artifactState);
    }

    public String getRefContext() {
        return refContext;
    }
}


