package org.broadinstitute.hellbender.tools.walkers.readorientation;

import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableUtils;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.IOException;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

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

    public static class PosteriorAltF1R2TableWriter extends TableWriter<PosteriorAltF1R2> {
        public PosteriorAltF1R2TableWriter(final Path output, final String sample) throws IOException {
            super(output, PosteriorAltF1R2TableColumn.COLUMNS);
            writeMetadata(TableUtils.SAMPLE_METADATA_TAG, sample);
        }

        @Override
        protected void composeLine(final PosteriorAltF1R2 posteriorAltF1R2, final DataLine dataLine) {
            dataLine.set(PosteriorAltF1R2TableColumn.CONTEXT.toString(), posteriorAltF1R2.getRefContext())
                    .set(PosteriorAltF1R2TableColumn.ARTIFACT_STATE.toString(), getRCContext())
                    .set(PosteriorAltF1R2TableColumn.ALPHA.toString(),
                            getPi(ArtifactState.F1R2_A))
                    .set(PosteriorAltF1R2TableColumn.BETA.toString(),
                            getPi(ArtifactState.F1R2_C))
        }
    }

    private enum PosteriorAltF1R2TableColumn {
        CONTEXT("context"),
        ARTIFACT_STATE("artifact_state"),
        ALPHA("alpha"),
        BETA("beta"),
        MEAN("mean"),
        VARIANCE("variance");

        private String columnName;
        private ArtifactState state;

        PosteriorAltF1R2TableColumn(final String columnName) {
            this.columnName = columnName;
        }

        @Override
        public String toString() {
            return columnName;
        }
        
        public static final TableColumnCollection COLUMNS = new TableColumnCollection((Object[]) values());
    }
}


