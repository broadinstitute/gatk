package org.broadinstitute.hellbender.tools.walkers.readorientation;

import htsjdk.samtools.util.SequenceUtil;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;
import scala.collection.Seq;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

/**
 * Container for the artifact prior probabilities for the read orientation model
 */
public class ArtifactPrior {
    private final String referenceContext;
    private final double[] pi;
    private final int numExamples;
    private final int numAltExamples;

    public ArtifactPrior(final String referenceContext, final double[] pi, final int numExamples, final int numAltExamples) {
        this.referenceContext = referenceContext;
        this.pi = pi;
        this.numExamples = numExamples;
        this.numAltExamples = numAltExamples;
    }

    public double getPi(final ArtifactState state) {
        return pi[state.ordinal()];
    }

    public double[] getPi() {
        return pi;
    }

    public String getReferenceContext() {
        return referenceContext;
    }

    public String getRCContext() {
        return SequenceUtil.reverseComplement(referenceContext);
    }


    public ArtifactPrior getReverseComplement(){
        final double[] revCompPi = new double[F1R2FilterConstants.NUM_STATES];
        final String revCompRefContext = SequenceUtil.reverseComplement(referenceContext);
        for (final ArtifactState s : ArtifactState.values()){
            revCompPi[s.ordinal()] = pi[s.getRevCompState().ordinal()];
        }
        return new ArtifactPrior(revCompRefContext, revCompPi, numExamples, numAltExamples);
    }

    public int getNumExamples() { return numExamples; }

    public int getNumAltExamples() { return numAltExamples; }

    public static class ArtifactPriorTableWriter extends TableWriter<ArtifactPrior> {
        public ArtifactPriorTableWriter(final File output) throws IOException {
            super(output, ArtifactPriorTableColumn.COLUMNS);
        }

        @Override
        protected void composeLine(final ArtifactPrior artifactPrior, final DataLine dataLine) {
            dataLine.set(ArtifactPriorTableColumn.CONTEXT.toString(), artifactPrior.getReferenceContext())
                    .set(ArtifactPriorTableColumn.REV_COMP.toString(), artifactPrior.getRCContext())
                    .set(ArtifactPriorTableColumn.F1R2_A.toString(),
                            artifactPrior.getPi(ArtifactState.F1R2_A))
                    .set(ArtifactPriorTableColumn.F1R2_C.toString(),
                            artifactPrior.getPi(ArtifactState.F1R2_C))
                    .set(ArtifactPriorTableColumn.F1R2_G.toString(),
                            artifactPrior.getPi(ArtifactState.F1R2_G))
                    .set(ArtifactPriorTableColumn.F1R2_T.toString(),
                            artifactPrior.getPi(ArtifactState.F1R2_T))
                    .set(ArtifactPriorTableColumn.F2R1_A.toString(),
                            artifactPrior.getPi(ArtifactState.F2R1_A))
                    .set(ArtifactPriorTableColumn.F2R1_C.toString(),
                            artifactPrior.getPi(ArtifactState.F2R1_C))
                    .set(ArtifactPriorTableColumn.F2R1_G.toString(),
                            artifactPrior.getPi(ArtifactState.F2R1_G))
                    .set(ArtifactPriorTableColumn.F2R1_T.toString(),
                            artifactPrior.getPi(ArtifactState.F2R1_T))
                    .set(ArtifactPriorTableColumn.HOM_REF.toString(),
                            artifactPrior.getPi(ArtifactState.HOM_REF))
                    .set(ArtifactPriorTableColumn.GERMLINE_HET.toString(),
                            artifactPrior.getPi(ArtifactState.GERMLINE_HET))
                    .set(ArtifactPriorTableColumn.SOMATIC_HET.toString(),
                            artifactPrior.getPi(ArtifactState.SOMATIC_HET))
                    .set(ArtifactPriorTableColumn.HOM_VAR.toString(),
                            artifactPrior.getPi(ArtifactState.HOM_VAR))
                    .set(ArtifactPriorTableColumn.N.toString(), artifactPrior.getNumExamples())
                    .set(ArtifactPriorTableColumn.N_ALT.toString(), artifactPrior.getNumAltExamples());
        }
    }

    public static class ArtifactPriorTableReader extends TableReader<ArtifactPrior> {
        public ArtifactPriorTableReader(final File table) throws IOException {
            super(table);
        }

        @Override
        protected ArtifactPrior createRecord(final DataLine dataLine) {
            final String referenceContext = dataLine.get(ArtifactPriorTableColumn.CONTEXT);
            final double[] pi = new double[ArtifactState.values().length];
            for (ArtifactPriorTableColumn column : ArtifactPriorTableColumn.getArtifactStateColumns()){
                pi[column.getState().ordinal()] = Double.valueOf(dataLine.get(column));
            }

            final int numExamples = Integer.parseInt(dataLine.get(ArtifactPriorTableColumn.N));
            final int numAltExamples = Integer.parseInt(dataLine.get(ArtifactPriorTableColumn.N_ALT));
            return new ArtifactPrior(referenceContext, pi, numExamples, numAltExamples);
        }
    }

    private enum ArtifactPriorTableColumn {
        CONTEXT("context"),
        REV_COMP("rev_comp"),
        F1R2_A("f1r2_a", ArtifactState.F1R2_A),
        F1R2_C("f1r2_c", ArtifactState.F1R2_C),
        F1R2_G("f1r2_g", ArtifactState.F1R2_G),
        F1R2_T("f1r2_t", ArtifactState.F1R2_T),
        F2R1_A("f2r1_a", ArtifactState.F2R1_A),
        F2R1_C("f2r1_c", ArtifactState.F2R1_C),
        F2R1_G("f2r1_g", ArtifactState.F2R1_G),
        F2R1_T("f2r1_t", ArtifactState.F2R1_T),
        HOM_REF("hom_ref", ArtifactState.HOM_REF),
        GERMLINE_HET("germline_het", ArtifactState.GERMLINE_HET),
        SOMATIC_HET("somatic_het", ArtifactState.SOMATIC_HET),
        HOM_VAR("hom_var", ArtifactState.HOM_VAR),
        N("num_examples"),
        N_ALT("num_alt_examples");

        private String columnName;
        private ArtifactState state;

        ArtifactPriorTableColumn(final String columnName) {
            this.columnName = columnName;
        }

        ArtifactPriorTableColumn(final String columnName, final ArtifactState state) {
            this.columnName = columnName;
            this.state = state;
        }

        @Override
        public String toString() {
            return columnName;
        }

        // If used on N or N_ALT, this method will return null
        public ArtifactState getState() { return state; }

        public static List<ArtifactPriorTableColumn> getArtifactStateColumns(){
            final List<ArtifactPriorTableColumn> nonArtifactColumns = Arrays.asList(CONTEXT, REV_COMP, N, N_ALT);
            return Arrays.stream(ArtifactPriorTableColumn.values())
                    .filter(column -> ! nonArtifactColumns.contains(column))
                    .collect(Collectors.toList());
        }

        public static final TableColumnCollection COLUMNS = new TableColumnCollection((Object[]) values());
    }
}

