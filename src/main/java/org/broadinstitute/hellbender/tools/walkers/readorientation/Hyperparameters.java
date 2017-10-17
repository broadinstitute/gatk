package org.broadinstitute.hellbender.tools.walkers.readorientation;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;
import org.broadinstitute.hellbender.tools.walkers.readorientation.LearnHyperparametersEngine.ArtifactState;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Created by tsato on 9/5/17.
 */
public class Hyperparameters {
    private final String referenceContext;
    private final double[] pi;
    private final int numExamples;
    private final int numAltExamples;

    public Hyperparameters(final String referenceContext, final double[] pi, final int numExamples, final int numAltExamples) {
        this.referenceContext = referenceContext;
        this.pi = pi;
        this.numExamples = numExamples;
        this.numAltExamples = numAltExamples;
    }

    public double[] getPi() {
        return pi;
    }

    public double getPi(final LearnHyperparametersEngine.ArtifactState state) {
        return pi[state.ordinal()];
    }

    public String getReferenceContext() {
        return referenceContext;
    }

    public int getNumExamples() { return numExamples; }

    public int getNumAltExamples() { return numAltExamples; }

    /** Reading and writing the learned hyperparameters of the model **/

    /** Code for writing hyperparameters to a table **/
    // TODO: I can abstract this, pretty much copied the code from PileupSummary
    private static class HyperparameterTableWriter extends TableWriter<Hyperparameters> {
        private HyperparameterTableWriter(final File output) throws IOException {
            super(output, HyperparameterTableColumn.COLUMNS);
        }

        @Override
        protected void composeLine(final Hyperparameters hps, final DataLine dataLine) {
            // it'd be nice to set() less manually...
            // Note that allele fraction f is not allele-specific, thus the same f array will be printed
            // four times for each context
            dataLine.set(HyperparameterTableColumn.CONTEXT.toString(), hps.getReferenceContext())
                    .set(HyperparameterTableColumn.F1R2_A.toString(),
                            hps.getPi(ArtifactState.F1R2_A))
                    .set(HyperparameterTableColumn.F1R2_C.toString(),
                            hps.getPi(ArtifactState.F1R2_C))
                    .set(HyperparameterTableColumn.F1R2_G.toString(),
                            hps.getPi(LearnHyperparametersEngine.ArtifactState.F1R2_G))
                    .set(HyperparameterTableColumn.F1R2_T.toString(),
                            hps.getPi(ArtifactState.F1R2_T))
                    .set(HyperparameterTableColumn.F2R1_A.toString(),
                            hps.getPi(LearnHyperparametersEngine.ArtifactState.F2R1_A))
                    .set(HyperparameterTableColumn.F2R1_C.toString(),
                            hps.getPi(ArtifactState.F2R1_C))
                    .set(HyperparameterTableColumn.F2R1_G.toString(),
                            hps.getPi(ArtifactState.F2R1_G))
                    .set(HyperparameterTableColumn.F2R1_T.toString(),
                            hps.getPi(ArtifactState.F2R1_T))
                    .set(HyperparameterTableColumn.HOM_REF.toString(),
                            hps.getPi(ArtifactState.HOM_REF))
                    .set(HyperparameterTableColumn.GERMLINE_HET.toString(),
                            hps.getPi(ArtifactState.GERMLINE_HET))
                    .set(HyperparameterTableColumn.SOMATIC_HET.toString(),
                            hps.getPi(ArtifactState.SOMATIC_HET))
                    .set(HyperparameterTableColumn.HOM_VAR.toString(),
                            hps.getPi(LearnHyperparametersEngine.ArtifactState.HOM_VAR))
                    .set(HyperparameterTableColumn.N.toString(), hps.getNumExamples())
                    .set(HyperparameterTableColumn.N_ALT.toString(), hps.getNumAltExamples());
        }
    }

    /**
     *  Converts a double array to a comma-separated string without surrounding brackets
     *  Compare to Arrays.toString(), which includes brackets and therefore is inferior
     */
    public static String doubleArrayToString(final double[] xs){
        Utils.validateArg(xs.length > 0, "xs must not be an empty (uninitialized?) array");
        StringBuilder sb = new StringBuilder(String.valueOf(xs[0]));
        for (int i = 1; i < xs.length ; i++){
            sb.append("," + String.valueOf(xs[i]));
        }
        return sb.toString();
    }

    public static void writeHyperparameters(final List<Hyperparameters> hyperparameters, final File outputTable) {
        try (HyperparameterTableWriter writer = new HyperparameterTableWriter(outputTable)) {
            writer.writeAllRecords(hyperparameters);
        } catch (IOException e) {
            throw new UserException(String.format("Encountered an IO exception while reading from %s.", outputTable), e);
        }
    }

    /** Code for reading hyperparameters from a table **/
    public static List<Hyperparameters> readHyperparameters(final File table) {
        try (HyperParameterTableReader reader = new HyperParameterTableReader(table)) {
            return reader.toList();
        } catch (IOException e) {
            throw new UserException(String.format("Encountered an IO exception while reading from %s.", table), e);
        }
    }

    /** Read a row for specific hyperparameters from the table **/
    public static Hyperparameters readHyperparameters(final File table, final String referenceContext) {
        Utils.validateArg(referenceContext.length() == 3, "reference context must be 3 bases long");
        try (HyperParameterTableReader reader = new HyperParameterTableReader(table)) {
            // TODO: might be worth revisiting if this is the right approach
            for (Hyperparameters hyp : reader){
                if (hyp.getReferenceContext().equals(referenceContext)) {
                    return hyp;
                }
            }

            throw new UserException(String.format("Reference context %s does not exist in the hyperparameter table"));
        } catch (IOException e) {
            throw new UserException(String.format("Encountered an IO exception while reading from %s.", table), e);
        }
    }


    private static class HyperParameterTableReader extends TableReader<Hyperparameters> {
        private HyperParameterTableReader(final File table) throws IOException {
            super(table);
        }

        @Override
        protected Hyperparameters createRecord(final DataLine dataLine) {
            final String referenceContext = dataLine.get(HyperparameterTableColumn.CONTEXT);
            final double[] pi = new double[ArtifactState.values().length];
            for (HyperparameterTableColumn column : HyperparameterTableColumn.getArtifactStateColumns()){
                pi[column.getState().ordinal()] = Double.valueOf(dataLine.get(column));
            }

            final int numExamples = Integer.parseInt(dataLine.get(HyperparameterTableColumn.N));
            final int numAltExamples = Integer.parseInt(dataLine.get(HyperparameterTableColumn.N_ALT));
            return new Hyperparameters(referenceContext, pi, numExamples, numAltExamples);
        }
    }

    private enum HyperparameterTableColumn {
        CONTEXT("context"),
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
        HOM_VAR("hom_var", LearnHyperparametersEngine.ArtifactState.HOM_VAR),
        N("num_examples"),
        N_ALT("num_alt_examples");

        private String columnName;
        private LearnHyperparametersEngine.ArtifactState state;

        HyperparameterTableColumn(final String columnName) {
            this.columnName = columnName;
        }

        HyperparameterTableColumn(final String columnName, final ArtifactState state) {
            this.columnName = columnName;
            this.state = state;
        }

        @Override
        public String toString() {
            return columnName;
        }

        // If used on N or N_ALT, this method will return null
        public ArtifactState getState() { return state; }

        public static List<HyperparameterTableColumn> getArtifactStateColumns(){
            final List<HyperparameterTableColumn> nonArtifactColumns = Arrays.asList(CONTEXT, N, N_ALT);
            return Arrays.stream(HyperparameterTableColumn.values())
                    .filter(column -> ! nonArtifactColumns.contains(column))
                    .collect(Collectors.toList());
        }

        public static final TableColumnCollection COLUMNS = new TableColumnCollection((Object[]) values());
    }

}

