package org.broadinstitute.hellbender.tools.walkers.readorientation;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import htsjdk.variant.variantcontext.Allele;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;

/**
 * Created by tsato on 9/5/17.
 */
public class Hyperparameters {
    private final String referenceContext;
    private final double[] pi;
    private final double[] f;
    private final double[] theta;
    private final int numExamples;
    private final int numAltExamples;

    public Hyperparameters(final String referenceContext, final double[] pi, final double[] f, final double[] theta,
                           final int numExamples, final int numAltExamples) {
        this.referenceContext = referenceContext;
        this.pi = pi;
        this.f = f;
        this.theta = theta;
        this.numExamples = numExamples;
        this.numAltExamples = numAltExamples;
    }

    public double[] getPi() {
        return pi;
    }

    public double[] getF() {
        return f;
    }

    public double[] getTheta() {
        return theta;
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
                    .set(HyperparameterTableColumn.PI.toString(), doubleArrayToString(hps.getPi()))
                    .set(HyperparameterTableColumn.F.toString(), doubleArrayToString(hps.getF()))
                    .set(HyperparameterTableColumn.THETA.toString(), doubleArrayToString(hps.getTheta()))
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
            final double[] pi = Arrays.stream(dataLine.get(HyperparameterTableColumn.PI).split(","))
                    .mapToDouble(Double::parseDouble).toArray();
            final double[] f = Arrays.stream(dataLine.get(HyperparameterTableColumn.F).split(","))
                    .mapToDouble(Double::parseDouble).toArray();
            final double[] theta = Arrays.stream(dataLine.get(HyperparameterTableColumn.THETA).split(","))
                    .mapToDouble(Double::parseDouble).toArray();
            final int numExamples = Integer.parseInt(dataLine.get(HyperparameterTableColumn.N));
            final int numAltExamples = Integer.parseInt(dataLine.get(HyperparameterTableColumn.N_ALT));
            return new Hyperparameters(referenceContext, pi, f, theta, numExamples, numAltExamples);
        }
    }

    private enum HyperparameterTableColumn {
        CONTEXT("context"),
        PI("pi"),
        F("allele_fraction"),
        THETA("alt_f1r2_fraction"),
        N("num_examples"),
        N_ALT("num_alt_examples");

        private String columnName;

        HyperparameterTableColumn(final String columnName) {
            this.columnName = columnName;
        }

        @Override
        public String toString() {
            return columnName;
        }

        public static final TableColumnCollection COLUMNS = new TableColumnCollection((Object[]) values());
    }

}

