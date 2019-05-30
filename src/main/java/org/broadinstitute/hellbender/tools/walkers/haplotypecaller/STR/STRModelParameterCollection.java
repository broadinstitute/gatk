package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.STR;

import org.broadinstitute.hellbender.utils.exceptions.UserException;

import java.io.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Collection of STR parameter sets per covariate combination.
 */
public final class STRModelParameterCollection {

    public static final STRModelParameterCollection NULL = new STRModelParameterCollection(new STRModelCovariate.List(Collections.emptyList()), Collections.emptyList());

    STRModelParameterCollection(final STRModelCovariate.List covariates, final List<STRModelParameters> elements) {
        this.covariates = covariates;
        this.elements = elements;
    }

    /**
     * Writes the collection out into a text file.
     * @param file the output file.
     * @throws IOException if there is some problem when writting the collection into the file provided.
     */
    public void write(final File file) throws IOException {
        if (file == null) {
            throw new IllegalArgumentException("the input file cannot be null");
        }
        final PrintWriter writer = new PrintWriter(new FileWriter(file));
        writer.println("# STR Model Parameters ");
        for (final STRModelCovariate covariate : covariates) {
            writer.println("Covariate " + covariate.toString());
        }
        final List<String> columnNames = new ArrayList<>();
        columnNames.add("Set");
        columnNames.addAll(covariates.stream().map(Object::getClass).map(Class::getSimpleName).collect(Collectors.toList()));
        columnNames.add("Size");
        columnNames.add("Pi");
        columnNames.add("Lambda");
        writer.println(columnNames.stream().collect(Collectors.joining("\t", "#", "")));
        int index = 0;
        for (final STRModelParameters params : elements) {
            final List<String> values = new ArrayList<>();
            final int setIndex = index;
            values.add("Set");
            values.add("" + index);
            values.addAll(covariates.stream().map(cov -> "" + covariates.valueAt(setIndex, cov)).collect(Collectors.toList()));
            values.add("" + params.size);
            values.add("" + params.pi);
            values.add("" + params.lambda);
            writer.println(values.stream().collect(Collectors.joining("\t")));
            index++;
        }
        writer.close();
    }

    public static STRModelParameterCollection parse(final File file)
        throws IOException {
        final BufferedReader reader = new BufferedReader(new FileReader(file));
        final List<STRModelCovariate> covariates = new ArrayList<>();
        final List<STRModelParameters> sets = new ArrayList<>();
        String line;
        while ((line = reader.readLine()) != null)  {
            if (line.matches("^\\s*$") || line.matches("^#.*$")) { // skip comments and empty lines.
               // nothing.
            } else if (line.matches("^Covariate(s)?.*$")) {
                final String[] covariateSpecs = line.split("\\s+");
                for (int i = 1; i < covariateSpecs.length; i++) {
                    covariates.add(STRModelCovariate.parse(covariateSpecs[i]));
                }
            } else if (line.matches("^Set.*$")) {
                sets.add(STRModelParameters.parse(line));
            } else {
                throw new UserException.BadInput("invalid line: " + line);
            }
        }
        final STRModelCovariate.List covariateList = new STRModelCovariate.List(covariates);
        final int expectedSetCount = covariateList.combinationsCount();
        if (sets.size() != expectedSetCount) {
            throw new IllegalArgumentException("covariates dimension and number of set mismatch: " + sets.size() + " != " + expectedSetCount);
        }
        return new STRModelParameterCollection(covariateList, sets);
    }

    private STRModelCovariate.List covariates;

    public List<STRModelParameters> elements;

    public STRModelParameters fromContext(final STRContext context) {
        final int index = covariates.indexFor(context);
        return elements.size() > index ? elements.get(index) : null;
    }

}
