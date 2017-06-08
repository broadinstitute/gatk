package org.broadinstitute.hellbender.tools.archive;

import au.com.bytecode.opencsv.CSVWriter;
import org.apache.commons.math3.linear.DiagonalMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollection;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollectionUtils;
import org.broadinstitute.hellbender.utils.SparkToggleCommandLineProgram;
import org.broadinstitute.hellbender.utils.svd.SVD;
import org.broadinstitute.hellbender.utils.svd.SVDFactory;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;


@CommandLineProgramProperties(
        summary = "(EXPERIMENTAL) Performs a SVD",
        oneLineSummary = "(EXPERIMENTAL) Run a SVD on Spark and write the matrix to a text file (tsv).  The output file has no header and is just the raw numbers.  This tool is unsupported.",
        programGroup = CopyNumberProgramGroup.class
)
public final class DecomposeSingularValues extends SparkToggleCommandLineProgram {
    private static final long serialVersionUID = 1L;

    protected static final String INPUT_FILE_SHORT_NAME = "i";
    protected static final String INPUT_FILE_LONG_NAME = "inputFile";

    protected static final String OUTPUT_V_FILE_SHORT_NAME = "ov";
    protected static final String OUTPUT_V_FILE_LONG_NAME = "outputFileV";
    protected static final String OUTPUT_S_FILE_SHORT_NAME = "os";
    protected static final String OUTPUT_S_FILE_LONG_NAME = "outputFileS";
    protected static final String OUTPUT_U_FILE_SHORT_NAME = "ou";
    protected static final String OUTPUT_U_FILE_LONG_NAME = "outputFileU";

    @Argument(doc = "Input tsv file to be SVD'd.",
            fullName = INPUT_FILE_LONG_NAME,
            shortName = INPUT_FILE_SHORT_NAME,
            optional = false
    )
    protected File inputFile;

    @Argument(doc = "Output tsv file for the SVD U.  No headers, just tab-separated values",
            fullName = OUTPUT_U_FILE_LONG_NAME,
            shortName = OUTPUT_U_FILE_SHORT_NAME,
            optional = false
    )
    protected File outputFileU;

    @Argument(doc = "Output tsv file for the SVD S.  No headers, just tab-separated values",
            fullName = OUTPUT_S_FILE_LONG_NAME,
            shortName = OUTPUT_S_FILE_SHORT_NAME,
            optional = false
    )
    protected File outputFileS;

    @Argument(doc = "Output tsv file for the SVD U.  No headers, just tab-separated values",
            fullName = OUTPUT_V_FILE_LONG_NAME,
            shortName = OUTPUT_V_FILE_SHORT_NAME,
            optional = false
    )
    protected File outputFileV;


    @Override
    protected void runPipeline(final JavaSparkContext ctx) {
        try {
            final ReadCountCollection rcc = ReadCountCollectionUtils.parse(inputFile);
            final SVD svd = SVDFactory.createSVD(rcc.counts(), ctx);

            writeMatrix(svd.getV(), outputFileV);
            writeMatrix(svd.getU(), outputFileU);
            writeMatrix(new DiagonalMatrix(svd.getSingularValues()), outputFileS);
        } catch (final IOException ioe) {
            throw new UserException.CouldNotReadInputFile(inputFile, ioe.getMessage());
        }
    }

    private void writeMatrix(final RealMatrix m, final File outputFilename) throws IOException {
        final List<String []> textTable = new ArrayList<>();
        for (int i = 0; i < m.getRowDimension(); i ++){
            textTable.add(Arrays.stream(m.getRow(i)).mapToObj(Double::toString).toArray(String[]::new));
        }
        final FileWriter fw = new FileWriter(outputFilename);
        CSVWriter csvWriter = new CSVWriter(fw, '\t', CSVWriter.NO_QUOTE_CHARACTER);
        csvWriter.writeAll(textTable);
        csvWriter.flush();
        csvWriter.close();
    }
}
