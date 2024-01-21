package org.broadinstitute.hellbender.tools.walkers.featuremapping;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.util.SequenceUtil;
import ml.dmlc.xgboost4j.java.Booster;
import ml.dmlc.xgboost4j.java.DMatrix;
import ml.dmlc.xgboost4j.java.XGBoost;
import ml.dmlc.xgboost4j.java.XGBoostError;
import org.apache.commons.lang3.StringUtils;
import org.apache.hadoop.shaded.org.apache.http.util.Asserts;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.walkers.groundtruth.SeriesStats;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.json.JSONObject;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import static htsjdk.samtools.SAMUtils.MAX_PHRED_SCORE;

public class SvnQualReadProcessor {

    private static final Logger logger = LogManager.getLogger(SvnQualReadProcessor.class);
    public static final int JOINT_PROB_WARNING_TIMES = 5;
    public static final String CONF_KEY_QUALITY_INTERPOLATION_DATA = "quality_interpolation_data";

    final private ApplySNVQRArgumentCollection aqArgs;
    final private  List<SnvqrFeature> snvqrFeatures = new LinkedList<>();
    final private Booster  booster;
    final private SeriesStats bqStats = new SeriesStats();
    final private SeriesStats snvStats = new SeriesStats();
    final private Interpolator qualityInterpolator;

    final private SeriesStats zeroCountStats = new SeriesStats();

    private int jointProbWarnCount = 0;

    public SvnQualReadProcessor(ApplySNVQRArgumentCollection aqArgs, final JSONObject conf) {

        // save params
        this.aqArgs = aqArgs;

        try {
            // load model
            booster = !StringUtils.isEmpty(aqArgs.model) ? XGBoost.loadModel(aqArgs.model) : null;

            // load features
            if ( booster != null ) {
                for (final Object name : booster.getFeatureNames()) {
                    snvqrFeatures.add(SnvqrFeatureFactory.getFeature(name.toString(), conf));
                }
            } else {
                snvqrFeatures.add(SnvqrFeatureFactory.getFeature("X_SCORE", conf));
            }
        } catch(XGBoostError e){
            throw new GATKException("", e);
        }

        // initialize interpolator
        if ( conf != null && conf.has(CONF_KEY_QUALITY_INTERPOLATION_DATA) ) {
            qualityInterpolator = new Interpolator(conf.getJSONObject(CONF_KEY_QUALITY_INTERPOLATION_DATA));
        } else {
            qualityInterpolator = null;
        }
    }

    public void incorporateReadFeatures(GATKRead read, List<MappedFeature> features, SAMFileHeader header) {

        try {
            if ( booster != null ) {
                // create DMatrix from read and features
                DMatrix data[] = buildReadDMatrices(read, features, SequenceUtil.VALID_BASES_UPPER, header);

                float[][][] predicts = new float[data.length][][];
                for (int i = 0; i < data.length; i++) {

                    // predict
                    predicts[i] = booster.predict(data[i]);

                    // incorporate predictions back to read
                    incorporatePredicts(read, predicts[i], SequenceUtil.VALID_BASES_UPPER[i], true);
                }
            } else {

                // create simulated predicts from read and features
                float[][] predicts[] = buildReadNoModelPredicts(read, features, SequenceUtil.VALID_BASES_UPPER, header);

                for (int i = 0; i < predicts.length; i++) {

                    // incorporate predictions back to read
                    incorporatePredicts(read, predicts[i], SequenceUtil.VALID_BASES_UPPER[i], false);
                }

            }

            buildNewQualString(read);

            verifyCorrectness(read);

        } catch (XGBoostError e) {
            throw new GATKException("", e);
        } catch (IOException e) {
            throw new GATKException("", e);
        }
    }

    private void verifyCorrectness(GATKRead read) {

        // get snvqr strings
        byte[][] qX = new byte[SequenceUtil.VALID_BASES_UPPER.length][];
        for ( int i = 0 ; i < qX.length ; i++  ) {
            qX[i] = read.getAttributeAsString(attrNameForNonCalledBase(SequenceUtil.VALID_BASES_UPPER[i])).getBytes();
        }

        // check that all lengths are the same
        for ( int i = 1 ; i < qX.length ; i++ ) {
            if ( qX[0].length != qX[i].length ) {
                throw new GATKException(String.format("%s: qX[0] length is different than qX[%d]", read.getName(), i));
            }
        }

        // check that there is exactly one zero in each base position - or at least accumulate stats
        for ( int ofs = 0 ; ofs < qX[0].length ; ofs++ ) {
            int zeroCount = 0;
            for ( int i = 0 ; i < qX.length ; i++ ) {
                if ( (qX[i][ofs]  - 33) == 0 ) {
                    zeroCount++;
                }
            }
            zeroCountStats.add(zeroCount);
            zeroCountStats.aux(zeroCount, qX[0].length - ofs - 1);
        }

    }

    private DMatrix[] buildReadDMatrices(GATKRead read, List<MappedFeature> features, final byte[] bases, SAMFileHeader header) throws XGBoostError, IOException {

        // allocate data
        float[][] data = new float[bases.length][];
        for ( int i = 0 ; i < bases.length ; i++ ) {
            data[i] = new float[features.size() * snvqrFeatures.size()];
        }

        // create data, row by row, optimize for value which are not base dependent
        int row_col = 0;
        int feature_i = 0;
        for ( MappedFeature feature : features ) {
            int snvqr_feature_i = 0;
            for ( SnvqrFeature snvqrFeature : snvqrFeatures ) {
                for ( int i = 0 ; i < bases.length ; i++ ) {
                    if ( feature_i == 0 || !snvqrFeature.isMappedFeatureIndependent() ) {
                        if (i == 0 || !snvqrFeature.isNonCalledBasedIndependent()) {
                            data[i][row_col] = snvqrFeature.getValue(feature, bases[i], header);
                        } else {
                            data[i][row_col] = data[0][row_col];
                        }
                    } else {
                        data[i][row_col] = data[0][snvqr_feature_i];
                    }
                }
                snvqr_feature_i++;
                row_col++;
            }
            feature_i++;
        }
        final int nRows = features.size();
        final int nCols = snvqrFeatures.size();

        // log?
        if ( aqArgs.debugReadName.size() > 0 && aqArgs.debugReadName.contains(read.getName()) ) {
            logInputMatrices(aqArgs, read, bases, nRows, nCols, data);
        }

        // create matrices
        DMatrix[] matrices = new DMatrix[bases.length];
        for ( int i = 0 ; i < bases.length ; i++ ) {
            matrices[i] = new DMatrix(data[i], nRows, nCols, 0.0f);
        }

        return matrices;
    }

    private float[][][] buildReadNoModelPredicts(GATKRead read, List<MappedFeature> features, final byte[] bases, SAMFileHeader header) throws XGBoostError, IOException {

        // allocate data
        float[][][] data = new float[bases.length][][];
        for ( int i = 0 ; i < bases.length ; i++ ) {
            data[i] = new float[features.size() * snvqrFeatures.size()][];
            for ( int j = 0 ; j < data[i].length ; j++ ) {
                data[i][j] = new float[1];
            }
        }

        // create data, row by row, optimize for value which are not base dependent
        int row_col = 0;
        int feature_i = 0;
        for ( MappedFeature feature : features ) {
            int snvqr_feature_i = 0;
            for ( SnvqrFeature snvqrFeature : snvqrFeatures ) {
                for ( int i = 0 ; i < bases.length ; i++ ) {
                    if ( feature_i == 0 || !snvqrFeature.isMappedFeatureIndependent() ) {
                        if (i == 0 || !snvqrFeature.isNonCalledBasedIndependent()) {
                            data[i][row_col][0] = 10.0f * snvqrFeature.getValue(feature, bases[i], header);
                        } else {
                            data[i][row_col] = data[0][row_col];
                        }
                    } else {
                        data[i][row_col] = data[0][snvqr_feature_i];
                    }
                }
                snvqr_feature_i++;
                row_col++;
            }
            feature_i++;
        }
        final int nRows = features.size();
        final int nCols = snvqrFeatures.size();

        // log?
        if ( aqArgs.debugReadName.size() > 0 && aqArgs.debugReadName.contains(read.getName()) ) {
            logInputMatrices(aqArgs, read, bases, nRows, nCols, data);
        }

        return data;
    }

    private void incorporatePredicts(GATKRead read, float[][] predicts, byte nonCalledBase, boolean predictsAreErrorProb) throws IOException {

        // log?
        if ( aqArgs.debugReadName.size() > 0 && aqArgs.debugReadName.contains(read.getName()) ) {
            logOutputPredictions(aqArgs, read, predicts, nonCalledBase);
        }

        // build quality string
        StringBuilder qualString = new StringBuilder();
        for ( int row = 0 ; row < predicts.length ; row++ ) {
            float v = predicts[row][0];
            byte qual = predictsAreErrorProb ? QualityUtils.errorProbToQual(v) : (byte)v;
            if ( qualityInterpolator != null ) {
                qual = (byte)qualityInterpolator.value(qual);
            }
            char ch = SAMUtils.phredToFastq(Math.max(0, Math.min(qual, MAX_PHRED_SCORE)));
            qualString.append(ch);
        }

        // verify snvqr is of the same length as the original quality string
        Asserts.check(qualString.length() == read.getBaseQualitiesNoCopy().length,
                "generated SNVQR length %d is different than quality length %d", qualString.length(), read.getBaseQualitiesNoCopy().length);

        // create attribute, named qa, qc, etc
        read.setAttribute(attrNameForNonCalledBase(nonCalledBase), qualString.toString());
    }

    private String attrNameForNonCalledBase(byte nonCalledBase) {
        return ("Q" + (char)nonCalledBase).toLowerCase();
    }

    private void buildNewQualString(GATKRead read) {

        // get snvqr strings
        byte[][] qX = new byte[SequenceUtil.VALID_BASES_UPPER.length][];
        for ( int i = 0 ; i < qX.length ; i++  ) {
            qX[i] = read.getAttributeAsString(attrNameForNonCalledBase(SequenceUtil.VALID_BASES_UPPER[i])).getBytes();
        }

        // build new quality string (BQ)
        byte[] bq = new byte[qX[0].length];
        byte[] readBases = read.getBasesNoCopy();
        for ( int offset = 0 ; offset < bq.length ; offset++ ) {

            // accumulate for all non-called-bases
            double v = 0.0;
            for ( int i = 0 ; i < qX.length ; i++  ) {
                if ( SequenceUtil.VALID_BASES_UPPER[i] != readBases[offset] ) {

                    int phred = SAMUtils.fastqToPhred((char)qX[i][offset]);
                    double prob = QualityUtils.qualToErrorProb(phred);
                    v += prob;

                    snvStats.add(phred);
                }
            }

            // limiting prob to range - should be removed?
            if ( v < 0.0f || v > 1.0f ) {

                // only warn if significantly out of range
                if ( v > 1.0001 ) {

                    if (jointProbWarnCount < JOINT_PROB_WARNING_TIMES) {
                        logger.warn("joint probability seems to be out of range: " + v);
                    }
                    jointProbWarnCount++;
                    if (jointProbWarnCount == JOINT_PROB_WARNING_TIMES) {
                        logger.warn("joint probability warning disabled after " + JOINT_PROB_WARNING_TIMES + " times");
                    }
                }

                v = Math.min(1.0f, Math.max(0.0f, v));
            }

            byte phred = QualityUtils.errorProbToQual(v);
            if ( phred <= 0 ) {
                logger.warn(String.format("%s: offset %d, BQ phred ZERO!", read.getName(), offset));
            }
            char ch = SAMUtils.phredToFastq(phred);
            bq[offset] = (byte)ch;

            bqStats.add(phred);
        }

        // verify bq is of the same length as the original quality string
        Asserts.check(bq.length == read.getBaseQualitiesNoCopy().length,
                "generated BQ length %d is different than quality length %d", bq.length, read.getBaseQualitiesNoCopy().length);

        if ( !aqArgs.replaceQualityMode ) {
            read.setAttribute("BQ", new String(bq));
        } else{
            read.setAttribute("OQ", new String(read.getBaseQualitiesNoCopy()));
            read.setBaseQualities(bq);
        }

        // log qualities
        // log?
        if ( aqArgs.debugReadName.size() > 0 && aqArgs.debugReadName.contains(read.getName()) ) {
            logQualities(read);
        }

    }

    public void writeStats(String path) throws IOException {

        // if a folder name provided, make sure it ends with /
        if ( (new File(path)).isDirectory() && !path.endsWith("/") )
            path += "/";

        final String bqPath = path + "_bq.csv";
        logger.info("writing BQ stats into " + bqPath);
        bqStats.csvWrite(bqPath);

        final String snvPath = path + "_snv.csv";
        logger.info("writing SNV stats into " + snvPath);
        snvStats.csvWrite(snvPath);
    }

    public void logZeroCounts() {
        logger.info("zeroCount report (details are location offset from end of read):");
        zeroCountStats.getBins().forEach((k, v) -> {
            logger.info(String.format("%d: %d - %s", k.intValue(), v.get(),
                    (k.intValue() != 1) ? zeroCountStats.getAuxBins().get(k).toDigest() : ""));
        });
    }

    private void logInputMatrices(final ApplySNVQRArgumentCollection aqArgs, final GATKRead read, final byte[] bases,
                                  final int nRows, final int nCols, float[][][] data) throws IOException, XGBoostError {

        // reduce matrix dimensions
        float[][] data2 = new float[data.length][];
        for ( int i = 0 ; i < data2.length ; i++ ) {
            data2[i] = new float[data[i].length];
            for ( int j = 0 ; j < data[i].length ; j++ ) {
                data2[i][j] = data[i][j][0];
            }
        }

        logInputMatrices(aqArgs, read, bases, nRows, nCols, data2);

    }
    private void logInputMatrices(final ApplySNVQRArgumentCollection aqArgs, final GATKRead read, final byte[] bases,
                                  final int nRows, final int nCols, float[][] data) throws IOException, XGBoostError {

        for ( int i = 0 ; i < bases.length ; i++ ) {

            // create input csv
            final String path = aqArgs.debugReadFolder + "/" + read.getName() + "." + (char)bases[i] + ".in";
            final String path2 = aqArgs.debugReadFolder + "/" + read.getName() + "_CAT." + (char)bases[i] + ".in";
            PrintWriter pw = new PrintWriter(path);
            PrintWriter pw2 = new PrintWriter(path2);
            logger.info("logging into " + path);
            logger.info("logging into " + path2);
            if ( booster != null ) {
                pw.print(StringUtils.join(booster.getFeatureNames(), ",") + "\n");
                pw2.print(StringUtils.join(booster.getFeatureNames(), ",") + "\n");
            }
            for (int row = 0; row < nRows; row++) {
                List<String> rowData = new LinkedList<>();
                List<String> rowData2 = new LinkedList<>();
                for (int col = 0; col < nCols; col++) {
                    float v = data[i][row * nCols + col];
                    rowData.add(Float.toString(v));
                    rowData2.add(snvqrFeatures.get(col).getCategoryValue(v).toString());
                }
                pw.print(StringUtils.join(rowData, ",") + "\n");
                pw2.print(StringUtils.join(rowData2, ",") + "\n");
            }
            pw.close();
            pw2.close();
        }
    }

    private void logOutputPredictions(ApplySNVQRArgumentCollection aqArgs, GATKRead read, float[][] predicts, byte nonCalledBase) throws IOException {

        // create input csv
        final String path = aqArgs.debugReadFolder + "/" + read.getName() +  "." + (char)nonCalledBase + ".out";
        PrintWriter pw = new PrintWriter(path);
        logger.info("logging into " + path);
        List<String> header = new LinkedList<>();
        for ( int col = 0 ; col < predicts[0].length ; col++ ) {
            header.add("P" + col);
        }
        pw.print(StringUtils.join(header, ",") + "\n");
        for ( int row = 0 ; row < predicts.length ; row++ ) {
            List<String> rowData = new LinkedList<>();
            for ( int col = 0 ; col < predicts[0].length ; col++ ) {
                rowData.add(Float.toString(predicts[row][col]));
            }
            pw.print(StringUtils.join(rowData, ",") + "\n");
        }
        pw.close();;
    }

    private void logQualities(GATKRead read) {
        logger.info("Read: " + read.getName());
        for ( final String name : "Q,qa,qc,qg,qt,BQ,OQ".split(",") ) {
            byte[] bytes;
            int floor = 33;
            if ( name.equals("Q") ) {
                bytes = read.getBaseQualitiesNoCopy();
                floor = 0;
            } else if ( !read.hasAttribute(name) ) {
                continue;
            } else {
                bytes = read.getAttributeAsByteArray(name);
            }

            List<Integer> list = new LinkedList<>();
            for ( final byte b : bytes )
                list.add(b - floor);

            logger.info(String.format("%2s: %s", name, StringUtils.join(list, ",")));
        }
    }
}
