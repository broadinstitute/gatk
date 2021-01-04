package org.broadinstitute.hellbender.tools.sv;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.StructuralVariantType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.VariantContextGetters;

import java.util.ArrayList;
import java.util.List;

public class SVGenotypeEngineFromModel extends SVGenotypeEngine {

    private static final String COLUMN_SEPARATOR = "\t";
    private static final String FIRST_DIM_SEPARATOR = ";";
    private static final String SECOND_DIM_SEPARATOR = ",";

    public static boolean isDepthOnlyVariant(final VariantContext variant) {
        if (!variant.hasAttribute(GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE)) {
            throw new GATKException("Variant record is missing attribute: " + GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE);
        }
        final List<String> algorithms = variant.getAttributeAsStringList(GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE, "");
        return algorithms.size() == 1 && algorithms.contains(GATKSVVCFConstants.DEPTH_ALGORITHM);
    }

    public static List<VCFHeaderLine> getVcfHeaderMetadata() {
        final List<VCFHeaderLine> lines = new ArrayList<>();
        lines.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_QUALITY_KEY, true));
        lines.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_PL_KEY, true));
        lines.add(new VCFFormatHeaderLine(GATKSVVCFConstants.COPY_NUMBER_FIELD, 1, VCFHeaderLineType.Integer, "Copy number"));
        lines.add(new VCFInfoHeaderLine(GATKSVVCFConstants.PAIRED_END_PROB_FIELD, 1, VCFHeaderLineType.Float, "Paired-end read support probability"));
        lines.add(new VCFInfoHeaderLine(GATKSVVCFConstants.FIRST_SPLIT_READ_PROB_FIELD, 1, VCFHeaderLineType.Float, "First split read support probability"));
        lines.add(new VCFInfoHeaderLine(GATKSVVCFConstants.SECOND_SPLIT_READ_PROB_FIELD, 1, VCFHeaderLineType.Float, "Second read support probability"));
        lines.add(new VCFInfoHeaderLine(GATKSVVCFConstants.PAIRED_END_BACKGROUND_FIELD, 1, VCFHeaderLineType.Float, "Paired-end read median background rate"));
        lines.add(new VCFInfoHeaderLine(GATKSVVCFConstants.FIRST_SPLIT_READ_BACKGROUND_FIELD, 1, VCFHeaderLineType.Float, "First split read median background rate"));
        lines.add(new VCFInfoHeaderLine(GATKSVVCFConstants.SECOND_SPLIT_READ_BACKGROUND_FIELD, 1, VCFHeaderLineType.Float, "Second split read median background rate"));
        lines.add(new VCFInfoHeaderLine(GATKSVVCFConstants.PAIRED_END_MEDIAN_BIAS_FIELD, 1, VCFHeaderLineType.Float, "Paired-end read median bias"));
        lines.add(new VCFInfoHeaderLine(GATKSVVCFConstants.FIRST_SPLIT_READ_MEDIAN_BIAS_FIELD, 1, VCFHeaderLineType.Float, "First split read median bias"));
        lines.add(new VCFInfoHeaderLine(GATKSVVCFConstants.SECOND_SPLIT_READ_MEDIAN_BIAS_FIELD, 1, VCFHeaderLineType.Float, "Second split read median bias"));
        lines.add(new VCFInfoHeaderLine(GATKSVVCFConstants.MEDIAN_HARDY_WEINBERG_Q_FIELD, 1, VCFHeaderLineType.Float, "Median bi-allelic Hardy-Weinberg probability"));
        lines.add(new VCFInfoHeaderLine(GATKSVVCFConstants.MEDIAN_HARDY_WEINBERG_R_FIELD, 1, VCFHeaderLineType.Float, "Median bi-allelic Hardy-Weinberg probability"));

        lines.add(new VCFInfoHeaderLine(GATKSVVCFConstants.PAIRED_END_BACKGROUND_IQR_FIELD, 1, VCFHeaderLineType.Float, "Paired-end read median background rate IQR"));
        lines.add(new VCFInfoHeaderLine(GATKSVVCFConstants.FIRST_SPLIT_READ_BACKGROUND_IQR_FIELD, 1, VCFHeaderLineType.Float, "First split read median background rate IQR"));
        lines.add(new VCFInfoHeaderLine(GATKSVVCFConstants.SECOND_SPLIT_READ_BACKGROUND_IQR_FIELD, 1, VCFHeaderLineType.Float, "Second split read median background rate IQR"));
        lines.add(new VCFInfoHeaderLine(GATKSVVCFConstants.PAIRED_END_BIAS_IQR_FIELD, 1, VCFHeaderLineType.Float, "Paired-end read bias IQR"));
        lines.add(new VCFInfoHeaderLine(GATKSVVCFConstants.FIRST_SPLIT_READ_BIAS_IQR_FIELD, 1, VCFHeaderLineType.Float, "First split read bias IQR"));
        lines.add(new VCFInfoHeaderLine(GATKSVVCFConstants.SECOND_SPLIT_READ_BIAS_IQR_FIELD, 1, VCFHeaderLineType.Float, "Second split read bias IQR"));
        lines.add(new VCFInfoHeaderLine(GATKSVVCFConstants.HARDY_WEINBERG_Q_IQR_FIELD, 1, VCFHeaderLineType.Float, "Bi-allelic Hardy-Weinberg probability IQR"));
        lines.add(new VCFInfoHeaderLine(GATKSVVCFConstants.HARDY_WEINBERG_R_IQR_FIELD, 1, VCFHeaderLineType.Float, "Bi-allelic Hardy-Weinberg probability IQR"));

        return lines;
    }

    public VariantContext genotypeFromModel(final VariantContext variant, final String modelOutputLine, final List<String> modelSampleList) {
        Utils.nonNull(variant);
        Utils.nonNull(modelOutputLine);

        final VariantOutput modelOutput = parseVariantOutput(modelOutputLine, modelSampleList.size());
        if (!modelOutput.getId().equals(variant.getID())) {
            throw new UserException.BadInput("Model and VCF record IDs did not match: " + modelOutput.getId() + ", " + variant.getID());
        }

        final VariantContextBuilder builder = new VariantContextBuilder(variant);
        builder.attribute(GATKSVVCFConstants.PAIRED_END_PROB_FIELD, modelOutput.getP_m_pe());
        builder.attribute(GATKSVVCFConstants.FIRST_SPLIT_READ_PROB_FIELD, modelOutput.getP_m_sr1());
        builder.attribute(GATKSVVCFConstants.SECOND_SPLIT_READ_PROB_FIELD, modelOutput.getP_m_sr2());

        builder.attribute(GATKSVVCFConstants.PAIRED_END_BACKGROUND_FIELD, modelOutput.getEps_pe_median());
        builder.attribute(GATKSVVCFConstants.FIRST_SPLIT_READ_BACKGROUND_FIELD, modelOutput.getEps_sr1_median());
        builder.attribute(GATKSVVCFConstants.SECOND_SPLIT_READ_BACKGROUND_FIELD, modelOutput.getEps_sr2_median());
        builder.attribute(GATKSVVCFConstants.PAIRED_END_MEDIAN_BIAS_FIELD, modelOutput.getPhi_pe_median());
        builder.attribute(GATKSVVCFConstants.FIRST_SPLIT_READ_MEDIAN_BIAS_FIELD, modelOutput.getPhi_sr1_median());
        builder.attribute(GATKSVVCFConstants.SECOND_SPLIT_READ_MEDIAN_BIAS_FIELD, modelOutput.getPhi_sr2_median());
        builder.attribute(GATKSVVCFConstants.MEDIAN_HARDY_WEINBERG_Q_FIELD, modelOutput.getQ_hw_median());
        builder.attribute(GATKSVVCFConstants.MEDIAN_HARDY_WEINBERG_R_FIELD, modelOutput.getR_hw_median());

        builder.attribute(GATKSVVCFConstants.PAIRED_END_BACKGROUND_IQR_FIELD, modelOutput.getEps_pe_iqr());
        builder.attribute(GATKSVVCFConstants.FIRST_SPLIT_READ_BACKGROUND_IQR_FIELD, modelOutput.getEps_sr1_iqr());
        builder.attribute(GATKSVVCFConstants.SECOND_SPLIT_READ_BACKGROUND_IQR_FIELD, modelOutput.getEps_sr2_iqr());
        builder.attribute(GATKSVVCFConstants.PAIRED_END_BIAS_IQR_FIELD, modelOutput.getPhi_pe_iqr());
        builder.attribute(GATKSVVCFConstants.FIRST_SPLIT_READ_BIAS_IQR_FIELD, modelOutput.getPhi_sr1_iqr());
        builder.attribute(GATKSVVCFConstants.SECOND_SPLIT_READ_BIAS_IQR_FIELD, modelOutput.getPhi_sr2_iqr());
        builder.attribute(GATKSVVCFConstants.HARDY_WEINBERG_Q_IQR_FIELD, modelOutput.getQ_hw_iqr());
        builder.attribute(GATKSVVCFConstants.HARDY_WEINBERG_R_IQR_FIELD, modelOutput.getR_hw_iqr());

        final int numSamples = modelSampleList.size();
        final StructuralVariantType svType = variant.getStructuralVariantType();
        final List<Genotype> newGenotypes = new ArrayList<>(variant.getNSamples());
        for (int sampleIndex = 0; sampleIndex < numSamples; sampleIndex++) {
            final String sample = modelSampleList.get(sampleIndex);
            final Genotype genotype = variant.getGenotype(sample);
            if (!sample.equals(genotype.getSampleName())) {
                throw new UserException.BadInput("Model and VCF samples do not match");
            }
            final double[] genotypeProbs = modelOutput.getSampleFrequencies(sampleIndex);
            newGenotypes.add(genotypeFromGivenProbs(genotype, svType, genotypeProbs));
        }
        builder.genotypes(newGenotypes);
        final double log10ProbNoVariant = SVGenotypeEngine.calculateLog10PNoError(newGenotypes);
        builder.log10PError(log10ProbNoVariant);
        return builder.make();
    }

    public static VariantOutput parseVariantOutput(final String line, final int numSamples) {
        final String[] values = line.trim().split(COLUMN_SEPARATOR);
        final String id = values[0];
        final String[] freqStringArray = values[1].split(FIRST_DIM_SEPARATOR);
        if (freqStringArray.length != numSamples) {
            throw new UserException.BadInput("Genotype frequencies did not match sample list size");
        }
        final List<double[]> freqList = new ArrayList<>(freqStringArray.length);
        for (int i = 0; i < numSamples; i++) {
            final String[] sampleFreqStringArray = freqStringArray[i].split(SECOND_DIM_SEPARATOR);
            final double[] sampleFreq = new double[sampleFreqStringArray.length];
            for (int j = 0; j < sampleFreq.length; j++) {
                sampleFreq[j] = Double.parseDouble(sampleFreqStringArray[j]);
            }
            freqList.add(sampleFreq);
        }

        final double p_m_pe = Double.parseDouble(values[2]);
        final double p_m_sr1 = Double.parseDouble(values[3]);
        final double p_m_sr2 = Double.parseDouble(values[4]);

        final double eps_pe_median = Double.parseDouble(values[5]);
        final double eps_sr1_median = Double.parseDouble(values[6]);
        final double eps_sr2_median = Double.parseDouble(values[7]);
        final double phi_pe_median = Double.parseDouble(values[8]);
        final double phi_sr1_median = Double.parseDouble(values[9]);
        final double phi_sr2_median = Double.parseDouble(values[10]);
        final double q_hw_median = Double.parseDouble(values[11]);
        final double r_hw_median = Double.parseDouble(values[12]);

        final double eps_pe_iqr = Double.parseDouble(values[13]);
        final double eps_sr1_iqr = Double.parseDouble(values[14]);
        final double eps_sr2_iqr = Double.parseDouble(values[15]);
        final double phi_pe_iqr = Double.parseDouble(values[16]);
        final double phi_sr1_iqr = Double.parseDouble(values[17]);
        final double phi_sr2_iqr = Double.parseDouble(values[18]);
        final double q_hw_iqr = Double.parseDouble(values[19]);
        final double r_hw_iqr = Double.parseDouble(values[20]);

        return new VariantOutput(
                id,
                freqList,
                p_m_pe,
                p_m_sr1,
                p_m_sr2,
                eps_pe_median,
                eps_sr1_median,
                eps_sr2_median,
                phi_pe_median,
                phi_sr1_median,
                phi_sr2_median,
                q_hw_median,
                r_hw_median,
                eps_pe_iqr,
                eps_sr1_iqr,
                eps_sr2_iqr,
                phi_pe_iqr,
                phi_sr1_iqr,
                phi_sr2_iqr,
                q_hw_iqr,
                r_hw_iqr
        );
    }

    private static final class VariantOutput {
        private final String id;
        private final List<double[]> frequencies;
        private final double p_m_pe;
        private final double p_m_sr1;
        private final double p_m_sr2;
        private final double eps_pe_median;
        private final double eps_sr1_median;
        private final double eps_sr2_median;
        private final double phi_pe_median;
        private final double phi_sr1_median;
        private final double phi_sr2_median;
        private final double q_hw_median;
        private final double r_hw_median;
        private final double eps_pe_iqr;
        private final double eps_sr1_iqr;
        private final double eps_sr2_iqr;
        private final double phi_pe_iqr;
        private final double phi_sr1_iqr;
        private final double phi_sr2_iqr;
        private final double q_hw_iqr;
        private final double r_hw_iqr;

        public VariantOutput(final String id,
                             final List<double[]> frequencies,
                             final double p_m_pe,
                             final double p_m_sr1,
                             final double p_m_sr2,
                             final double eps_pe_median,
                             final double eps_sr1_median,
                             final double eps_sr2_median,
                             final double phi_pe_median,
                             final double phi_sr1_median,
                             final double phi_sr2_median,
                             final double q_hw_median,
                             final double r_hw_median,
                             final double eps_pe_iqr,
                             final double eps_sr1_iqr,
                             final double eps_sr2_iqr,
                             final double phi_pe_iqr,
                             final double phi_sr1_iqr,
                             final double phi_sr2_iqr,
                             final double q_hw_iqr,
                             final double r_hw_iqr) {
            this.id = id;
            this.frequencies = frequencies;
            this.p_m_pe = p_m_pe;
            this.p_m_sr1 = p_m_sr1;
            this.p_m_sr2 = p_m_sr2;
            this.eps_pe_median = eps_pe_median;
            this.eps_sr1_median = eps_sr1_median;
            this.eps_sr2_median = eps_sr2_median;
            this.phi_pe_median = phi_pe_median;
            this.phi_sr1_median = phi_sr1_median;
            this.phi_sr2_median = phi_sr2_median;
            this.q_hw_median = q_hw_median;
            this.r_hw_median = r_hw_median;
            this.eps_pe_iqr = eps_pe_iqr;
            this.eps_sr1_iqr = eps_sr1_iqr;
            this.eps_sr2_iqr = eps_sr2_iqr;
            this.phi_pe_iqr = phi_pe_iqr;
            this.phi_sr1_iqr = phi_sr1_iqr;
            this.phi_sr2_iqr = phi_sr2_iqr;
            this.q_hw_iqr = q_hw_iqr;
            this.r_hw_iqr = r_hw_iqr;
        }

        public int getNumGenotypes() {
            if (frequencies.isEmpty()) return 0;
            final int max = frequencies.stream().mapToInt(f -> f.length).max().getAsInt();
            final int min = frequencies.stream().mapToInt(f -> f.length).min().getAsInt();
            if (max != min) {
                throw new UserException.BadInput("Genotype frequency arrays not of uniform size");
            }
            return max;
        }

        public String getId() {
            return id;
        }

        public double[] getSampleFrequencies(final int sampleIndex) {
            Utils.validateArg(sampleIndex >= 0 && sampleIndex < frequencies.size(), "Invalid sample index: " + sampleIndex);
            return frequencies.get(sampleIndex);
        }

        public double getP_m_pe() {
            return p_m_pe;
        }

        public double getP_m_sr1() {
            return p_m_sr1;
        }

        public double getP_m_sr2() {
            return p_m_sr2;
        }

        public double getEps_pe_median() {
            return eps_pe_median;
        }

        public double getEps_sr1_median() {
            return eps_sr1_median;
        }

        public double getEps_sr2_median() {
            return eps_sr2_median;
        }

        public double getPhi_pe_median() {
            return phi_pe_median;
        }

        public double getPhi_sr1_median() {
            return phi_sr1_median;
        }

        public double getPhi_sr2_median() {
            return phi_sr2_median;
        }

        public double getQ_hw_median() {
            return q_hw_median;
        }

        public double getR_hw_median() {
            return r_hw_median;
        }

        public double getEps_pe_iqr() {
            return eps_pe_iqr;
        }

        public double getEps_sr1_iqr() {
            return eps_sr1_iqr;
        }

        public double getEps_sr2_iqr() {
            return eps_sr2_iqr;
        }

        public double getPhi_pe_iqr() {
            return phi_pe_iqr;
        }

        public double getPhi_sr1_iqr() {
            return phi_sr1_iqr;
        }

        public double getPhi_sr2_iqr() {
            return phi_sr2_iqr;
        }

        public double getQ_hw_iqr() {
            return q_hw_iqr;
        }

        public double getR_hw_iqr() {
            return r_hw_iqr;
        }
    }
}
