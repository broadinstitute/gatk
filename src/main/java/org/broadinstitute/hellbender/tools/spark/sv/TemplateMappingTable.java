package org.broadinstitute.hellbender.tools.spark.sv;

import org.apache.commons.collections4.CollectionUtils;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVHaplotype;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.Serializable;
import java.util.*;
import java.util.stream.IntStream;

/**
 * Created by valentin on 5/18/17.
 */
public final class TemplateMappingTable implements Serializable {

    private static final long serialVersionUID = 1L;

    private final TemplateMapping[][] mappingInfo;

    private final double[][] bestMappingScorePerFragment;
    private final double[][] bestMappingFragmentMQ;
    private final double[][] worstMappingScorePerFragment;
    private final double[][] worstMappingScoresFragmentMQ;
    private boolean needsToCalculateBestMappingScores = true;

    private final double[][] values;

    private final List<Template> templates;

    private final Map<String, Integer> templateIndexByName;

    private final List<SVHaplotype> haplotypes;

    public TemplateMappingTable(final Iterable<Template> templates, final Iterable<SVHaplotype> haplotypes)
    {
        this.templates = CollectionUtils.collect(templates, t -> t, new ArrayList<>(1000));
        this.haplotypes = CollectionUtils.collect(haplotypes, t -> t, new ArrayList<>());
        values = new double[this.haplotypes.size()][this.templates.size()];
        mappingInfo = new TemplateMapping[this.haplotypes.size()][this.templates.size()];
        bestMappingScorePerFragment = new double[this.templates.size()][2];
        bestMappingFragmentMQ = new double[this.templates.size()][2];
        worstMappingScorePerFragment = new double[this.templates.size()][2];
        worstMappingScoresFragmentMQ = new double[this.templates.size()][2];
        this.templateIndexByName = composeTemplateIndex(this.templates);
    }

    private Map<String,Integer> composeTemplateIndex(final List<Template> templates) {
        final Map<String, Integer> result = new HashMap<>(templates.size());
        for (int i = 0; i < templates.size(); i++) {
            result.put(templates.get(i).name(), i);
        }
        return result;
    }

    public double get(final SVHaplotype allele, final Template template) {
        return get(indexOf(allele), indexOf(template));
    }

    public int indexOf(final Template template) {
        Utils.nonNull(template);
        return templateIndexByName.getOrDefault(template.name(), -1);
    }

    public int indexOf(final SVHaplotype allele) {
        Utils.nonNull(allele);
        int i = 0;
        while (i < haplotypes.size()) {
            if (Arrays.equals(haplotypes.get(i).getBases(), allele.getBases())) {
                return i;
            }
            i++;
        }
        return -1;
    }

    public double get(final int alleleIndex, final int templateIndex) {
        Utils.validIndex(alleleIndex, haplotypes.size());
        Utils.validIndex(templateIndex, templates.size());
        return values[alleleIndex][templateIndex];
    }

    public void set(final int alleleIndex, final int templateIndex, final double value) {
        Utils.validIndex(alleleIndex, haplotypes.size());
        Utils.validIndex(templateIndex, templates.size());
        values[alleleIndex][templateIndex] = value;
    }

    public TemplateMapping getMappingInfo(final int alleleIndex, final int templateIndex) {
        Utils.validIndex(alleleIndex, haplotypes.size());
        Utils.validIndex(templateIndex, templates.size());
        return Objects.requireNonNull(mappingInfo[alleleIndex][templateIndex],  "requested null mapping information for template " + templates.get(templateIndex) + " on haplotype " + haplotypes.get(alleleIndex));
    }

    public void setMapping(final int alleleIndex, final int templateIndex, final TemplateMapping mappingInformation) {
        Utils.validIndex(alleleIndex, haplotypes.size());
        Utils.validIndex(templateIndex, templates.size());
        Objects.requireNonNull(mappingInfo[alleleIndex][templateIndex] = mappingInformation, "requested null mapping information for template " + templates.get(templateIndex) + " on haplotype " + haplotypes.get(alleleIndex));
    }

    public int numberOfTemplates() {
        return templates.size();
    }

    public int numberOfHaplotypes() {
        return haplotypes.size();
    }

    public List<Template> templates() {
        return templates;
    }

    public List<SVHaplotype> haplotypes() {
        return haplotypes;
    }


    public int[] informativeTemplateIndexes() {
        return IntStream.range(0, templates.size())
                .filter(i -> {
                    final double firstLikelihood = values[0][i];
                    if (Double.isNaN(firstLikelihood)) {
                        return false;
                    } else {
                        boolean foundDifference = false;
                        for (int j = 1; j < values.length; j++) {
                            if (Double.isNaN(values[j][i])) {
                                return false;
                            } else if (values[j][i] != firstLikelihood) {
                                foundDifference = true;
                            }
                        }
                        return foundDifference;
                    }
                }).toArray();
    }

    @SuppressWarnings("unused") // may be useful in the future.
    public void dropUninformativeTemplates() {
        final int[] informativeIndexes = informativeTemplateIndexes();
        if (informativeIndexes.length == 0) {
            templates.clear();
            templateIndexByName.clear();
            for (int j = 0; j < values.length; j++) {
                values[j] = new double[0];
            }
        } else if (informativeIndexes.length != templates.size()) {
            final List<Template> newTemplates = new ArrayList<>(informativeIndexes.length);
            templateIndexByName.clear();
            for (final int informativeIndexe : informativeIndexes) {
                newTemplates.add(templates.get(informativeIndexe));
                templateIndexByName.put(newTemplates.get(newTemplates.size() - 1).name(), newTemplates.size() - 1);
            }
            templates.clear();
            templates.addAll(newTemplates);
            for (int j = 0; j < haplotypes.size(); j++) {
                final double[] newValues = new double[informativeIndexes.length];
                final TemplateMapping[] newMappingInfo = new TemplateMapping[informativeIndexes.length];
                for (int i = 0; i < informativeIndexes.length; i++) {
                    newValues[i] = values[j][informativeIndexes[i]];
                    newMappingInfo[i] = mappingInfo[j][informativeIndexes[i]];
                }
                values[j] = newValues;
                mappingInfo[j] = newMappingInfo;
            }
        } // else {...} no changes.
    }

    public void calculateBestMappingScores() {
        if (!needsToCalculateBestMappingScores) {
            return;
        }
        for (int i = 0; i < templates.size(); i++) {
            for (int j = 0; j < 2; j++) {
                double best = Double.NaN;
                double bestMq = Double.NaN;
                double worst = Double.NaN;
                double worstMq = Double.NaN;
                for (int k = 0; k < haplotypes.size(); k++) {
                    final OptionalDouble score = (j == 0 ? (getMappingInfo(k,i).firstAlignmentScore) : (getMappingInfo(k,i).secondAlignmentScore));
                    if (score.isPresent() && (Double.isNaN(best) || score.getAsDouble() > best)) {
                        best = score.getAsDouble();
                        if (haplotypes.get(k).isContig()) {
                            bestMq = haplotypes.get(k).mappingQuality();
                        }
                    }
                    if (score.isPresent() && (Double.isNaN(worst) || score.getAsDouble() < worst)) {
                        worst = score.getAsDouble();
                        if (haplotypes.get(k).isContig()) {
                            worstMq = haplotypes.get(k).mappingQuality();
                        }
                    }
                }
                bestMappingScorePerFragment[i][j] = best;
                bestMappingFragmentMQ[i][j] = bestMq;
                worstMappingScorePerFragment[i][j] = worst;
                worstMappingScoresFragmentMQ[i][j] = worstMq;
            }
        }
        needsToCalculateBestMappingScores = false;
    }

    public double getWorstAlignmentScore(final int templateIndex, final int fragmentIndex) {
        calculateBestMappingScores();
        final double result = worstMappingScorePerFragment[templateIndex][fragmentIndex];
        return Double.isNaN(result) ? Double.NEGATIVE_INFINITY : result;
    }

    public double getBestAlignmentScore(final int templateIndex, final int fragmentIndex) {
        calculateBestMappingScores();
        final double result = bestMappingScorePerFragment[templateIndex][fragmentIndex];
        return Double.isNaN(result) ? Double.NEGATIVE_INFINITY : result;
    }


    public void applyMissingAlignmentScore(final int template, final int fragment, final double missingAlignmentScore) {
        final OptionalDouble score = OptionalDouble.of(missingAlignmentScore);
        for (int h = 0; h < haplotypes.size(); h++) {
            if (fragment == 0 && !getMappingInfo(h, template).firstAlignmentScore.isPresent())
                getMappingInfo(h, template).firstAlignmentScore = score;
            else if ( fragment == 1 && !getMappingInfo(h, template).secondAlignmentScore.isPresent())
                getMappingInfo(h, template).secondAlignmentScore = score;
        }
    }

    public String toString(final InsertSizeDistribution distr, final int[] refBreakPoints, final int[] altBreakPoints, final Set<String> filterDown) {
        final StringBuilder builder = new StringBuilder((400 * templates.size()) + 20 * (templates.size() * haplotypes().size()));
        builder.append("template");
        for (final SVHaplotype haplotype : haplotypes()) {
            builder.append('\t').append(haplotype.getName());
            if (haplotype.isContig()) {
                final SVContig contig = (SVContig) haplotype;
                final double refScore = contig.getReferenceHaplotypeScore();
                final double altScore = contig.getAlternativeHaplotypeScore();
                final String call = refScore < altScore ? "altHaplotype" : ((altScore < refScore) ? "refHaplotype" : ".");
                final double qual = call.equals("altHaplotype") ? (altScore - refScore) : (call.equals("refHaplotype") ? refScore - altScore : Double.NaN);
                builder.append(':').append(call).append(Double.isNaN(qual) ? "" : String.format(":%.2f", qual));
            }
        }
        for (final Template template : templates) {
            final int templateIndex = indexOf(template);
            if (filterDown != null) {
                if (!filterDown.contains(template.name())) {
                    continue;
                }
            }
            builder.append('\n').append(template.name());
            for (final SVHaplotype haplotype : haplotypes()) {
                final int haplotypeIndex = indexOf(haplotype);
                builder.append('\t');
                final TemplateMapping info = getMappingInfo(haplotypeIndex, templateIndex);
                builder.append(String.format("%.2f", get(haplotypeIndex, templateIndex)))
                        .append(':').append(info.pairOrientation)
                        .append(':').append(info.insertSize.orElse(-1))
                        .append(':').append(String.format("%.2f", info.insertSize.isPresent() ? distr.logProbability(info.insertSize.getAsInt()) : Double.NaN))
                        .append(':').append(String.format("%.2f", info.firstAlignmentScore.orElse(Double.NaN)))
                        .append(':').append(String.format("%.2f", info.secondAlignmentScore.orElse(Double.NaN)));
                if (haplotype.getName().equals("refHaplotype")) {
                    builder.append(':').append(info.crossesBreakPoint(refBreakPoints) ? "1" : "0");
                } else if (haplotype.getName().equals("altHaplotype")) {
                    builder.append(':').append(info.crossesBreakPoint(altBreakPoints) ? "1" : "0");
                }
            }
        }
        return builder.toString();
    }

    public boolean containsNulls() {
        return Arrays.stream(mappingInfo).flatMap(a -> Arrays.stream(a)).anyMatch(Objects::isNull);
    }
}
