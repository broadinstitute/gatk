package org.broadinstitute.hellbender.tools.spark.sv;

import com.google.cloud.dataflow.sdk.repackaged.com.google.common.base.Functions;
import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import org.apache.commons.collections4.CollectionUtils;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeLikelihoodCalculator;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeLikelihoodCalculators;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.IndexedAlleleList;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;

import java.io.Serializable;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Created by valentin on 5/18/17.
 */
public class TemplateHaplotypeScoreTable implements Serializable {

    private static final long serialVersionUID = 1L;

    private final TemplateMappingInformation[][] mappingInfo;

    private final double[][] values;

    private final List<Template> templates;

    private final Map<String, Integer> templateIndex;

    private final List<Haplotype> haplotypes;


    public int maximumInsertSize() {
        int result = 0;
        for (int i = 0; i < mappingInfo.length; i++) {
            for (int j = 0; j < mappingInfo[i].length; j++) {
                if (mappingInfo[i][j] != null && mappingInfo[i][j].insertSize.isPresent() && mappingInfo[i][j].insertSize.getAsInt() > result) {
                    result =  mappingInfo[i][j].insertSize.getAsInt();
                }
            }
        }
        return result;
    }

    public double minimumAlignmentScore() {
        double result = 0;
        for (final TemplateMappingInformation[] row : mappingInfo) {
            for (final TemplateMappingInformation ti : row) {
                if (ti == null)
                    continue;
                if (ti.firstAlignmentScore.isPresent() && ti.firstAlignmentScore.getAsDouble() > result)
                    result = ti.firstAlignmentScore.getAsDouble();
                if (ti.secondAlignmentScore.isPresent() && ti.secondAlignmentScore.getAsDouble() > result)
                    result = ti.secondAlignmentScore.getAsDouble();
            }
        }
        return result;
    }

    public TemplateHaplotypeScoreTable(final Iterable<Template> templates, final Iterable<Haplotype> haplotypes)
    {
        this.templates = CollectionUtils.collect(templates, t -> t, new ArrayList<>(1000));
        this.haplotypes = CollectionUtils.collect(haplotypes, t -> t, new ArrayList<>());
        values = new double[this.haplotypes.size()][this.templates.size()];
        mappingInfo = new TemplateMappingInformation[this.haplotypes.size()][this.templates.size()];
        this.templateIndex = composeTemplateIndex(this.templates);
    }

    private Map<String,Integer> composeTemplateIndex(final List<Template> templates) {
        final Map<String, Integer> result = new HashMap<>(templates.size());
        for (int i = 0; i < templates.size(); i++) {
            result.put(templates.get(i).name(), i);
        }
        return result;
    }

    public double get(final Haplotype allele, final Template template) {
        return get(indexOf(allele), indexOf(template));
    }

    public int indexOf(final Template template) {
        Utils.nonNull(template);
        return templateIndex.getOrDefault(template.name(), -1);
    }

    public int indexOf(final Haplotype allele) {
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

    public TemplateMappingInformation getMappingInfo(final int alleleIndex, final int templateIndex) {
        Utils.validIndex(alleleIndex, haplotypes.size());
        Utils.validIndex(templateIndex, templates.size());
        return mappingInfo[alleleIndex][templateIndex];
    }

    public void setMappingInfo(final int alleleIndex, final int templateIndex, final TemplateMappingInformation mappingInformation) {
        Utils.validIndex(alleleIndex, haplotypes.size());
        Utils.validIndex(templateIndex, templates.size());
        mappingInfo[alleleIndex][templateIndex] = mappingInformation;
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

    public List<Haplotype> haplotypes() {
        return haplotypes;
    }

    public String toString() {
        final StringBuilder sb = new StringBuilder();
        for (int i = 0; i < numberOfHaplotypes(); i++) {
            sb.append(Arrays.toString(values[i]));
            sb.append(", ");
        }
        if (sb.length() > 0) {
            sb.setLength(sb.length() - 1);
        }
        return sb.toString();
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

    public void dropUninformativeTemplates() {
        final int[] informativeIndexes = informativeTemplateIndexes();
        if (informativeIndexes.length == 0) {
            templates.clear();
            templateIndex.clear();
            for (int j = 0; j < values.length; j++) {
                values[j] = new double[0];
            }
        } else if (informativeIndexes.length != templates.size()) {
            final List<Template> newTemplates = new ArrayList<>(informativeIndexes.length);
            templateIndex.clear();
            for (final int informativeIndexe : informativeIndexes) {
                newTemplates.add(templates.get(informativeIndexe));
                templateIndex.put(newTemplates.get(newTemplates.size() - 1).name(), newTemplates.size() - 1);
            }
            templates.clear();
            templates.addAll(newTemplates);
            for (int j = 0; j < haplotypes.size(); j++) {
                final double[] newValues = new double[informativeIndexes.length];
                for (int i = 0; i < informativeIndexes.length; i++) {
                    newValues[i] = values[j][informativeIndexes[i]];
                }
                values[j] = newValues;
            }
        } // else {...} no changes.
    }

    public GenotypeLikelihoods calculateGenotypeLikelihoods(final int ploidy) {
        final ReadLikelihoods<Haplotype> likelihoods = new ReadLikelihoods<>(SampleList.singletonSampleList("the-sample"),
                new IndexedAlleleList<>(haplotypes), Collections.singletonMap("the-sample", templates.stream().map(t -> t.fragments().get(0).toUnmappedRead(null, false)).collect(Collectors.toList())));
        final LikelihoodMatrix<Haplotype> matrix = likelihoods.sampleMatrix(0);
        matrix.setAll(values);
        final GenotypeLikelihoodCalculator calculator = new GenotypeLikelihoodCalculators().getInstance(ploidy, haplotypes.size());
        return calculator.genotypeLikelihoods(matrix);
    }

    public double[] getRow(final int i) {
        return values[Utils.validIndex(i, values.length)].clone();
    }
}
