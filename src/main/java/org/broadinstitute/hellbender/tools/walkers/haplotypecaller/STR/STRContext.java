package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.STR;

import com.google.gson.Gson;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.utils.GenomeLoc;
import org.broadinstitute.hellbender.utils.GenomeLocParser;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.collections.Pair;
import org.broadinstitute.hellbender.utils.genotyper.*;
import org.broadinstitute.hellbender.utils.sam.GATKSAMRecord;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.List;
import java.util.Map;
import java.util.function.Function;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

/**
 * Captures relevant information about a variant STR site.
 */
public final class STRContext {

    private final VariantContext originalVariantContext;
    private int[] alleleDepth;
    private final STRAlleleSet alleles;
    private final ReadLikelihoods<STRAllele> likelihoods;
    private double gCContent;

    public GenomeLoc getRepeatLocus() {
        return locus;
    }

    private GenomeLoc locus;

    STRContext(final GenomeLoc locus, final VariantContext variantContext, final STRAlleleSet alleles, final ReadLikelihoods<STRAllele> likelihoods) {
        this.locus = locus;
        this.likelihoods = likelihoods;
        this.alleleDepth = likelihoods != null ? calculateAlleleDepth(alleles, likelihoods) :
                (variantContext != null ? calculateAlleleDepth(alleles, variantContext) : null);
        this.originalVariantContext = variantContext;
        this.alleles = alleles;
    }

    private int[] calculateAlleleDepth(final STRAlleleSet alleles, final VariantContext variantContext) {
        if (variantContext.getAlternateAlleles().size() == 1
                && variantContext.getAlternateAllele(0).equals(GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE)) {
            final int[] result = new int[alleles.size()];
            result[alleles.alleleIndex(alleles.getReference())] = variantContext.getGenotype(0).getDP();
            return result;
        } else {
            final int[] result = new int[alleles.size()];
            final int[] ad = variantContext.getGenotype(0).getAD();
            if (ad == null) {
                return null;
            } else {
                final List<Allele> vcAlleles = variantContext.getAlleles();
                for (int i = 0; i < vcAlleles.size(); i++) {
                    final int index = alleles.indexOf(vcAlleles.get(i));
                    if (index >= 0) {
                        result[index] = ad[i];
                    }
                }
                return result;
            }
        }
    }

    private static int[] calculateAlleleDepth(final STRAlleleSet alleles, final ReadLikelihoods<STRAllele> likelihoods) {
        if (likelihoods == null) {
            return null;
        }
        final int[] result = new int[alleles.size()];
        for (final ReadLikelihoods.BestAllele<STRAllele> best : likelihoods.bestAlleles()) {
            if (best.isInformative()) {
                result[alleles.indexOf(best.allele)]++;
            }
        }
        return result;
    }

    public VariantContext getOriginalVariantContext() {
        return originalVariantContext;
    }

    public String toString() {
        return String.format("[STR@%s %s AD='%s']", locus, alleles, getAlleleDepthString());
    }

    public String toLongString() {
        return String.format("[STR@%s %s AD='%s', VC='%s', LK='%s']", locus, alleles, getAlleleDepthString(),
                originalVariantContext, likelihoods);
    }

    private static final Pattern LONG_STRING_PATTERN = Pattern.compile("^\\s*\\[STR@(\\S+)\\s+(\\S+)\\s+AD='(.*?)',\\s+VC='(.*?)', LK='(.*?)'\\]\\s*$");

    public static STRContext fromLongString(final String str,
                                            final GenomeLocParser locParser,
                                            final Function<String, GATKSAMRecord> readFactory) {
        final Matcher matcher = LONG_STRING_PATTERN.matcher(str);
        if (!matcher.matches()) {
            throw new IllegalArgumentException("wrong format: " + str);
        } else {
            final STRAlleleSet alleleSet = STRAlleleSet.fromString(matcher.group(2));
            final String locString = matcher.group(1);
            final GenomeLoc loc = locParser.parseGenomeLoc(locString);
            final ReadLikelihoods<STRAllele> likelihoods = composeReadLikelihoods(matcher.group(5), alleleSet, readFactory);
            final STRAlleleSet alleles = likelihoods.alleles().get(0).set;
            final STRAlleleSet alleles2 = STRAlleleSet.fromString(matcher.group(2));
            if (!alleles.equals(alleles2)) {
                throw new IllegalArgumentException("incompatible allele elements in input: " + alleles + " vs " + alleles2);
            }
            return new STRContext(loc, null, alleles, likelihoods);
        }
    }

    private static ReadLikelihoods<STRAllele> composeReadLikelihoods(final String str, final STRAlleleSet alleleSet, final Function<String, GATKSAMRecord> readFactory) {
        final Map<String, Object> json = new Gson().fromJson(str, Map.class);
        final List<String> alleleStrings = (List<String>) json.get("ALLELES");
        final List<Allele> alleles = alleleStrings.stream()
                .map(s -> {
                    if (s.contains("*")) {
                        return Allele.create(s.replace("*",""), true);
                    } else {
                        return Allele.create(s, false);
                    }
                })
                .collect(Collectors.toList());
        final int[] alleleMap = alleles.stream().mapToInt(alleleSet::indexOf).toArray();
        final Map<String, Map<String, List<Object>>> sampleReadLikelihoods = (Map<String, Map<String, List<Object>>>) json.get("SAMPLES");
        final SampleList samples = new IndexedSampleList(sampleReadLikelihoods.keySet());
        final AlleleList<STRAllele> resultAlleles = new IndexedAlleleList<>(alleleSet);
        final Map<String, List<GATKSAMRecord>> sampleReadMap = sampleReadLikelihoods.entrySet().stream()
                .map(e -> new Pair<>(e.getKey(), e.getValue().keySet().stream().map(readFactory).collect(Collectors.toList())))
                .collect(Collectors.toMap(Pair::getFirst, Pair::getSecond));
        final ReadLikelihoods<STRAllele> result = new ReadLikelihoods<>(samples, resultAlleles, sampleReadMap);
        for (final String sample : sampleReadMap.keySet()) {
            final ReadLikelihoods.Matrix<STRAllele> sampleMatrix = result.sampleMatrix(sample);
            final Map<String, List<Object>> sampleLikelihoods = sampleReadLikelihoods.get(sample);
            for (final GATKSAMRecord read : sampleMatrix.reads()) {
                final List<Object> likelihoods = sampleLikelihoods.get(read.getPairedReadName());
                final int readIndex = sampleMatrix.readIndex(read);
                for (int i = 0; i < resultAlleles.alleleCount(); i++) {
                    final double likelihood = Double.parseDouble(String.valueOf(likelihoods.get(i)));
                    sampleMatrix.set(alleleMap[i], readIndex, likelihood);
                }
            }
        }
        return result;
    }

    public STRAlleleSet getAlleles() {
        return alleles;
    }

    public String getAlleleDepthString() { return Utils.join(",", alleleDepth); };

    public int[] getAlleleDepths() {
        return alleleDepth.clone();
    }

    public ReadLikelihoods<STRAllele> getLikelihoods() {
        return likelihoods;
    }

    public double getMajorAlleleBalance() {
        final long sum = MathUtils.sum(alleleDepth);
        final int max = MathUtils.arrayMax(alleleDepth);
        return max / (double) sum;
    }

    public void setAlleleDepths(int[] alleleDepths) {
        this.alleleDepth = alleleDepths;
    }

    public void setGCContent(final double value) {
        this.gCContent = value;
    }

    public double getGCContent() {
        return gCContent;
    }

    public int getMajorAlleleRepeatCount() {
        return alleles.get(MathUtils.maxElementIndex(alleleDepth)).repeatCount;
    }

    public int getMinimumRepeatCount() {
        return alleles.getMinimumRepeatCount();
    }

    public int getMaximumRepeatCount() {
        return alleles.getMaximumRepeatCount();
    }
}
