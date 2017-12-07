package org.broadinstitute.hellbender.tools;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.MultiVariantWalker;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import scala.Tuple2;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;

@CommandLineProgramProperties(summary = "Combine multiple 10x VCF files (normalized by NormalizeTenxVCFs",
        oneLineSummary = "Combine multiple 10x VCF files (normalized by NormalizeTenxVCFs",
        programGroup = VariantProgramGroup.class)
public class CombineTenxSVVCF extends MultiVariantWalker {

    @Argument(fullName= StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName=StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="The combined VCF output file", optional=false)
    private File outputFile;

    @Argument(fullName="breakpointMergeThreshold", shortName="breakpointMergeThreshold", doc = "Distance at which to merge coherent breakpoints")
    private int breakpointMergeThreshold = 500;

    LinkedList<VariantContext> currentBreakpoints = new LinkedList<>();

    Map<String, String> idMappings = new HashMap<>();

    private LinkedList<VariantContext> beforeForwardVariants = new LinkedList<>();
    private LinkedList<VariantContext> beforeRevVariants = new LinkedList<>();
    private LinkedList<VariantContext> afterForwardVariants = new LinkedList<>();
    private LinkedList<VariantContext> afterRevVariants = new LinkedList<>();
    private VariantContextWriter writer;

    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    public void onTraversalStart() {
        final SortedSet<String> samples = getSamplesForVariants();

        final VCFHeader inputVCFHeader = new VCFHeader(getHeaderForVariants().getMetaDataInInputOrder(), samples);

        writer = createVCFWriter(outputFile);

        final Set<String> sampleNameSet = new IndexedSampleList(samples).asSetOfSamples();
        final VCFHeader vcfHeader = new VCFHeader(inputVCFHeader.getMetaDataInInputOrder(), new TreeSet<>(sampleNameSet));
        vcfHeader.addMetaDataLine(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_COUNT_KEY));
        writer.writeHeader(vcfHeader);

    }

    @Override
    public void apply(final VariantContext variant, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext) {
        if (variant.getAlternateAlleles().size() != 1) {
            throw new GATKException("Can't handle vc with more than one alt");
        }

        if (!variant.getAttribute("SVTYPE").equals("BND")) {
            return;
        }

        System.err.println("processing variant " + variant);
        final BreakendAdjacency bnd = parseBreakendAllele(variant.getAlternateAllele(0).getDisplayString());

        processOutOfScopeVariantLists(variant);

        final LinkedList<VariantContext> variantList = getVariantList(bnd);

        variantList.addLast(variant);

        // if this vc does match a current


    }

    private void processOutOfScopeVariantLists(final VariantContext variant) {
        System.err.println("collecting out of scope variants");
        final List<VariantContext> cliques = new ArrayList<>();
        if (!beforeForwardVariants.isEmpty() && outOfScope(breakpointMergeThreshold, beforeForwardVariants.getLast(), variant)) {
            cliques.addAll(emitCliques(beforeForwardVariants, breakpointMergeThreshold));
        }
        if (!beforeRevVariants.isEmpty() && outOfScope(breakpointMergeThreshold, beforeRevVariants.getLast(), variant)) {
            cliques.addAll(emitCliques(beforeRevVariants, breakpointMergeThreshold));
        }
        if (!afterForwardVariants.isEmpty() && outOfScope(breakpointMergeThreshold, afterForwardVariants.getLast(), variant)) {
            cliques.addAll(emitCliques(afterForwardVariants, breakpointMergeThreshold));
        }
        if (!afterRevVariants.isEmpty() && outOfScope(breakpointMergeThreshold, afterRevVariants.getLast(), variant)) {
            cliques.addAll(emitCliques(afterRevVariants, breakpointMergeThreshold));
        }
        System.err.println("variants to emit: " + cliques.size());
        cliques.sort(IntervalUtils.getDictionaryOrderComparator(getReferenceDictionary()));
        cliques.forEach(System.err::println);
        for (VariantContext vc : cliques) {
            writer.add(vc);
        }
    }

    @Override
    public Object onTraversalSuccess() {
        final List<VariantContext> cliques = new ArrayList<>();
        cliques.addAll(emitCliques(beforeForwardVariants, breakpointMergeThreshold));
        cliques.addAll(emitCliques(beforeRevVariants, breakpointMergeThreshold));
        cliques.addAll(emitCliques(afterForwardVariants, breakpointMergeThreshold));
        cliques.addAll(emitCliques(afterRevVariants, breakpointMergeThreshold));
        cliques.sort(IntervalUtils.getDictionaryOrderComparator(getReferenceDictionary()));
        for (VariantContext vc : cliques) {
            writer.add(vc);
        }
        writer.close();
        return null;
    }

    private List<VariantContext> emitCliques(final LinkedList<VariantContext> variantList, final int breakpointMergeThreshold) {
        final List<VariantContext> cliques = new ArrayList<>();
        while (! variantList.isEmpty()) {
            Set<VariantContext> maxClique = getMaxClique(variantList, breakpointMergeThreshold, getReferenceDictionary());
            cliques.add(emitClique(maxClique));
            variantList.removeAll(maxClique);
        }
        return cliques;
    }

    private VariantContext emitClique(final Set<VariantContext> maxClique) {
        System.err.println("emitting clique " + maxClique);
        if (maxClique.size() == 1) {
            return maxClique.iterator().next();
        }



        final List<VariantContext> sortedByStart = maxClique.stream().sorted(Comparator.comparing(VariantContext::getStart)).collect(Collectors.toList());
        final VariantContext medianStart = sortedByStart.get(Math.round(sortedByStart.size() / 2));
        final List<VariantContext> sortedByBNDPos = maxClique.stream().sorted(Comparator.comparing(v -> parseBreakendAllele(v.getAlternateAllele(0).getDisplayString()).position)).collect(Collectors.toList());
        final VariantContext medianBndPos = sortedByBNDPos.get(Math.round(sortedByBNDPos.size() / 2));
        final BreakendAdjacency bnd = parseBreakendAllele(medianBndPos.getAlternateAllele(0).getDisplayString());
        final String newAltAlleleString = new BreakendAdjacency(medianStart.getReference().getBaseString(), bnd.contig, bnd.position, bnd.before, bnd.revComp).toBNDString();
        final List<Genotype> gts = new ArrayList<>(maxClique.size());

        final List<Allele> alleles = new ArrayList<>(2);
        alleles.add(medianStart.getReference());
        alleles.add(Allele.create(newAltAlleleString));
        VariantContextBuilder builder = new VariantContextBuilder(medianStart)
                .alleles(alleles);

        int ac = 0;

        for (final VariantContext v : maxClique) {
            Genotype oldGenotype = v.getGenotype(0);

            Allele newAllele1 = alleles.get(v.getAlleles().indexOf(oldGenotype.getAllele(0)));
            Allele newAllele2 = alleles.get(v.getAlleles().indexOf(oldGenotype.getAllele(1)));
            gts.add(new GenotypeBuilder(oldGenotype)
                    .alleles(Arrays.asList(newAllele1,
                            newAllele2))
                     .make());
            if (newAllele1.isNonReference()) ac++;
            if (newAllele2.isNonReference()) ac++;
         }

        builder = builder.genotypes(gts);
        builder = builder.attribute(VCFConstants.ALLELE_COUNT_KEY, ac);
        final VariantContext vc = builder.make();
        System.err.println("adding " + vc);

        return vc;
    }

    static Set<VariantContext> getMaxClique(final LinkedList<VariantContext> variantList, final int breakpointMergeThreshold, final SAMSequenceDictionary referenceDictionary) {

        // rectangle graph sweep algorithm as described in http://www.cs.technion.ac.il/~minati/inplace-clique-DAM.pdf
        // imagine the variants are rectangles defined by the uncertainty intervals at the breakends on the x and y axes
//        final SVIntervalTree<VariantContext> xIntervals = new SVIntervalTree<>();
//        for (final VariantContext variantContext : variantList) {
//            xIntervals.put(getXInterval(variantContext, getReferenceDictionary()), variantContext);
//        }
        final Comparator<Locatable> dictionaryOrderComparator = IntervalUtils.getDictionaryOrderComparator(referenceDictionary);
        final SortedMap<Tuple2<Boolean, SimpleInterval>,List<VariantContext>> yBoundaries =
                new TreeMap<>(Comparator.comparing(Tuple2::_2, Collections.reverseOrder(dictionaryOrderComparator)));

        variantList.forEach(v -> {
            final Tuple2<Tuple2<Boolean, SimpleInterval>, Tuple2<Boolean, SimpleInterval>> variantYBoundaries = getYBoundaries(v, breakpointMergeThreshold);
            addToBoundaryMap(yBoundaries, v, variantYBoundaries._1);
            addToBoundaryMap(yBoundaries, v, variantYBoundaries._2);
        });

        final Set<VariantContext> maxClique = getMaxClique(dictionaryOrderComparator, yBoundaries, breakpointMergeThreshold);
        return maxClique;
    }

    private static void addToBoundaryMap(final SortedMap<Tuple2<Boolean, SimpleInterval>, List<VariantContext>> map, final VariantContext v, final Tuple2<Boolean, SimpleInterval> key) {
        if (! map.containsKey(key)) {
            map.put(key, new ArrayList<>());
        }
        map.get(key).add(v);
    }

    private static Set<VariantContext> getMaxClique(final Comparator<Locatable> dictionaryOrderComparator,
                                                    final SortedMap<Tuple2<Boolean, SimpleInterval>, List<VariantContext>> yBoundaries,
                                                    final int breakpointMergeThreshold) {
        Set<VariantContext> maxClique = Collections.emptySet();
        final Set<VariantContext> activeVariants = new HashSet<>();

        for (final Map.Entry<Tuple2<Boolean, SimpleInterval>, List<VariantContext>> next : yBoundaries.entrySet()) {
            // boundary is top
            if (!next.getKey()._1) {
                activeVariants.addAll(next.getValue());
            } else {
                for (final VariantContext ending : next.getValue()) {
                    final SimpleInterval xInterval = getXInterval(ending, breakpointMergeThreshold);
                    TreeMap<Tuple2<Boolean, SimpleInterval>, List<VariantContext>> overlappingXIntervalEndpoints =
                            new TreeMap<>(Comparator.comparing(Tuple2::_2, dictionaryOrderComparator));
                    for (final VariantContext active : activeVariants) {
                        if (active == ending) continue;
                        final SimpleInterval activeX = getXInterval(active, breakpointMergeThreshold);
                        if (activeX.overlaps(xInterval)) {
                            final int start = Math.max(xInterval.getStart(), activeX.getStart());
                            final int end = Math.min(xInterval.getEnd(), activeX.getEnd());
                            final SimpleInterval newStart =
                                    new SimpleInterval(activeX.getContig(), start, start);
                            final SimpleInterval newEnd =
                                    new SimpleInterval(activeX.getContig(), end, end);
                            final Tuple2<Boolean, SimpleInterval> startKey = new Tuple2<>(true, newStart);
                            final Tuple2<Boolean, SimpleInterval> endKey = new Tuple2<>(false, newEnd);
                            addToBoundaryMap(overlappingXIntervalEndpoints, active, startKey);
                            addToBoundaryMap(overlappingXIntervalEndpoints, active, endKey);
                        }
                    }
                    final Set<VariantContext> currentClique = new HashSet<>();
                    currentClique.add(ending);
                    if (currentClique.size() > maxClique.size()) {
                        maxClique = new HashSet<>(currentClique);
                    }
                    for (final Map.Entry<Tuple2<Boolean, SimpleInterval>, List<VariantContext>> xBoundaries : overlappingXIntervalEndpoints.entrySet()) {
                        if (xBoundaries.getKey()._1) {
                            // adding to the clique
                            currentClique.addAll(xBoundaries.getValue());
                            if (currentClique.size() > maxClique.size()) {
                                maxClique = new HashSet<>(currentClique);
                            }
                        } else {
                            currentClique.removeAll(xBoundaries.getValue());
                        }
                    }
                    activeVariants.remove(ending);
                }
            }
        }
        return maxClique;
    }

    private static SimpleInterval getXInterval(final VariantContext variantContext, final int breakpointMergeThreshold) {
        return new SimpleInterval(variantContext.getContig(),
                Math.max(1, variantContext.getStart() - breakpointMergeThreshold),
                variantContext.getStart() + breakpointMergeThreshold);
    }

    private static Tuple2<Tuple2<Boolean, SimpleInterval>, Tuple2<Boolean, SimpleInterval>> getYBoundaries(final VariantContext variantContext, final int breakpointMergeThreshold) {
        final BreakendAdjacency breakendAdjacency = parseBreakendAllele(variantContext.getAlternateAllele(0).getDisplayString());
        final SimpleInterval yStart = new SimpleInterval(breakendAdjacency.contig,
                Math.max(1,breakendAdjacency.position - breakpointMergeThreshold),
                Math.max(1,breakendAdjacency.position - breakpointMergeThreshold));
        final SimpleInterval yEnd = new SimpleInterval(breakendAdjacency.contig, breakendAdjacency.position + breakpointMergeThreshold, breakendAdjacency.position + breakpointMergeThreshold);

        return new Tuple2<>(new Tuple2<>(true, yStart), new Tuple2<>(false, yEnd));
    }

    private boolean outOfScope(final int breakpointMergeThreshold, final VariantContext last, final VariantContext variant) {
        return !last.getContig().equals(variant.getContig()) || Math.abs(last.getStart() - variant.getStart()) > 2 * breakpointMergeThreshold;
    }

    private LinkedList<VariantContext> getVariantList(final BreakendAdjacency bnd) {
        if (bnd.before) {
            if (bnd.revComp) {
                return beforeRevVariants;
            } else {
                return beforeForwardVariants;
            }
        } else {
            if (bnd.revComp) {
                return afterRevVariants;
            } else {
                return afterForwardVariants;
            }
        }
    }

    static boolean concordant(final int breakpointMergeThreshold, final VariantContext vc1, final VariantContext vc2) {
        if (vc1.getAttribute("SVTYPE").equals("BND") && vc2.getAttribute("SVTYPE").equals("BND")) {
            // todo
            Utils.validate(vc1.getAlternateAlleles().size() == 1 && vc2.getAlternateAlleles().size() == 1, "can't handle BND variants with more than one alt");

            final String vc1Contig = vc1.getContig();
            final int vc1Pos = vc1.getStart();
            final BreakendAdjacency bnd1 = parseBreakendAllele(vc1.getAlternateAllele(0).getDisplayString());

            final String vc2Contig = vc2.getContig();
            final int vc2Pos = vc2.getStart();
            final BreakendAdjacency bnd2 = parseBreakendAllele(vc2.getAlternateAllele(0).getDisplayString());

            // todo: local bases?
            return (vc1Contig.equals(vc2Contig) &&
                    Math.abs(vc1Pos - vc2Pos) <= breakpointMergeThreshold &&
                    bnd1.before == bnd2.before &&
                    bnd1.revComp == bnd2.revComp &&
                    bnd1.contig.equals(bnd2.contig) &&
                    Math.abs(bnd1.position - bnd2.position) <= breakpointMergeThreshold);

        } else {
            return false;
        }
    }


    static BreakendAdjacency parseBreakendAllele(final String breakendAllele) {
        final String localBases;
        final boolean before;
        final boolean revComp;
        final SimpleInterval adjacentPos;

        final String[] fields = breakendAllele.split("(?<=(\\[|\\])|(?=(\\[|\\])))");
        if (fields[0].equals("[") || fields[0].equals("]")) {
            before = true;
            if (fields[0].equals("]")) {
                revComp = false;
            } else {
                revComp = true;
            }
            adjacentPos = new SimpleInterval(fields[1]);
            localBases = fields[fields.length - 1];
        } else if (fields[fields.length -1].equals("[") || fields[fields.length -1].equals("]")) {
            before = false;
            if (fields[fields.length - 1].equals("[")) {
                revComp = false;
            } else {
                revComp = true;
            }
            adjacentPos = new SimpleInterval(fields[fields.length - 2]);
            localBases = fields[0];
        } else {
            // todo: unpaired breakends
            throw new GATKException("Don't know how to parse BND allele " + breakendAllele);
        }
        return new BreakendAdjacency(localBases, adjacentPos.getContig(), adjacentPos.getStart(), before, revComp);
    }

    static class BreakendAdjacency {
        final String localBases;
        final String contig;
        final int position;
        final boolean before;
        final boolean revComp;

        public BreakendAdjacency(final String localBases, final String contig, final int position, final boolean before, final boolean revComp) {
            this.localBases = localBases;
            this.contig = contig;
            this.position = position;
            this.before = before;
            this.revComp = revComp;
        }

        public String toBNDString() {
            if (before) {
                final String bracket;
                if (revComp) {
                    bracket = "[";
                } else {
                    bracket = "]";
                }
                return bracket +
                        contig +
                        ":" +
                        position +
                        bracket +
                        localBases;
            } else {
                final String bracket;
                if (revComp) {
                    bracket = "]";
                } else {
                    bracket = "[";
                }
                return localBases +
                        bracket +
                        contig +
                        ":" +
                        position +
                        bracket;
            }
        }
    }

}
