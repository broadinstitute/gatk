package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.tribble.readers.TabixReader;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.DiscordantPairEvidence;
import org.broadinstitute.hellbender.tools.sv.SVCallRecordUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Reader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.EnumSet;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * Refines complex structural variants and translocations using PE and RD evidence
 */
@CommandLineProgramProperties(
        summary = "Reassess PE and RD support for complex structural variants and translocations",
        oneLineSummary = "Refine complex structural variants",
        programGroup = StructuralVariantDiscoveryProgramGroup.class
)
@DocumentedFeature
@BetaFeature
public final class RefineComplexVariants extends VariantWalker {

    static final String BATCH_SAMPLE_LISTS_LONG_NAME = "batch-sample-lists";
    static final String DISCORDANT_PAIR_FILES_LONG_NAME = "discordant-pairs-files";
    static final String DEPTH_DEL_BEDS_LONG_NAME = "depth-del-beds";
    static final String DEPTH_DUP_BEDS_LONG_NAME = "depth-dup-beds";
    static final String MIN_PE_CPX_LONG_NAME = "min-pe-cpx";
    static final String MIN_PE_CTX_LONG_NAME = "min-pe-ctx";

    static final String UNRESOLVED_FILTER = "UNRESOLVED";
    static final String SOURCE_ATTRIBUTE = "SOURCE";

    private static final int PE_FLANK_BACK = 1000;
    private static final int PE_FLANK_FRONT = 100;
    private static final int LARGE_DEPTH_INTERVAL_THRESHOLD = 5000;
    private static final double MIN_DEPTH_COVERAGE = 0.5;

    private static final EnumSet<GATKSVVCFConstants.ComplexVariantSubtype> STRUCTURALLY_UNRESOLVED_SUBTYPES =
            EnumSet.of(GATKSVVCFConstants.ComplexVariantSubtype.piDUP_FR,
                    GATKSVVCFConstants.ComplexVariantSubtype.piDUP_RF);

    private static final EnumSet<GATKSVVCFConstants.ComplexVariantSubtype> INVERSION_CNV_SUBTYPES =
            EnumSet.of(GATKSVVCFConstants.ComplexVariantSubtype.delINV,
                    GATKSVVCFConstants.ComplexVariantSubtype.INVdel,
                    GATKSVVCFConstants.ComplexVariantSubtype.dupINV,
                    GATKSVVCFConstants.ComplexVariantSubtype.INVdup,
                    GATKSVVCFConstants.ComplexVariantSubtype.delINVdel,
                    GATKSVVCFConstants.ComplexVariantSubtype.delINVdup,
                    GATKSVVCFConstants.ComplexVariantSubtype.dupINVdel,
                    GATKSVVCFConstants.ComplexVariantSubtype.dupINVdup);

    @Argument(
            doc = "Output VCF",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private GATKPath outputVcf;

    @Argument(
            doc = "Batch sample membership files. Order must match the PE and depth inputs.",
            fullName = BATCH_SAMPLE_LISTS_LONG_NAME,
            suppressFileExpansion = true
    )
    private List<GATKPath> batchSampleLists;

    @Argument(
            doc = "Batch-level discordant-pair evidence files (*.pe.txt.gz). Order must match the sample lists.",
            fullName = DISCORDANT_PAIR_FILES_LONG_NAME,
            suppressFileExpansion = true
    )
    private List<GATKPath> discordantPairFiles;

    @Argument(
            doc = "Batch-level raw depth DEL BED files. Order must match the sample lists.",
            fullName = DEPTH_DEL_BEDS_LONG_NAME,
            suppressFileExpansion = true
    )
    private List<GATKPath> depthDelBeds;

    @Argument(
            doc = "Batch-level raw depth DUP BED files. Order must match the sample lists.",
            fullName = DEPTH_DUP_BEDS_LONG_NAME,
            suppressFileExpansion = true
    )
    private List<GATKPath> depthDupBeds;

    @Argument(
            doc = "Minimum discordant pairs required at each breakpoint for CPX and INS-with-INV refinement.",
            fullName = MIN_PE_CPX_LONG_NAME,
            optional = true,
            minValue = 0
    )
    private int minPeCpx = 3;

    @Argument(
            doc = "Minimum discordant pairs required at each breakpoint for CTX refinement.",
            fullName = MIN_PE_CTX_LONG_NAME,
            optional = true,
            minValue = 0
    )
    private int minPeCtx = 3;

    private VariantContextWriter writer;
    private SAMSequenceDictionary dictionary;
    private Map<String, Integer> sampleToBatchIndex;
    private List<FeatureDataSource<DiscordantPairEvidence>> discordantPairSources;
    private List<TabixReader> delDepthReaders;
    private List<TabixReader> dupDepthReaders;

    @Override
    public void onTraversalStart() {
        validateInputs();

        final VCFHeader inputHeader = getHeaderForVariants();
        dictionary = getBestAvailableSequenceDictionary();
        Utils.validate(!inputHeader.getSampleNamesInOrder().isEmpty(), "Input VCF must contain sample genotypes");

        sampleToBatchIndex = loadSampleToBatchIndex(batchSampleLists);
        discordantPairSources = discordantPairFiles.stream()
            .map(path -> new FeatureDataSource<DiscordantPairEvidence>(path.toString(), null,
                0, DiscordantPairEvidence.class))
                .collect(Collectors.toList());
        delDepthReaders = createEmptyDepthReaderList(depthDelBeds.size());
        dupDepthReaders = createEmptyDepthReaderList(depthDupBeds.size());

        writer = createVCFWriter(outputVcf);
        writer.writeHeader(createHeader(inputHeader));
    }

    @Override
    public void apply(final VariantContext variant,
                      final ReadsContext readsContext,
                      final ReferenceContext referenceContext,
                      final FeatureContext featureContext) {
        // Each record follows the same workflow: determine the PE/RD evidence plan, revise carrier
        // genotypes if needed, then rewrite the record into the workflow-compatible CPX representation.
        final EvaluationPlan plan = createEvaluationPlan(variant, dictionary);
        if (plan.variantType == EvaluatedVariantType.NONE && !requiresFormattingTransform(variant)) {
            writer.add(variant);
            return;
        }

        final List<String> carrierSamples = getCarrierSamples(variant);
        final Set<String> updatedFilters = new LinkedHashSet<>(variant.getFilters());
        final boolean preexistingUnresolved = updatedFilters.contains(UNRESOLVED_FILTER);
        GenotypesContext updatedGenotypes = GenotypesContext.copy(variant.getGenotypes());

        if (plan.structuralUnresolved) {
            updatedFilters.add(UNRESOLVED_FILTER);
        } else if ((!preexistingUnresolved || shouldReevaluatePreexistingUnresolved(variant))
            && !plan.queries.isEmpty() && !carrierSamples.isEmpty()) {
            final Map<String, int[]> peCounts = collectDiscordantPairCounts(plan.queries, carrierSamples);
            final Map<String, List<Boolean>> depthSupport = plan.requiresDepthEvidence
                    ? collectDepthSupport(variant, carrierSamples)
                    : Collections.emptyMap();

            int unresolvedCarrierCount = 0;
            final ArrayList<Genotype> revisedGenotypes = new ArrayList<>(updatedGenotypes.size());
            for (final Genotype genotype : updatedGenotypes) {
                if (!SVCallRecordUtils.isAltGenotype(genotype)) {
                    revisedGenotypes.add(genotype);
                    continue;
                }

                final int[] counts = peCounts.getOrDefault(genotype.getSampleName(), new int[]{0, 0});
                final CarrierRefinement refinement = evaluateCarrierSupport(
                        counts[0],
                        counts[1],
                        plan.variantType == EvaluatedVariantType.CTX ? minPeCtx : minPeCpx,
                        depthSupport.getOrDefault(genotype.getSampleName(), Collections.emptyList()));
                revisedGenotypes.add(refinement.reviseGenotype ? makeNoCallGenotype(genotype) : genotype);
                if (refinement.countTowardsUnresolved) {
                    unresolvedCarrierCount++;
                }
            }

            updatedGenotypes = GenotypesContext.create(revisedGenotypes);
            if (plan.variantType == EvaluatedVariantType.CPX
                    && unresolvedCarrierCount * 2 >= carrierSamples.size()) {
                updatedFilters.add(UNRESOLVED_FILTER);
            }
        }

        final Map<String, Object> updatedAttributes = new LinkedHashMap<>(variant.getAttributes());
        final Allele rewrittenAlternateAllele = applyFormattingTransform(variant, updatedAttributes);
        if (rewrittenAlternateAllele != null) {
            updatedGenotypes = remapGenotypesToAlternateAllele(updatedGenotypes, rewrittenAlternateAllele);
        }

        final VariantContextBuilder builder = new VariantContextBuilder(variant)
                .attributes(updatedAttributes)
                .genotypes(updatedGenotypes);
        if (updatedFilters.isEmpty()) {
            builder.passFilters();
        } else {
            builder.filters(updatedFilters);
        }
        if (rewrittenAlternateAllele != null) {
            builder.alleles(Arrays.asList(variant.getReference(), rewrittenAlternateAllele));
        }
        writer.add(builder.make());
    }

    @Override
    public void closeTool() {
        if (writer != null) {
            writer.close();
        }
        if (discordantPairSources != null) {
            discordantPairSources.forEach(FeatureDataSource::close);
        }
        if (delDepthReaders != null) {
            delDepthReaders.stream().filter(reader -> reader != null).forEach(TabixReader::close);
        }
        if (dupDepthReaders != null) {
            dupDepthReaders.stream().filter(reader -> reader != null).forEach(TabixReader::close);
        }
    }

    private void validateInputs() {
        Utils.nonEmpty(batchSampleLists, "At least one batch sample list is required");
        Utils.nonEmpty(discordantPairFiles, "At least one PE evidence file is required");
        Utils.nonEmpty(depthDelBeds, "At least one DEL raw depth BED is required");
        Utils.nonEmpty(depthDupBeds, "At least one DUP raw depth BED is required");

        final int batchCount = batchSampleLists.size();
        Utils.validate(batchCount == discordantPairFiles.size(),
                "PE file count must match the number of batch sample lists");
        Utils.validate(batchCount == depthDelBeds.size(),
                "DEL depth BED count must match the number of batch sample lists");
        Utils.validate(batchCount == depthDupBeds.size(),
                "DUP depth BED count must match the number of batch sample lists");
    }

    private VCFHeader createHeader(final VCFHeader inputHeader) {
        final VCFHeader header = new VCFHeader(inputHeader);
        if (!header.hasFilterLine(UNRESOLVED_FILTER)) {
            header.addMetaDataLine(new VCFFilterHeaderLine(
                    UNRESOLVED_FILTER,
                    "Variant remained unresolved after complex-variant refinement"));
        }
        return header;
    }

    private Map<String, Integer> loadSampleToBatchIndex(final List<GATKPath> sampleLists) {
        final Map<String, Integer> mapping = new HashMap<>();
        for (int batchIndex = 0; batchIndex < sampleLists.size(); batchIndex++) {
            try (BufferedReader reader = openBufferedReader(sampleLists.get(batchIndex).toString())) {
                String line;
                while ((line = reader.readLine()) != null) {
                    final String sample = line.trim();
                    if (sample.isEmpty()) {
                        continue;
                    }
                    final Integer previous = mapping.put(sample, batchIndex);
                    if (previous != null) {
                        throw new UserException.BadInput(
                                "Sample " + sample + " appears in multiple batch sample lists");
                    }
                }
            } catch (final IOException e) {
                throw new UserException.CouldNotReadInputFile(sampleLists.get(batchIndex).toString(), e);
            }
        }
        return mapping;
    }

    private List<TabixReader> createEmptyDepthReaderList(final int batchCount) {
        return new ArrayList<>(Collections.nCopies(batchCount, null));
    }

    private TabixReader getOrOpenDepthReader(final List<TabixReader> readers,
                                             final List<GATKPath> paths,
                                             final int batchIndex,
                                             final String label) {
        TabixReader reader = readers.get(batchIndex);
        if (reader != null) {
            return reader;
        }

        final String path = paths.get(batchIndex).toString();
        try {
            reader = new TabixReader(path);
        } catch (final IOException e) {
            throw new UserException.CouldNotReadInputFile(
                    path,
                    "Could not open indexed " + label + " BED file (expected bgzip + .tbi): " + path,
                    e);
        }
        readers.set(batchIndex, reader);
        return reader;
    }

    private double queryDepthCoverage(final TabixReader reader,
                                      final String sample,
                                      final String contig,
                                      final int queryStart,
                                      final int queryEnd,
                                      final String pathForError) {
        final int regionStart = Math.max(1, queryStart - 1); // TabixReader uses 0-based half-open coordinates
        final int regionEnd = Math.max(regionStart, queryEnd);
        final String region = contig + ":" + regionStart + "-" + regionEnd;

        final List<DepthInterval> intervals = new ArrayList<>();
        try {
            final TabixReader.Iterator iterator = reader.query(region);
            if (iterator == null) {
                return 0.0;
            }

            String line;
            while ((line = iterator.next()) != null) {
                if (line.isEmpty() || line.charAt(0) == '#') {
                    continue;
                }
                final String[] fields = line.split("\t", -1);
                if (fields.length < 5) {
                    throw new UserException.BadInput("Malformed raw depth BED line in " + pathForError + ": " + line);
                }
                if (!sample.equals(fields[4])) {
                    continue;
                }

                final int start = Integer.parseInt(fields[1]);
                final int end = Integer.parseInt(fields[2]);
                if (end > queryStart && start < queryEnd) {
                    intervals.add(new DepthInterval(start, end));
                }
            }
        } catch (final IOException e) {
            throw new UserException.CouldNotReadInputFile(pathForError, e);
        } catch (final NumberFormatException e) {
            throw new UserException.BadInput("Malformed coordinate in raw depth BED line from " + pathForError, e);
        }

        intervals.sort(Comparator.comparingInt(interval -> interval.start));
        return computeCoverageFraction(intervals, regionStart, regionEnd);
    }

    private static BufferedReader openBufferedReader(final String path) throws IOException {
        final Reader reader = IOUtils.makeReaderMaybeGzipped(IOUtils.getPath(path));
        return new BufferedReader(reader);
    }

    private Map<String, int[]> collectDiscordantPairCounts(final List<DiscordantPairQuery> queries,
                                                           final List<String> carrierSamples) {
        final Map<String, int[]> counts = new HashMap<>();
        final Map<Integer, Set<String>> carriersByBatch = new HashMap<>();
        for (final String sample : carrierSamples) {
            final Integer batchIndex = sampleToBatchIndex.get(sample);
            if (batchIndex == null) {
                throw new UserException.BadInput(
                        "Carrier sample " + sample + " is missing from the batch sample lists");
            }
            carriersByBatch.computeIfAbsent(batchIndex, ignored -> new HashSet<>()).add(sample);
            counts.put(sample, new int[queries.size()]);
        }

        for (final Map.Entry<Integer, Set<String>> entry : carriersByBatch.entrySet()) {
            final FeatureDataSource<DiscordantPairEvidence> source = discordantPairSources.get(entry.getKey());
            for (int queryIndex = 0; queryIndex < queries.size(); queryIndex++) {
                final Map<String, Integer> queryCounts = countMatchingEvidence(source, queries.get(queryIndex), entry.getValue());
                for (final String sample : entry.getValue()) {
                    counts.get(sample)[queryIndex] = queryCounts.getOrDefault(sample, 0);
                }
            }
        }
        return counts;
    }

    private Map<String, Integer> countMatchingEvidence(final FeatureDataSource<DiscordantPairEvidence> source,
                                                       final DiscordantPairQuery query,
                                                       final Set<String> carrierSamples) {
        final Map<String, Integer> counts = new HashMap<>();
        final SimpleInterval interval = makeQueryInterval(query.startContig, query.startMin, query.startMax);
        source.query(interval).forEachRemaining(evidence -> {
            if (!carrierSamples.contains(evidence.getSample())) {
                return;
            }
            if (!evidence.getContig().equals(query.startContig) || !evidence.getEndContig().equals(query.endContig)) {
                return;
            }
            if (evidence.getStartStrand() != query.startStrand || evidence.getEndStrand() != query.endStrand) {
                return;
            }
            if (evidence.getStart() < query.startMin || evidence.getStart() > query.startMax) {
                return;
            }
            if (evidence.getEndPosition() <= query.endLowerExclusive || evidence.getEndPosition() >= query.endUpperExclusive) {
                return;
            }
            counts.merge(evidence.getSample(), 1, Integer::sum);
        });
        return counts;
    }

    private SimpleInterval makeQueryInterval(final String contig, final int start, final int end) {
        final SAMSequenceRecord record = dictionary.getSequence(contig);
        if (record == null) {
            throw new UserException.BadInput(
                    "Contig " + contig + " is not present in the input sequence dictionary");
        }
        final int trimmedStart = Math.max(1, start);
        final int trimmedEnd = Math.min(record.getSequenceLength(), Math.max(trimmedStart, end));
        return new SimpleInterval(contig, trimmedStart, trimmedEnd);
    }

    private Map<String, List<Boolean>> collectDepthSupport(final VariantContext variant,
                                                           final List<String> carrierSamples) {
        final List<SVSegment> depthIntervals = collectLargeDepthIntervals(variant);
        if (depthIntervals.isEmpty() || carrierSamples.isEmpty()) {
            return Collections.emptyMap();
        }

        final Map<String, List<Boolean>> support = new HashMap<>();
        for (final String sample : carrierSamples) {
            final Integer batchIndex = sampleToBatchIndex.get(sample);
            if (batchIndex == null) {
                throw new UserException.BadInput(
                        "Carrier sample " + sample + " is missing from the batch sample lists");
            }

            final List<Boolean> intervalSupport = new ArrayList<>(depthIntervals.size());
            for (final SVSegment interval : depthIntervals) {
                final boolean isDeletion = interval.getIntervalSVType() == GATKSVVCFConstants.StructuralVariantAnnotationType.DEL;
                final List<TabixReader> readers = isDeletion ? delDepthReaders : dupDepthReaders;
                final List<GATKPath> paths = isDeletion ? depthDelBeds : depthDupBeds;
                final TabixReader reader = getOrOpenDepthReader(readers, paths, batchIndex, isDeletion ? "DEL" : "DUP");
                final double coverage = queryDepthCoverage(
                    reader,
                        sample,
                        interval.getContig(),
                        interval.getStart(),
                    interval.getEnd(),
                    paths.get(batchIndex).toString());
                intervalSupport.add(coverage > MIN_DEPTH_COVERAGE);
            }
            support.put(sample, intervalSupport);
        }
        return support;
    }

    /**
     * Collects large DEL and DUP intervals from complex events to be evaluated for depth support
     */
    static List<SVSegment> collectLargeDepthIntervals(final VariantContext variant) {
        final List<SVSegment> intervals = new ArrayList<>();
        addLargeDepthIntervals(variant.getAttributeAsStringList(GATKSVVCFConstants.CPX_INTERVALS, null), intervals);

        final String source = variant.getAttributeAsString(SOURCE_ATTRIBUTE, null);
        if (source != null && !source.equals(VCFConstants.MISSING_VALUE_v4)) {
            addLargeDepthIntervals(Arrays.asList(source.split(",")), intervals);
        }
        return intervals;
    }

    private static void addLargeDepthIntervals(final List<String> encodedIntervals, final List<SVSegment> intervals) {
        for (final SVSegment segment : SVAnnotateEngine.parseComplexIntervals(encodedIntervals)) {
            final GATKSVVCFConstants.StructuralVariantAnnotationType intervalType = segment.getIntervalSVType();
            if ((intervalType == GATKSVVCFConstants.StructuralVariantAnnotationType.DEL
                    || intervalType == GATKSVVCFConstants.StructuralVariantAnnotationType.DUP)
                    && segment.getEnd() - segment.getStart() >= LARGE_DEPTH_INTERVAL_THRESHOLD) {
                intervals.add(segment);
            }
        }
    }

    /**
     * Adjudicates carrier support based on PE and depth evidence
     */
    static CarrierRefinement evaluateCarrierSupport(final int firstBreakpointCount,
                                                    final int secondBreakpointCount,
                                                    final int minPe,
                                                    final List<Boolean> depthSupport) {
        final boolean hasDepthEvidence = !depthSupport.isEmpty();
        final boolean allDepthIntervalsPass = hasDepthEvidence && depthSupport.stream().allMatch(Boolean::booleanValue);
        final boolean anyDepthIntervalFails = hasDepthEvidence && depthSupport.stream().anyMatch(pass -> !pass);
        final boolean firstHasAny = firstBreakpointCount > 0;
        final boolean secondHasAny = secondBreakpointCount > 0;
        final boolean firstPassesThreshold = firstBreakpointCount >= minPe;
        final boolean secondPassesThreshold = secondBreakpointCount >= minPe;

        if (firstPassesThreshold && secondPassesThreshold) {
            return new CarrierRefinement(anyDepthIntervalFails, false);
        }
        if (allDepthIntervalsPass) {
            return CarrierRefinement.KEEP;
        }

        return new CarrierRefinement(true, firstHasAny ^ secondHasAny);
    }

    /**
     * Converts insertions with inverted source to dDUP
     */
    static Allele applyFormattingTransform(final VariantContext variant, final Map<String, Object> attributes) {
        final GATKSVVCFConstants.StructuralVariantAnnotationType variantType =
                SVCallRecordUtils.inferStructuralVariantType(variant);
        final GATKSVVCFConstants.ComplexVariantSubtype complexSubtype = SVCallRecordUtils.getComplexSubtype(variant);
        final String source = variant.getAttributeAsString(SOURCE_ATTRIBUTE, null);

        if (variantType == GATKSVVCFConstants.StructuralVariantAnnotationType.CPX
                && (complexSubtype == GATKSVVCFConstants.ComplexVariantSubtype.dDUP
                || complexSubtype == GATKSVVCFConstants.ComplexVariantSubtype.dDUP_iDEL
                || complexSubtype == GATKSVVCFConstants.ComplexVariantSubtype.INS_iDEL)) {
            if (complexSubtype == GATKSVVCFConstants.ComplexVariantSubtype.INS_iDEL && source != null) {
                final List<String> intervals = new ArrayList<>(variant.getAttributeAsStringList(GATKSVVCFConstants.CPX_INTERVALS, null));
                intervals.add(source);
                attributes.put(GATKSVVCFConstants.CPX_INTERVALS, intervals);
            }
            attributes.remove(SOURCE_ATTRIBUTE);
            return null;
        }

        if (variantType != GATKSVVCFConstants.StructuralVariantAnnotationType.INS || source == null || !source.contains("INV")) {
            return null;
        }

        final List<SVSegment> sourceSegments = SVAnnotateEngine.parseComplexIntervals(Arrays.asList(source));
        final SVSegment inversionSegment = getRequiredSegment(
                sourceSegments,
                GATKSVVCFConstants.StructuralVariantAnnotationType.INV,
                "Insertion-with-inversion record is missing an INV source interval");
        final int deletionLength = variant.getEnd() - variant.getStart();
        final GATKSVVCFConstants.ComplexVariantSubtype rewrittenSubtype = deletionLength < 50
                ? GATKSVVCFConstants.ComplexVariantSubtype.dDUP
                : GATKSVVCFConstants.ComplexVariantSubtype.dDUP_iDEL;

        final List<String> rewrittenIntervals = new ArrayList<>();
        rewrittenIntervals.add(new SVCallRecord.ComplexEventInterval(
                GATKSVVCFConstants.StructuralVariantAnnotationType.DUP,
                inversionSegment.getInterval()).encode());
        rewrittenIntervals.add(new SVCallRecord.ComplexEventInterval(
                GATKSVVCFConstants.StructuralVariantAnnotationType.INV,
                inversionSegment.getInterval()).encode());
        if (rewrittenSubtype == GATKSVVCFConstants.ComplexVariantSubtype.dDUP_iDEL) {
            rewrittenIntervals.add(new SVCallRecord.ComplexEventInterval(
                    GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
                    new SimpleInterval(variant)).encode());
        }

        attributes.put(GATKSVVCFConstants.SVTYPE, GATKSVVCFConstants.StructuralVariantAnnotationType.CPX);
        attributes.put(GATKSVVCFConstants.CPX_TYPE, SVCallRecordUtils.getComplexSubtypeString(rewrittenSubtype));
        attributes.put(GATKSVVCFConstants.CPX_INTERVALS, rewrittenIntervals);
        attributes.remove(SOURCE_ATTRIBUTE);
        return Allele.create("<CPX>", false);
    }

    /**
     * Builds the per-record evidence plan to re-evaluate carrier genotypes
     */
    static EvaluationPlan createEvaluationPlan(final VariantContext variant, final SAMSequenceDictionary dictionary) {
        final GATKSVVCFConstants.StructuralVariantAnnotationType svType =
                SVCallRecordUtils.inferStructuralVariantType(variant);
        final GATKSVVCFConstants.ComplexVariantSubtype complexSubtype = SVCallRecordUtils.getComplexSubtype(variant);

        if (svType == GATKSVVCFConstants.StructuralVariantAnnotationType.CPX) {
            return createComplexEvaluationPlan(variant, complexSubtype, dictionary);
        }
        if (svType == GATKSVVCFConstants.StructuralVariantAnnotationType.CTX) {
            return createTranslocationEvaluationPlan(variant, complexSubtype, dictionary);
        }
        if (svType == GATKSVVCFConstants.StructuralVariantAnnotationType.INS) {
            final String source = variant.getAttributeAsString(SOURCE_ATTRIBUTE, null);
            if (source != null && source.contains("INV")) {
                return createInsertionWithInversionEvaluationPlan(variant, dictionary, SVAnnotateEngine.parseComplexIntervals(Arrays.asList(source)));
            }
        }
        return EvaluationPlan.noEvaluation();
    }

    private static EvaluationPlan createComplexEvaluationPlan(final VariantContext variant,
                                                              final GATKSVVCFConstants.ComplexVariantSubtype complexSubtype,
                                                              final SAMSequenceDictionary dictionary) {
        if (complexSubtype == null) {
            throw new UserException.BadInput("Complex variant " + variant.getID() + " is missing CPX_TYPE");
        }
        if (STRUCTURALLY_UNRESOLVED_SUBTYPES.contains(complexSubtype)) {
            return EvaluationPlan.structuralUnresolved(EvaluatedVariantType.CPX);
        }

        final List<SVSegment> complexIntervals = SVAnnotateEngine.parseComplexIntervals(
                variant.getAttributeAsStringList(GATKSVVCFConstants.CPX_INTERVALS, null));
        switch (complexSubtype) {
            case dDUP:
            case dDUP_iDEL:
                return createDispersedDuplicationPlan(variant, complexIntervals, dictionary);
            case INS_iDEL:
                return createInsertionDeletionPlan(variant, complexIntervals, dictionary);
            case CTX_PP_QQ:
            case CTX_PQ_QP:
            case CTX_INV:
                return createTranslocationEvaluationPlan(variant, complexSubtype, dictionary);
            default:
                if (INVERSION_CNV_SUBTYPES.contains(complexSubtype)) {
                    return createInversionCnvPlan(variant, complexSubtype, complexIntervals);
                }
                return EvaluationPlan.structuralUnresolved(EvaluatedVariantType.CPX);
        }
    }

    private static EvaluationPlan createDispersedDuplicationPlan(final VariantContext variant,
                                                                 final List<SVSegment> complexIntervals,
                                                                 final SAMSequenceDictionary dictionary) {
        final SVSegment sourceDuplicationSegment = getRequiredSegment(
            complexIntervals,
            GATKSVVCFConstants.StructuralVariantAnnotationType.DUP,
            "dDUP complex variant is missing a DUP interval");
        final boolean sourceHasInversion =
            findFirstSegment(complexIntervals, GATKSVVCFConstants.StructuralVariantAnnotationType.INV) != null;
        final SVSegment sinkDeletionSegment =
            findFirstSegment(complexIntervals, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL);
        final SimpleInterval sinkInterval = sinkDeletionSegment == null
            ? new SimpleInterval(variant)
            : sinkDeletionSegment.getInterval();
        return createSourceSinkCpxEvaluationPlan(
            variant.getID(),
            sinkInterval,
            sourceDuplicationSegment.getInterval(),
            sourceHasInversion,
            dictionary);
    }

    private static EvaluationPlan createInsertionDeletionPlan(final VariantContext variant,
                                                              final List<SVSegment> complexIntervals,
                                                              final SAMSequenceDictionary dictionary) {
        final List<SVSegment> sourceSegments = SVAnnotateEngine.parseComplexIntervals(
            Arrays.asList(variant.getAttributeAsString(SOURCE_ATTRIBUTE, null)));
        final SVSegment sourceInsertionSegment = getRequiredSegment(
            sourceSegments,
            GATKSVVCFConstants.StructuralVariantAnnotationType.INS,
            "INS_iDEL variant is missing an INS source interval");
        final SVSegment sinkDeletionSegment = getRequiredSegment(
            complexIntervals,
            GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
            "INS_iDEL variant is missing a DEL interval in CPX_INTERVALS");
        final boolean sourceHasInversion =
            findFirstSegment(sourceSegments, GATKSVVCFConstants.StructuralVariantAnnotationType.INV) != null;
        return createSourceSinkCpxEvaluationPlan(
            variant.getID(),
            sinkDeletionSegment.getInterval(),
            sourceInsertionSegment.getInterval(),
            sourceHasInversion,
            dictionary);
    }

    private static EvaluationPlan createInsertionWithInversionEvaluationPlan(final VariantContext variant,
                                                                             final SAMSequenceDictionary dictionary,
                                                                             final List<SVSegment> sourceSegments) {
        final SVSegment sourceInversionSegment = getRequiredSegment(
            sourceSegments,
            GATKSVVCFConstants.StructuralVariantAnnotationType.INV,
            "Insertion-with-inversion record is missing an INV source interval");
        return createSourceSinkCpxEvaluationPlan(
            variant.getID(),
            new SimpleInterval(variant),
            sourceInversionSegment.getInterval(),
            true,
            dictionary);
    }

    private static EvaluationPlan createInversionCnvPlan(final VariantContext variant,
                                                         final GATKSVVCFConstants.ComplexVariantSubtype complexSubtype,
                                                         final List<SVSegment> complexIntervals) {
        final int[] breakpoints = getInversionCnvBreakpoints(complexIntervals, complexSubtype);
        final String contig = variant.getContig();
        final List<DiscordantPairQuery> queries = new ArrayList<>(2);

        switch (complexSubtype) {
            case delINV:
                queries.add(createBreakpointQuery(contig, breakpoints[0], true, contig, breakpoints[2], true));
                queries.add(createBreakpointQuery(contig, breakpoints[1], false, contig, breakpoints[2], false));
                break;
            case INVdel:
                queries.add(createBreakpointQuery(contig, breakpoints[0], true, contig, breakpoints[1], true));
                queries.add(createBreakpointQuery(contig, breakpoints[0], false, contig, breakpoints[2], false));
                break;
            case dupINV:
                queries.add(createBreakpointQuery(contig, breakpoints[1], true, contig, breakpoints[2], true));
                queries.add(createBreakpointQuery(contig, breakpoints[0], false, contig, breakpoints[2], false));
                break;
            case INVdup:
                queries.add(createBreakpointQuery(contig, breakpoints[0], true, contig, breakpoints[2], true));
                queries.add(createBreakpointQuery(contig, breakpoints[0], false, contig, breakpoints[1], false));
                break;
            case delINVdel:
                queries.add(createBreakpointQuery(contig, breakpoints[0], true, contig, breakpoints[2], true));
                queries.add(createBreakpointQuery(contig, breakpoints[1], false, contig, breakpoints[3], false));
                break;
            case delINVdup:
                queries.add(createBreakpointQuery(contig, breakpoints[0], true, contig, breakpoints[2], true));
                queries.add(createBreakpointQuery(contig, breakpoints[1], false, contig, breakpoints[4], false));
                break;
            case dupINVdel:
                queries.add(createBreakpointQuery(contig, breakpoints[1], true, contig, breakpoints[2], true));
                queries.add(createBreakpointQuery(contig, breakpoints[0], false, contig, breakpoints[3], false));
                break;
            case dupINVdup:
                queries.add(createBreakpointQuery(contig, breakpoints[1], true, contig, breakpoints[3], true));
                queries.add(createBreakpointQuery(contig, breakpoints[0], false, contig, breakpoints[4], false));
                break;
            default:
                return EvaluationPlan.structuralUnresolved(EvaluatedVariantType.CPX);
        }

        return new EvaluationPlan(EvaluatedVariantType.CPX, queries, true, false);
    }

    private static EvaluationPlan createTranslocationEvaluationPlan(final VariantContext variant,
                                                                    final GATKSVVCFConstants.ComplexVariantSubtype complexSubtype,
                                                                    final SAMSequenceDictionary dictionary) {
        if (complexSubtype != GATKSVVCFConstants.ComplexVariantSubtype.CTX_PQ_QP
                && complexSubtype != GATKSVVCFConstants.ComplexVariantSubtype.CTX_PP_QQ) {
            return EvaluationPlan.structuralUnresolved(EvaluatedVariantType.CTX);
        }

        final SimpleInterval first = new SimpleInterval(variant.getContig(), variant.getStart(), variant.getEnd());
        final String secondContig = variant.getAttributeAsString(GATKSVVCFConstants.CONTIG2_ATTRIBUTE, null);
        final int secondPosition = variant.getAttributeAsInt(GATKSVVCFConstants.END2_ATTRIBUTE, 0);
        if (secondContig == null || secondPosition <= 0) {
            throw new UserException.BadInput(
                "CTX variant " + variant.getID() + " is missing required CHR2/END2 attributes");
        }
        final SimpleInterval second = new SimpleInterval(secondContig, secondPosition, secondPosition);
        final boolean firstComesFirst = compareContigs(first.getContig(), second.getContig(), dictionary) <= 0;

        final SimpleInterval left = firstComesFirst ? first : second;
        final SimpleInterval right = firstComesFirst ? second : first;
        final List<DiscordantPairQuery> queries = new ArrayList<>(2);
        if (complexSubtype == GATKSVVCFConstants.ComplexVariantSubtype.CTX_PQ_QP) {
            queries.add(createBreakpointQuery(left.getContig(), left.getStart(), true, right.getContig(), right.getStart(), true));
            queries.add(createBreakpointQuery(left.getContig(), left.getEnd(), false, right.getContig(), right.getEnd(), false));
        } else {
            queries.add(createBreakpointQuery(left.getContig(), left.getStart(), true, right.getContig(), right.getStart(), false));
            queries.add(createBreakpointQuery(left.getContig(), left.getEnd(), false, right.getContig(), right.getEnd(), true));
        }
        return new EvaluationPlan(EvaluatedVariantType.CTX, queries, false, false);
    }

    private static EvaluationPlan createSourceSinkCpxEvaluationPlan(final String variantId,
                                                                    final SimpleInterval sinkInterval,
                                                                    final SimpleInterval sourceInterval,
                                                                    final boolean sourceHasInversion,
                                                                    final SAMSequenceDictionary dictionary) {
        final SinkPositionRelativeToSource sinkPosition =
            getSinkPositionRelativeToSource(sinkInterval, sourceInterval, dictionary);
        final List<DiscordantPairQuery> queries = new ArrayList<>(2);

        switch (sinkPosition) {
            case BEFORE_SOURCE:
                addQueriesForSinkBeforeSource(queries, sinkInterval, sourceInterval, sourceHasInversion);
                break;
            case AFTER_SOURCE:
                addQueriesForSinkAfterSource(queries, sinkInterval, sourceInterval, sourceHasInversion);
                break;
            case WITHIN_SOURCE:
                addQueriesForSinkWithinSource(queries, sinkInterval, sourceInterval, sourceHasInversion);
                break;
            default:
                throw new IllegalStateException("Unhandled sink/source ordering: " + sinkPosition);
        }

        return new EvaluationPlan(EvaluatedVariantType.CPX, queries, true, false);
    }

    private static void addQueriesForSinkBeforeSource(final List<DiscordantPairQuery> queries,
                                                      final SimpleInterval sinkInterval,
                                                      final SimpleInterval sourceInterval,
                                                      final boolean sourceHasInversion) {
        if (sourceHasInversion) {
            queries.add(createBreakpointQuery(
                    sinkInterval.getContig(), sinkInterval.getStart(), true,
                    sourceInterval.getContig(), sourceInterval.getEnd(), true));
            queries.add(createBreakpointQuery(
                    sinkInterval.getContig(), sinkInterval.getEnd(), false,
                    sourceInterval.getContig(), sourceInterval.getStart(), false));
            return;
        }

        queries.add(createBreakpointQuery(
                sinkInterval.getContig(), sinkInterval.getStart(), true,
                sourceInterval.getContig(), sourceInterval.getStart(), false));
        queries.add(createBreakpointQuery(
                sinkInterval.getContig(), sinkInterval.getEnd(), false,
                sourceInterval.getContig(), sourceInterval.getEnd(), true));
    }

    private static void addQueriesForSinkAfterSource(final List<DiscordantPairQuery> queries,
                                                     final SimpleInterval sinkInterval,
                                                     final SimpleInterval sourceInterval,
                                                     final boolean sourceHasInversion) {
        if (sourceHasInversion) {
            queries.add(createBreakpointQuery(
                    sourceInterval.getContig(), sourceInterval.getEnd(), true,
                    sinkInterval.getContig(), sinkInterval.getStart(), true));
            queries.add(createBreakpointQuery(
                    sourceInterval.getContig(), sourceInterval.getStart(), false,
                    sinkInterval.getContig(), sinkInterval.getEnd(), false));
            return;
        }

        queries.add(createBreakpointQuery(
                sourceInterval.getContig(), sourceInterval.getEnd(), true,
                sinkInterval.getContig(), sinkInterval.getEnd(), false));
        queries.add(createBreakpointQuery(
                sourceInterval.getContig(), sourceInterval.getStart(), false,
                sinkInterval.getContig(), sinkInterval.getStart(), true));
    }

    private static void addQueriesForSinkWithinSource(final List<DiscordantPairQuery> queries,
                                                      final SimpleInterval sinkInterval,
                                                      final SimpleInterval sourceInterval,
                                                      final boolean sourceHasInversion) {
        final String contig = sinkInterval.getContig();
        if (sourceHasInversion) {
            queries.add(createBreakpointQuery(contig, sinkInterval.getStart(), false, contig, sourceInterval.getEnd(), false));
            queries.add(createBreakpointQuery(contig, sourceInterval.getStart(), true, contig, sinkInterval.getEnd(), true));
            return;
        }
        queries.add(createBreakpointQuery(contig, sourceInterval.getStart(), false, contig, sinkInterval.getStart(), true));
        queries.add(createBreakpointQuery(contig, sinkInterval.getEnd(), false, contig, sourceInterval.getEnd(), true));
    }

    private static DiscordantPairQuery createBreakpointQuery(final String startContig,
                                                             final int startBreakpoint,
                                                             final boolean startStrand,
                                                             final String endContig,
                                                             final int endBreakpoint,
                                                             final boolean endStrand) {
        final int startMin = startStrand ? startBreakpoint - PE_FLANK_BACK : startBreakpoint - PE_FLANK_FRONT;
        final int startMax = startStrand ? startBreakpoint + PE_FLANK_FRONT : startBreakpoint + PE_FLANK_BACK;
        final int endLowerExclusive = endStrand ? endBreakpoint - PE_FLANK_BACK : endBreakpoint - PE_FLANK_FRONT;
        final int endUpperExclusive = endStrand ? endBreakpoint + PE_FLANK_FRONT : endBreakpoint + PE_FLANK_BACK;
        return new DiscordantPairQuery(
                startContig,
                startMin,
                startMax,
                startStrand,
                endContig,
                endLowerExclusive,
                endUpperExclusive,
                endStrand);
    }

    private static int[] getInversionCnvBreakpoints(final List<SVSegment> complexIntervals,
                                                    final GATKSVVCFConstants.ComplexVariantSubtype complexSubtype) {
        final List<Integer> breakpoints = new ArrayList<>();
        final SVSegment firstSegment = complexIntervals.get(0);
        breakpoints.add(firstSegment.getStart());
        breakpoints.add(firstSegment.getEnd());

        if (complexSubtype != GATKSVVCFConstants.ComplexVariantSubtype.INVdup) {
            for (int i = 1; i < complexIntervals.size(); i++) {
                breakpoints.add(complexIntervals.get(i).getEnd());
            }
        } else {
            final SVSegment secondSegment = complexIntervals.get(1);
            if (secondSegment.getStart() < breakpoints.get(1)) {
                breakpoints.set(1, secondSegment.getStart());
                breakpoints.add(firstSegment.getEnd());
            }
        }

        if (complexSubtype == GATKSVVCFConstants.ComplexVariantSubtype.dupINVdup
                || complexSubtype == GATKSVVCFConstants.ComplexVariantSubtype.delINVdup) {
            final SVSegment thirdSegment = complexIntervals.get(2);
            breakpoints.add(thirdSegment.getStart());
            breakpoints.add(thirdSegment.getEnd());
        }

        return breakpoints.stream().mapToInt(Integer::intValue).toArray();
    }

    private static boolean requiresFormattingTransform(final VariantContext variant) {
        final GATKSVVCFConstants.StructuralVariantAnnotationType svType =
                SVCallRecordUtils.inferStructuralVariantType(variant);
        final GATKSVVCFConstants.ComplexVariantSubtype complexSubtype = SVCallRecordUtils.getComplexSubtype(variant);
        if (svType == GATKSVVCFConstants.StructuralVariantAnnotationType.CPX
                && (complexSubtype == GATKSVVCFConstants.ComplexVariantSubtype.dDUP
                || complexSubtype == GATKSVVCFConstants.ComplexVariantSubtype.dDUP_iDEL)) {
            return variant.hasAttribute(SOURCE_ATTRIBUTE);
        }
        if (svType == GATKSVVCFConstants.StructuralVariantAnnotationType.INS) {
            final String source = variant.getAttributeAsString(SOURCE_ATTRIBUTE, null);
            return source != null && source.contains("INV");
        }
        return false;
    }

    private static boolean shouldReevaluatePreexistingUnresolved(final VariantContext variant) {
        if (SVCallRecordUtils.inferStructuralVariantType(variant)
                != GATKSVVCFConstants.StructuralVariantAnnotationType.INS) {
            return false;
        }
        final String source = variant.getAttributeAsString(SOURCE_ATTRIBUTE, null);
        return source != null && source.contains("INV");
    }

    private static GenotypesContext remapGenotypesToAlternateAllele(final GenotypesContext genotypes,
                                                                    final Allele alternateAllele) {
        final ArrayList<Genotype> updated = new ArrayList<>(genotypes.size());
        for (final Genotype genotype : genotypes) {
            final List<Allele> newAlleles = genotype.getAlleles().stream()
                    .map(allele -> allele == null || allele.isReference() || allele.isNoCall() ? allele : alternateAllele)
                    .collect(Collectors.toList());
            updated.add(new GenotypeBuilder(genotype).alleles(newAlleles).make());
        }
        return GenotypesContext.create(updated);
    }

    private static Genotype makeNoCallGenotype(final Genotype genotype) {
        final int ploidy = genotype.getPloidy() > 0 ? genotype.getPloidy() : 2;
        return new GenotypeBuilder(genotype)
                .alleles(Collections.nCopies(ploidy, Allele.NO_CALL))
                .make();
    }

    private static List<String> getCarrierSamples(final VariantContext variant) {
        return variant.getGenotypes().stream()
                .filter(SVCallRecordUtils::isAltGenotype)
                .map(Genotype::getSampleName)
                .collect(Collectors.toList());
    }

    private static int compareContigs(final String left,
                                      final String right,
                                      final SAMSequenceDictionary dictionary) {
        return Integer.compare(dictionary.getSequenceIndex(left), dictionary.getSequenceIndex(right));
    }

    static SinkPositionRelativeToSource getSinkPositionRelativeToSource(final SimpleInterval sinkInterval,
                                                                        final SimpleInterval sourceInterval,
                                                                        final SAMSequenceDictionary dictionary) {
        final int contigComparison = compareContigs(sinkInterval.getContig(), sourceInterval.getContig(), dictionary);
        if (contigComparison < 0) {
            return SinkPositionRelativeToSource.BEFORE_SOURCE;
        }
        if (contigComparison > 0) {
            return SinkPositionRelativeToSource.AFTER_SOURCE;
        }
        if (sinkInterval.getStart() <= sourceInterval.getStart()) {
            return SinkPositionRelativeToSource.BEFORE_SOURCE;
        }
        if (sinkInterval.getEnd() >= sourceInterval.getEnd()) {
            return SinkPositionRelativeToSource.AFTER_SOURCE;
        }
        return SinkPositionRelativeToSource.WITHIN_SOURCE;
    }


    private static SVSegment findFirstSegment(final List<SVSegment> segments,
                                              final GATKSVVCFConstants.StructuralVariantAnnotationType type) {
        return segments.stream()
                .filter(segment -> segment.getIntervalSVType() == type)
                .findFirst()
                .orElse(null);
    }

    private static SVSegment getRequiredSegment(final List<SVSegment> segments,
                                                final GATKSVVCFConstants.StructuralVariantAnnotationType type,
                                                final String errorMessage) {
        final SVSegment segment = findFirstSegment(segments, type);
        if (segment == null) {
            throw new UserException.BadInput(errorMessage);
        }
        return segment;
    }

    enum EvaluatedVariantType {
        CPX,
        CTX,
        NONE
    }

    enum SinkPositionRelativeToSource {
        BEFORE_SOURCE,
        AFTER_SOURCE,
        WITHIN_SOURCE
    }

    static final class EvaluationPlan {
        final EvaluatedVariantType variantType;
        final List<DiscordantPairQuery> queries;
        final boolean requiresDepthEvidence;
        final boolean structuralUnresolved;

        EvaluationPlan(final EvaluatedVariantType variantType,
                       final List<DiscordantPairQuery> queries,
                       final boolean requiresDepthEvidence,
                       final boolean structuralUnresolved) {
            this.variantType = variantType;
            this.queries = queries;
            this.requiresDepthEvidence = requiresDepthEvidence;
            this.structuralUnresolved = structuralUnresolved;
        }

        static EvaluationPlan structuralUnresolved(final EvaluatedVariantType variantType) {
            return new EvaluationPlan(variantType, Collections.emptyList(), false, true);
        }

        static EvaluationPlan noEvaluation() {
            return new EvaluationPlan(EvaluatedVariantType.NONE, Collections.emptyList(), false, false);
        }
    }

    static final class DiscordantPairQuery {
        final String startContig;
        final int startMin;
        final int startMax;
        final boolean startStrand;
        final String endContig;
        final int endLowerExclusive;
        final int endUpperExclusive;
        final boolean endStrand;

        DiscordantPairQuery(final String startContig,
                            final int startMin,
                            final int startMax,
                            final boolean startStrand,
                            final String endContig,
                            final int endLowerExclusive,
                            final int endUpperExclusive,
                            final boolean endStrand) {
            this.startContig = startContig;
            this.startMin = startMin;
            this.startMax = startMax;
            this.startStrand = startStrand;
            this.endContig = endContig;
            this.endLowerExclusive = endLowerExclusive;
            this.endUpperExclusive = endUpperExclusive;
            this.endStrand = endStrand;
        }
    }

    static final class CarrierRefinement {
        static final CarrierRefinement KEEP = new CarrierRefinement(false, false);

        final boolean reviseGenotype;
        final boolean countTowardsUnresolved;

        CarrierRefinement(final boolean reviseGenotype, final boolean countTowardsUnresolved) {
            this.reviseGenotype = reviseGenotype;
            this.countTowardsUnresolved = countTowardsUnresolved;
        }
    }

    static double computeCoverageFraction(final List<DepthInterval> intervals,
                                          final int queryStart,
                                          final int queryEnd) {
        if (queryEnd <= queryStart || intervals.isEmpty()) {
            return 0.0;
        }

        int covered = 0;
        int currentStart = -1;
        int currentEnd = -1;
        for (final DepthInterval interval : intervals) {
            if (interval.end <= queryStart) {
                continue;
            }
            if (interval.start >= queryEnd) {
                break;
            }

            final int overlapStart = Math.max(queryStart, interval.start);
            final int overlapEnd = Math.min(queryEnd, interval.end);
            if (overlapStart >= overlapEnd) {
                continue;
            }
            if (currentStart < 0) {
                currentStart = overlapStart;
                currentEnd = overlapEnd;
            } else if (overlapStart <= currentEnd) {
                currentEnd = Math.max(currentEnd, overlapEnd);
            } else {
                covered += currentEnd - currentStart;
                currentStart = overlapStart;
                currentEnd = overlapEnd;
            }
        }

        if (currentStart >= 0) {
            covered += currentEnd - currentStart;
        }
        return covered / (double) (queryEnd - queryStart);
    }

    static final class DepthInterval {
        final int start;
        final int end;

        DepthInterval(final int start, final int end) {
            this.start = start;
            this.end = end;
        }
    }
}
