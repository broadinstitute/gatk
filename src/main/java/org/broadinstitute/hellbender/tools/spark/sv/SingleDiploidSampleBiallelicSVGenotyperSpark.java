package org.broadinstitute.hellbender.tools.spark.sv;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Sets;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
import org.apache.commons.io.FilenameUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.cmdline.Advanced;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.engine.datasources.ReferenceWindowFunctions;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.engine.spark.datasources.VariantsSparkSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeAssignmentMethod;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeLikelihoodCalculator;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeLikelihoodCalculators;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import scala.Tuple2;

import java.io.BufferedOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.*;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

/**
 * TODO: when we move to genotyping multiple types of sv, this should be moved to an engine class instead of base tool class.
 * TODO: right now we are dealing with single sample only, so we don't need a {@link org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypingData} alike,
 *       for the purpose of storing a {@link org.broadinstitute.hellbender.tools.walkers.genotyper.PloidyModel}
 *       and getting the appropriate {@link org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix} for the sample.
 *       it may make sense to refactor {@link ReadLikelihoods} for a more explicit sample-stratified structure
 */
abstract class SingleDiploidSampleBiallelicSVGenotyperSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    @Argument(doc       = "An URI to the directory where all called inversion breakends are stored as text files. " +
            "The stored file format should conform to the one specified in {@link ContigAligner.AssembledBreakpoint#toString()}",
            shortName = StandardArgumentDefinitions.VARIANT_SHORT_NAME,
            fullName  = StandardArgumentDefinitions.VARIANT_LONG_NAME,
            optional  = false)
    protected String VCFFile;

    @Argument(doc       = "An URI to the directory where all called inversion breakends are stored as text files. " +
            "The stored file format should conform to the one specified in {@link ContigAligner.AssembledBreakpoint#toString()}",
            shortName = "fastq",
            fullName  = "FASTQDir",
            optional  = false)
    protected String pathToFASTQFiles;

    @Argument(doc       = "An URI to the directory where all called inversion breakends are stored as text files. " +
            "The stored file format should conform to the one specified in {@link ContigAligner.AssembledBreakpoint#toString()}",
            shortName = "assembly",
            fullName  = "assembledFASTADir",
            optional  = false)
    protected String pathToAssembledFASTAFiles;

    @Argument(doc       = "An URI to the directory where all called inversion breakends are stored as text files. " +
            "The stored file format should conform to the one specified in {@link ContigAligner.AssembledBreakpoint#toString()}",
            shortName = "aln",
            fullName  = "contigAlignments",
            optional  = false)
    protected String pathToContigAlignments;

    @Argument(doc       = "An URI to the directory where all called inversion breakends are stored as text files. " +
            "The stored file format should conform to the one specified in {@link ContigAligner.AssembledBreakpoint#toString()}",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            optional = false)
    protected String outFileName;


    @Argument(doc = "FASTA formatted reference", shortName = "fastaReference",
            fullName = "fastaReference", optional = false)
    private String fastaReference;

    @Argument(doc = "Path to a file where debugging string will be writen out to", shortName = "dso",
            fullName = "debugstringoutpath", optional = true)
    @Advanced
    private String debugStringSavingOutputPath;


    private static final int ploidy = 2;
    @Override
    public final boolean requiresReference() {
        return true;
    }

    @Override
    public final boolean requiresReads() {
        return true;
    }

    // for developer performance debugging use
    static final boolean in_debug_state = true;
    private static final Logger logger = LogManager.getLogger(SingleDiploidSampleBiallelicSVGenotyperSpark.class);

    static final String testSampleName = "NA12878_PCR-_30X";

    // -----------------------------------------------------------------------------------------------
    // Manager
    // -----------------------------------------------------------------------------------------------

    @Override
    public void runTool(final JavaSparkContext ctx){

        // parse discovery VC calls
        final List<SVJunction> junctions = parseVCFInput(ctx, VCFFile, pathToAssembledFASTAFiles, pathToContigAlignments, getReference());
        filterJunctions(junctions, makeJunctionFilter());

        // gather and filter reads that will be used in genotyping
        final JavaPairRDD<SVJunction, List<GATKRead>> genotypingReads = loadGenotypingReads(ctx, junctions, pathToFASTQFiles).cache();
        final Set<SVJunction> sparkJunctions = genotypingReads.keys().collect().stream().collect(Collectors.toSet());
        if(Sets.symmetricDifference(sparkJunctions, junctions.stream().collect(Collectors.toSet())).size()!=0){
            throw new IllegalStateException("Junctions are different before and after read gathering!");
        }
        logger.info(".........DONE MERGING READS.........");

        final JavaPairRDD<SVJunction, ReadLikelihoods<SVDummyAllele>> readLikelihoodsJavaPairRDD = genotypingReads.mapToPair(this::computeReadLikelihoods);
        genotypingReads.unpersist();

        final JavaPairRDD<SVJunction, GenotypeLikelihoods> step3 = readLikelihoodsJavaPairRDD.mapValues(this::getGenotypeLikelihoods).cache();
        readLikelihoodsJavaPairRDD.unpersist();

        // turn back to VC and save
        final SAMSequenceDictionary referenceSequenceDictionary = new ReferenceMultiSource(getAuthenticatedGCSOptions(), fastaReference, ReferenceWindowFunctions.IDENTITY_FUNCTION).getReferenceSequenceDictionary(null);
        final List<VariantContext> genotypedVCs = collectAndSortGenotypedVC(step3, referenceSequenceDictionary);

        try (final OutputStream outputStream = new BufferedOutputStream(BucketUtils.createFile(outFileName, getAuthenticatedGCSOptions()))){
            output(genotypedVCs, referenceSequenceDictionary, outputStream);
            writeoutDebugReadLikelihoods(step3.coalesce(1).keys(), debugStringSavingOutputPath);
        } catch (final IOException ioex){
            throw new GATKException("Could not create output file", ioex);
        }
    }

    // -----------------------------------------------------------------------------------------------
    // To be overridden
    // -----------------------------------------------------------------------------------------------

    /**
     * Converting from a discovery {@link VariantContext} call to an appropriate {@link SVJunction} suitable for genotyping.
     */
    protected abstract SVJunction convertToSVJunction(final VariantContext vc,
                                                      final Broadcast<Map<Long, List<LocalAssemblyContig>>> assemblyID2assemblyContents,
                                                      final Broadcast<ReferenceMultiSource> referenceMultiSourceBroadcast);

    protected Predicate<SVJunction> makeJunctionFilter(){
        return new JunctionFilterBasedOnInterval();
    }

    /**
     * Go back to the original bam and gather reads around region identified by {@link SVJunction}.
     *
     * Reads for each junction are filtered based on subclass custom defined {@link #readSuitableForGenotypingJunction(SVJunction, GATKRead)}.
     *
     * Note that read names are changed by appending "/1" or "/2" to their original names so later when bam reads are
     * merged with FASTQ reads extracted by {@link #getFASTQReads(JavaPairRDD, JavaPairRDD)},
     * duplicate reads can be removed by read names.
     *
     * Note also that the interface seems redundant because we want to test this method (so calling getReads() internally would make that impossible).
     */
    @VisibleForTesting
    protected JavaPairRDD<SVJunction, List<GATKRead>> getBamReads(final JavaRDD<GATKRead> reads, final List<SVJunction> svJunctions){
        return reads
                .filter(read -> svJunctions.stream().anyMatch(junction -> readSuitableForGenotypingJunction(junction, read))) // pre-filer: if suitable for any junction
                .flatMapToPair( read -> svJunctions.stream()
                                                    .filter(junction -> readSuitableForGenotypingJunction(junction, read))
                                                    .map(junc -> new Tuple2<>(junc, read))
                                                    .collect(Collectors.toList()) )                                           // read -> { (junction, read) }
                .groupByKey()
                .mapValues( it -> {// bwa strips out the "/1" "/2" part, but we need these for merging with FASTQ reads
                            it.forEach(read -> read.setName(read.getName() + (read.isFirstOfPair() ? "/1" : "/2")));
                            return StreamSupport.stream(it.spliterator(), false).collect(Collectors.toList());}); // turn iterable to list
    }

    /**
     * Filter read that will be used for genotyping.
     * Assumes the logic applies to all junctions of the same type.
     * TODO: this interface implicitly works against pairedness.
     */
    protected abstract boolean readSuitableForGenotypingJunction(final SVJunction junction,
                                                                 final GATKRead read);

    /**
     * Update read based on result from its realignment to ref and alt contigs.
     * Input reads are deep copied for making the necessary changes.
     */
    protected abstract ReadLikelihoods<SVDummyAllele> updateReads(final SVJunction junction, final ReadLikelihoods<SVDummyAllele> reads);

    // -----------------------------------------------------------------------------------------------
    // Utilities: IO
    // -----------------------------------------------------------------------------------------------

    /**
     * Parse input (output from previous step in the pipeline),
     * and return a list of convenience struct for use by specific genotyping module.
     *
     * Depends on subclasses' {@link #convertToSVJunction(VariantContext, Broadcast, Broadcast)} for the actual mapping.
     *
     * Not tested because most all methods this depends on are tested else where already.
     */
    private List<SVJunction> parseVCFInput(final JavaSparkContext ctx,
                                           final String vcfInput,
                                           final String localAssemblyResultPath,
                                           final String alignmentResultPath,
                                           final ReferenceMultiSource reference){

        final List<VariantContext> listOfVCs = new VariantsSparkSource(ctx).getParallelVariantContexts(vcfInput, null).collect();

        final JavaPairRDD<Long, ContigsCollection> asmId2Contigs = ContigsCollection.loadContigsCollectionKeyedByAssemblyId(ctx, localAssemblyResultPath)
                .mapToPair( p -> new Tuple2<>(Long.parseLong(p._1().replace(SVConstants.FASTQ_OUT_PREFIX, "")), p._2()));

        final JavaPairRDD<Long, Iterable<AlignmentRegion>> asmId2Alignments = ctx.textFile(alignmentResultPath)
                .map(ContigsCollection::parseAlignedAssembledContigLine)
                .mapToPair(ar -> new Tuple2<>(Long.parseLong(ar.assemblyId.replace(SVConstants.FASTQ_OUT_PREFIX, "")), ar))
                .groupByKey();

        final Map<Long, List<LocalAssemblyContig>> assemblyID2assembleContents = mapAssemblyID2ItsAlignments(asmId2Contigs, asmId2Alignments, gatherAsmIDs(listOfVCs));
        final Broadcast<Map<Long, List<LocalAssemblyContig>>> assemblyID2assembleContentsBroadcast = ctx.broadcast(assemblyID2assembleContents);

        final Broadcast<ReferenceMultiSource> referenceMultiSourceBroadcast = ctx.broadcast(reference);

        if(in_debug_state) logger.info(".........BEGIN CONSTRUCTING JUNCTIONS FROM VC AND ALIGNMENTS.........");
        final List<SVJunction> junctions = listOfVCs.stream().map(vc -> convertToSVJunction(vc, assemblyID2assembleContentsBroadcast, referenceMultiSourceBroadcast)).collect(Collectors.toList());
        if(in_debug_state) logger.info(String.format(".........CONSTRUCTED %d JUNCTIONS FROM VC AND ALIGNMENTS.........", junctions.size()));
        assemblyID2assembleContentsBroadcast.destroy();referenceMultiSourceBroadcast.destroy();

        return junctions;
    }

    /**
     * @return a mapping from assembly id to the result of locally assembled contigs aligned back to ref,
     *         given the path to re-alignment result.
     */
    @VisibleForTesting
    static Map<Long, List<LocalAssemblyContig>> mapAssemblyID2ItsAlignments(final JavaPairRDD<Long, ContigsCollection> asmId2Contigs,
                                                                            final JavaPairRDD<Long, Iterable<AlignmentRegion>> asmId2Alignments,
                                                                            final Set<Long> assembliesToKeep){
        return asmId2Alignments.join(asmId2Contigs).mapValues( pair -> {
            final Map<String, ArrayList<AlignmentRegion>> contigId2Alignments = StreamSupport.stream(pair._1.spliterator(), false).collect(Collectors.groupingBy(ar -> ar.contigId))
                    .entrySet().stream().collect(Collectors.toMap( p -> p.getKey(), p -> new ArrayList<>(p.getValue())));
            final List<LocalAssemblyContig> collection = pair._2.getContents() ;
            collection.forEach(contig -> contig.setAlignmentRegions(contigId2Alignments.get(contig.contigID)));
            return collection;
        }).filter(p -> assembliesToKeep.contains(p._1)).collectAsMap();
    }

    @VisibleForTesting
    static Set<Long> gatherAsmIDs(final List<VariantContext> vcs){
        return vcs.stream()
                .flatMap(vc -> Arrays.stream(vc.getAttributeAsString(GATKSVVCFHeaderLines.ASSEMBLY_IDS, "").replace("[", "").replace("]", "").split(", ")).map(Long::valueOf)  )
                .collect(Collectors.toSet());
    }

    @VisibleForTesting
    final static class JunctionFilterBasedOnInterval implements Predicate<SVJunction>{
        private static final List<SimpleInterval> defaultRegionsToExclude = Arrays.asList(new SimpleInterval("2", 33141200, 33141700),
                new SimpleInterval("MT", 12867, 14977),
                new SimpleInterval("17", 41379304, 41401057),
                new SimpleInterval("17", 41381590, 41463575),
                new SimpleInterval("17", 41401193, 41466836));

        private final List<SimpleInterval> regionsToExclude;

        JunctionFilterBasedOnInterval(){ regionsToExclude = defaultRegionsToExclude; }

        JunctionFilterBasedOnInterval(final List<SimpleInterval> regionsToExclude) { this.regionsToExclude = regionsToExclude;}

        public boolean test(final SVJunction j){
            final Tuple2<SimpleInterval, SimpleInterval> windows = j.getReferenceWindows();
            final SimpleInterval left = windows._1();
            final SimpleInterval right = windows._2();
            return regionsToExclude.stream().anyMatch(region -> region.overlaps(left) || region.overlaps(right));
        }
    }

    @VisibleForTesting
    static void filterJunctions(final List<SVJunction> junctions, final Predicate<SVJunction> predicate){
        final List<SVJunction> junctionsToBeFiltered = junctions.stream().filter(predicate::test).collect(Collectors.toList());
        if(junctionsToBeFiltered.size()!=0){
            logger.warn(".........FOLLOWING JUNCTIONS TO BE EXCLUDED IN GENOTYPING BECAUSE OF HIGH COVERAGE.........");
            junctionsToBeFiltered.stream().map(j -> j.getOriginalVC().getID()).forEach(logger::warn);
            junctions.removeAll(junctionsToBeFiltered);
            logger.warn(String.format(".........%d JUNCTIONS LEFT TO BE GENOTYPED.........", junctions.size()));
        }
    }

    /**
     * Prepare reads that will be used for likelihood calculations.
     *
     * Reads are gathered from two resources:
     * 1) original FASTQ records pulled out by {@link FindBreakpointEvidenceSpark} that are associated with
     *      putative breakpoints, which in turn the appropriate {@link SVJunction} are associated with.
     * 2) original bam file. Exactly how reads are pulled in from the bam file is determined by {@link #getBamReads(JavaRDD, List)},
     *      whose exact behavior is overridable by subclasses.
     *
     * For a junction, reads that are gathered from both sources (bam, and reconstructed from FASTQ), will be de-duplicated
     * such that the bam version will be kept.
     * @param svJunctions           a list of sv junctions to be genotyped
     * @param pathToFASTQFiles      a string to the directory where the FASTQ files generated by {@link FindBreakpointEvidenceSpark} are stored
     */
    private JavaPairRDD<SVJunction, List<GATKRead>> loadGenotypingReads(final JavaSparkContext ctx,
                                                                        final List<SVJunction> svJunctions,
                                                                        final String pathToFASTQFiles){

        // gather reads from FASTQ files collected in previous steps in pipeline
        final JavaPairRDD<Long, String> fastq = SVFastqUtils.loadFASTQFiles(ctx, pathToFASTQFiles)
                .mapToPair(pair -> new Tuple2<>(Long.parseLong(FilenameUtils.getBaseName(pair._1()).replace(SVConstants.FASTQ_OUT_PREFIX, "")), pair._2())); // essentially extract asm id from FASTQ file name

        final JavaPairRDD<Long, List<SVJunction>> asmid2SVJunctions = ctx.parallelizePairs(getAssemblyToJunctionsMapping(svJunctions).entrySet().stream().map(e -> new Tuple2<>(e.getKey(), e.getValue())).collect(Collectors.toList()));

        final JavaPairRDD<SVJunction, List<GATKRead>> fastqReads = getFASTQReads(fastq, asmid2SVJunctions);

        // gather reads from bam file (will overlap with FASTQ reads)
        if(in_debug_state) logger.info(".........BEGIN EXTRACTING BAM READS.........");
        final JavaPairRDD<SVJunction, List<GATKRead>> bamReads = getBamReads(getReads(), svJunctions).cache();
        if(in_debug_state){
            bamReads.count(); logger.info(".........DONE EXTRACTING BAM READS.........");
        }

        // merge together the fastqReads and bamReads (deduplication to follow)
        if(in_debug_state) logger.info(".........BEGIN MERGING READS.........");
        final JavaPairRDD<SVJunction, Tuple2<List<GATKRead>, List<GATKRead>>> redundantReads = fastqReads.join(bamReads).cache();
        bamReads.unpersist();

        return redundantReads.mapValues(pair -> mergeReadsByName(pair._1, pair._2));
    }

    /**
     * Get FASTQ reads that were pulled out by {@link FindBreakpointEvidenceSpark}.
     * Reads are filtered by custom filters defined by sub classes (see {@link #readSuitableForGenotypingJunction(SVJunction, GATKRead) <ReferenceMultiSource>)}.
     * @return Reads associated with a particular breakpoint.
     */
    @VisibleForTesting
    JavaPairRDD<SVJunction, List<GATKRead>> getFASTQReads(final JavaPairRDD<Long, String> fastqContents,
                                                          final JavaPairRDD<Long, List<SVJunction>> asmid2SVJunctions){

        final Set<Long> interestingAsmIDs = new HashSet<>( asmid2SVJunctions.keys().collect() );

        // map assembly id (only those giving signal) to FASTQ contents then reconstruct read
        final JavaPairRDD<Long, List<GATKRead>> preFilteredReads = fastqContents.filter(pair -> interestingAsmIDs.contains(pair._1())) // filter out unused assemblies
                                                                                .mapValues(contents -> SVFastqUtils.convertToReads(SVFastqUtils.extractFASTQContents(contents))); // FASTQ-> GATKRead

        return preFilteredReads.join(asmid2SVJunctions)
                .mapValues(pair -> {
                    final List<GATKRead> reads = pair._1();
                    return pair._2.stream().map(j -> new Tuple2<>(j, reads.stream().filter(read -> readSuitableForGenotypingJunction(j, read)).collect(Collectors.toList())))
                                                                        .collect(Collectors.toList());
                }) // (long, {(junction, {reads})}
                .values()
                .flatMapToPair(ls -> ls) // (junction, {reads}) with duplicated keys
                .groupByKey()
                .mapValues(it -> StreamSupport.stream(it.spliterator(), false).flatMap(List::stream).distinct().collect(Collectors.toList()) );
    }

    @VisibleForTesting
    static Map<Long, List<SVJunction>> getAssemblyToJunctionsMapping(final List<SVJunction> svJunctions){
        // get mapping from assembly id to SVJunction's; see discussion on StackOverflow 38471056
        return svJunctions.stream()
                    .flatMap(svJunction -> svJunction.getAssemblyIDs().stream().map(asid -> new AbstractMap.SimpleEntry<>(asid, svJunction)))
                    .collect(Collectors.groupingBy(Map.Entry::getKey, Collectors.mapping(Map.Entry::getValue, Collectors.toList())));
    }

    @VisibleForTesting
    static List<GATKRead> mergeReadsByName(final List<GATKRead> fastqReads,
                                           final List<GATKRead> bamReads){
        // for reads appear in both, keep the one from bam
        final Set<String> common = Sets.intersection(fastqReads.stream().map(GATKRead::getName).collect(Collectors.toSet()),
                                                     bamReads.stream().map(GATKRead::getName).collect(Collectors.toSet()));
        final List<GATKRead> toReturn = new ArrayList<>(bamReads);
        toReturn.addAll( fastqReads.stream().filter(r -> !common.contains(r.getName())).collect(Collectors.toList()) );
        return toReturn;
    }

    @VisibleForTesting
    static List<VariantContext> collectAndSortGenotypedVC(final JavaPairRDD<SVJunction, GenotypeLikelihoods> genotypedJunctions,
                                                          final SAMSequenceDictionary referenceSequenceDictionary){
        return StreamSupport.stream(genotypedJunctions.map(pair -> pair._1().setPL(pair._2().getAsPLs())).collect().spliterator(), false)// set PL
                .map(SVJunction::getGenotypedVC) // add annotation and update VC
                .sorted((VariantContext v1, VariantContext v2) -> IntervalUtils.compareLocatables(v1, v2, referenceSequenceDictionary))
                .collect(Collectors.toList());
    }

    /**
     * TODO: should we emit reference confidence, we need more
     * Writes out new VCF FORMAT fields computed by this tool.
     */
    @VisibleForTesting
    static void output(final List<VariantContext> vcList,
                       final SAMSequenceDictionary referenceSequenceDictionary,
                       final OutputStream outputStream){

        try(final VariantContextWriter vcfWriter = new VariantContextWriterBuilder().clearOptions()
                .setOutputStream(outputStream)
                .setReferenceDictionary(referenceSequenceDictionary)
                .setOption(Options.ALLOW_MISSING_FIELDS_IN_HEADER)
                .setOptions(Arrays.stream(new Options[]{}).collect(Collectors.toCollection(()->EnumSet.noneOf(Options.class))))
                .setOption(Options.WRITE_FULL_FORMAT_FIELD)
                .build()){
            final VCFHeader header = new VCFHeader(GATKSVVCFHeaderLines.vcfHeaderLines.values().stream().collect(Collectors.toSet()), Collections.singletonList(testSampleName));
            header.setSequenceDictionary(referenceSequenceDictionary);
            header.addMetaDataLine( VCFStandardHeaderLines.getInfoLine(VCFConstants.END_KEY) );

            header.addMetaDataLine( VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_KEY) );
            header.addMetaDataLine( VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_QUALITY_KEY));
            header.addMetaDataLine( VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_PL_KEY) );

            vcfWriter.writeHeader(header);
            vcList.forEach(vc -> logger.debug(String.format(".........%s:%d-%d\t%s.........", vc.getContig(), vc.getStart(), vc.getEnd(), vc.getGenotype(0).getLikelihoodsString())));
            vcList.forEach(vcfWriter::add);
        }
    }

    // -----------------------------------------------------------------------------------------------
    // Utilities: delegation (to avoid lambda calls that requires too many other classes to be serialized)
    // -----------------------------------------------------------------------------------------------

    private Tuple2<SVJunction, ReadLikelihoods<SVDummyAllele>> computeReadLikelihoods(final Tuple2<SVJunction, List<GATKRead>> readsForAJunction){
        // send reads to read likelihood calculator, but first configure the calculator appropriately
        // TODO: think about how to configure calculator
        SVReadLikelihoodCalculator readLikelihoodCalculator = new InversionReadLikelihoodCalculator();
        readLikelihoodCalculator.initialize();
        readLikelihoodCalculator.configure();
        final List<GATKRead> preprocessedReads = readLikelihoodCalculator.preprocessReads(readsForAJunction._2());
        logger.debug(".........DONE PREPROCESSING READS.........");

        final SVJunction junction = readsForAJunction._1();
        final Map<String, List<GATKRead>> sample2Reads = new LinkedHashMap<>();

        final SampleList sampleList = SampleList.singletonSampleList(testSampleName);
        sample2Reads.put(sampleList.getSample(0), preprocessedReads);

        final ReadLikelihoods<SVDummyAllele> rll = readLikelihoodCalculator.computeReadLikelihoods(sampleList, preprocessedReads, sample2Reads, junction);
        readLikelihoodCalculator.close();
        logger.debug(".........RLC DONE.........");

        final ReadLikelihoods<SVDummyAllele> result = updateReads(junction, rll);
        logger.debug(".........POSTPROCESSING DONE.........");

        return new Tuple2<>(junction, result);
    }

    private GenotypeLikelihoods getGenotypeLikelihoods(final ReadLikelihoods<SVDummyAllele> rll){
        final GenotypeLikelihoodCalculator genotypeLikelihoodCalculator = new GenotypeLikelihoodCalculators().getInstance(2, 2); // TODO: diploid, biallelic assumption
        final GenotypeLikelihoods result = genotypeLikelihoodCalculator.genotypeLikelihoods(  rll.sampleMatrix(0)  );
        logger.debug(".........GLC DONE.........");
        return result;
    }

    /**
     * TODO: this is possible only because of single sample; full-blown cohort genotyping for SV is to be done (SNP and indel models should be learned)
     * Simply infer genotype of the single diploid sample from the computed PL.
     * @return a new genotype based on the PL of input genotype
     */
    @VisibleForTesting
    static Genotype inferGenotypeFromPL(final Genotype gtWithPL, final List<Allele> alleles){

        final GenotypeBuilder builder = new GenotypeBuilder(gtWithPL);
        if(!GATKVariantContextUtils.isInformative(gtWithPL.getLikelihoods().getAsVector())){
            builder.noPL();
        }

        final double[] ll = MathUtils.normalizeFromLog10(gtWithPL.getLikelihoods().getAsVector(), false, true);

        GATKVariantContextUtils.makeGenotypeCall(ploidy, builder, GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN, ll, alleles);

        return builder.make();
    }

    // -----------------------------------------------------------------------------------------------
    // Utilities: debugging
    // -----------------------------------------------------------------------------------------------

    private void writeoutDebugReadLikelihoods(final JavaRDD<SVJunction> readLikelihoodsJavaPairRDD,
                                              final String debugStringSavingPath){
        readLikelihoodsJavaPairRDD.map(j -> j.debugString).saveAsTextFile(debugStringSavingPath);
    }

}
