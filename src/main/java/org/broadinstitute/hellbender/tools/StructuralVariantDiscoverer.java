package org.broadinstitute.hellbender.tools;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.BasicReference;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigAlignmentsArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AnnotatedVariantProducer;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SvType;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignedContig;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AssemblyContigWithFineTunedAlignments;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AssemblyContigWithFineTunedAlignments.AlignmentSignatureBasicType;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.ContigChimericAlignmentIterativeInterpreter;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.CpxVariantCanonicalRepresentation;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.CpxVariantInducingAssemblyContig;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.CpxVariantInterpreter;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.NovelAdjacencyAndAltHaplotype;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.SegmentedCpxVariantSimpleVariantExtractor;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.SegmentedCpxVariantSimpleVariantExtractor.MultiSegmentsCpxVariantExtractor;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.SegmentedCpxVariantSimpleVariantExtractor.RelevantAttributes;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.SegmentedCpxVariantSimpleVariantExtractor.ZeroAndOneSegmentCpxVariantExtractor;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.SimpleChimera;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.SimpleNovelAdjacencyAndChimericAlignmentEvidence;
import org.broadinstitute.hellbender.tools.spark.sv.utils.CNVInputReader;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.utils.SVIntervalTree;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVVCFWriter;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Tuple2;

import java.util.*;


@BetaFeature
@CommandLineProgramProperties(
        oneLineSummary = "(Internal) Examines aligned contigs from local assemblies and calls structural variants or their breakpoints",
        summary =
            "This tool takes a file containing the alignments of assembled contigs and searches it for contigs with" +
            " split alignments or large gaps indicating the presence of structural variation breakpoints." +
            " Variations' types are determined by analyzing the signatures of the split alignments," +
            " and are written to a VCF file.",
        programGroup = StructuralVariantDiscoveryProgramGroup.class)
public class StructuralVariantDiscoverer extends ReadWalker {
    @Argument(doc = "Name of output VCF.", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = "outputVCFName")
    private String outputVCFName;

    @Argument(doc = "file containing non-canonical chromosome names (e.g chrUn_KI270588v1) in the reference, " +
            "human reference (hg19 or hg38) assumed when omitted", shortName = "alt-tigs",
            fullName = "non-canonical-contig-names-file", optional = true)
    private String nonCanonicalChromosomeNamesFile;

    @ArgumentCollection
    private final DiscoverVariantsFromContigAlignmentsArgumentCollection discoverStageArgs =
            new DiscoverVariantsFromContigAlignmentsArgumentCollection();

    private static final double SCORE_DIFF_TOLERANCE = 0.;

    private String sampleId;
    private SAMSequenceDictionary refDict;
    private BasicReference reference;
    private SVIntervalTree<VariantContext> cnvCalls;
    private Set<String> canonicalChromosomes;

    private String currentContigName = null;
    private final List<GATKRead> readsForCurrentContig = new ArrayList<>();
    private final Map<NovelAdjacencyAndAltHaplotype, SimpleNovelAdjacencyAndChimericAlignmentEvidence>
            simpleEvidenceForNovelAdjacencyMap = new HashMap<>(10000);
    private final Map<CpxVariantCanonicalRepresentation, List<CpxVariantInducingAssemblyContig>>
            complexEventContigsMap =
            new HashMap<>(1000);
    private final List<AlignedContig> complexContigs = new ArrayList<>(1000);

    @Override public boolean requiresReads() { return true; }
    @Override public boolean requiresReference() { return true; }

    @Override public List<ReadFilter> getDefaultReadFilters() {
        return Arrays.asList(ReadFilterLibrary.MAPPED, ReadFilterLibrary.NOT_SECONDARY_ALIGNMENT);
    }

    @Override public void onTraversalStart() {
        final SAMFileHeader header = getHeaderForReads();
        if ( header.getSortOrder() != SAMFileHeader.SortOrder.queryname ) {
            throw new UserException("This tool requires a queryname-sorted source of reads.");
        }
        sampleId = SVUtils.getSampleId(header);
        refDict = header.getSequenceDictionary();
        cnvCalls = discoverStageArgs.cnvCallsFile == null ? null :
                CNVInputReader.loadCNVCalls(discoverStageArgs.cnvCallsFile, header);
        canonicalChromosomes = SVUtils.getCanonicalChromosomes(nonCanonicalChromosomeNamesFile, refDict);
    }

    @Override public void apply( final GATKRead read,
                                 final ReferenceContext referenceContext,
                                 final FeatureContext featureContext ) {
        reference = referenceContext;
        final String contigName = read.getName();
        if ( !contigName.equals(currentContigName) ) {
            if ( !readsForCurrentContig.isEmpty() ) {
                processContigAlignments(readsForCurrentContig);
                readsForCurrentContig.clear();
            }
            currentContigName = contigName;
        }
        readsForCurrentContig.add(read);
    }

    @Override public Object onTraversalSuccess() {
        final Object result = super.onTraversalSuccess();

        if ( !readsForCurrentContig.isEmpty() ) {
            processContigAlignments(readsForCurrentContig);
        }

        final List<VariantContext> variants =
                new ArrayList<>(2 * simpleEvidenceForNovelAdjacencyMap.size());
        final Collection<SimpleNovelAdjacencyAndChimericAlignmentEvidence> simpleEvidence =
                simpleEvidenceForNovelAdjacencyMap.values();
        for ( final SimpleNovelAdjacencyAndChimericAlignmentEvidence novelAdjacencyAndEvidence : simpleEvidence ) {
            final List<SvType> svTypes =
                    novelAdjacencyAndEvidence.getNovelAdjacencyReferenceLocations().toSimpleOrBNDTypes(reference);
            variants.addAll(novelAdjacencyAndEvidence.toVariantContexts(svTypes, sampleId, refDict, cnvCalls));
        }
        final ZeroAndOneSegmentCpxVariantExtractor zeroAndOneSegmentCpxVariantExtractor =
                new ZeroAndOneSegmentCpxVariantExtractor();
        final List<VariantContext> multiSegmentVariants = new ArrayList<>(complexEventContigsMap.size());
        for ( final Map.Entry<CpxVariantCanonicalRepresentation, List<CpxVariantInducingAssemblyContig>> entry :
                complexEventContigsMap.entrySet() ) {
            final VariantContext variantContext =
                    CpxVariantInterpreter.toVariantContext(entry.getKey(), entry.getValue(), reference);
            final int refSegs =
                    SVUtils.getAttributeAsStringList(variantContext, GATKSVVCFConstants.CPX_SV_REF_SEGMENTS).size();
            if ( refSegs < 2 ) {
                variants.addAll(zeroAndOneSegmentCpxVariantExtractor.extract(variantContext, reference));
            } else {
                multiSegmentVariants.add(variantContext);
            }
        }
        final Map<String, RelevantAttributes> contigNameToCpxVariantAttributes =
                new HashMap<>(2 * multiSegmentVariants.size());
        for ( final VariantContext variantContext : multiSegmentVariants ) {
            final RelevantAttributes relevantAttributes = new RelevantAttributes(variantContext);
            for ( final String contigName :
                    SVUtils.getAttributeAsStringList(variantContext, GATKSVVCFConstants.CONTIG_NAMES) ) {
                contigNameToCpxVariantAttributes.put(contigName, relevantAttributes);
            }
        }
        final Map<NovelAdjacencyAndAltHaplotype, List<SimpleChimera>> novelAdjacenciesToReinterpret = new HashMap<>();
        for ( final AlignedContig alignedContig : complexContigs ) {
            if ( contigNameToCpxVariantAttributes.containsKey(alignedContig.getContigName()) &&
                    alignedContig.getAlignments().size() > 1 ) {
                final List<SimpleChimera> chimeras =
                        ContigChimericAlignmentIterativeInterpreter.parseOneContig(
                                alignedContig,
                                refDict,
                                true,
                                DiscoverVariantsFromContigAlignmentsArgumentCollection.DEFAULT_MIN_ALIGNMENT_LENGTH,
                                DiscoverVariantsFromContigAlignmentsArgumentCollection.CHIMERIC_ALIGNMENTS_HIGHMQ_THRESHOLD,
                                true);
                for ( final SimpleChimera simpleChimera : chimeras ) {
                    final NovelAdjacencyAndAltHaplotype novelAdjacency =
                            new NovelAdjacencyAndAltHaplotype(simpleChimera, alignedContig.getContigSequence(), refDict);
                    final List<SimpleChimera> chimerasForNovelAdjacency =
                            novelAdjacenciesToReinterpret.get(novelAdjacency);
                    if ( chimerasForNovelAdjacency != null ) {
                        chimerasForNovelAdjacency.add(simpleChimera);
                    } else {
                        final List<SimpleChimera> newList = new ArrayList<>(2);
                        newList.add(simpleChimera);
                        novelAdjacenciesToReinterpret.put(novelAdjacency, newList);
                    }
                }
            }
        }
        final List<VariantContext> reinterpretedVariants = new ArrayList<>(novelAdjacenciesToReinterpret.size());
        for ( final Map.Entry<NovelAdjacencyAndAltHaplotype, List<SimpleChimera>> entry :
                novelAdjacenciesToReinterpret.entrySet() ) {
            reinterpretedVariants.add(
                new SimpleNovelAdjacencyAndChimericAlignmentEvidence(entry.getKey(), entry.getValue())
                        .produceAnnotatedVcFromAssemblyEvidence(
                                ContigChimericAlignmentIterativeInterpreter
                                        .inferSimpleTypeFromNovelAdjacency(entry.getKey(), reference),
                                refDict, cnvCalls, sampleId).make());
        }
        final List<VariantContext> extractedMultiSegmentVariants = new ArrayList<>(multiSegmentVariants.size());
        final MultiSegmentsCpxVariantExtractor multiSegmentsCpxVariantExtractor =
                new MultiSegmentsCpxVariantExtractor();
        for ( final VariantContext variantContext : multiSegmentVariants ) {
            extractedMultiSegmentVariants.addAll(multiSegmentsCpxVariantExtractor.extract(variantContext, reference));
        }
        final List<VariantContext> consistentVariants =
                SegmentedCpxVariantSimpleVariantExtractor
                        .filterForConsistency(reinterpretedVariants, contigNameToCpxVariantAttributes, reference);
        variants.addAll(SegmentedCpxVariantSimpleVariantExtractor.removeDuplicates(extractedMultiSegmentVariants, consistentVariants));

        final List<VariantContext> filteredMergedVariantList =
                AnnotatedVariantProducer.filterMergedVariantList(variants, discoverStageArgs);
        SVVCFWriter.writeVCF(filteredMergedVariantList, outputVCFName, refDict, getDefaultToolVCFHeaderLines(), logger);
        return result;
    }

    private void processContigAlignments( final List<GATKRead> contigAlignments ) {
        final List<AlignmentInterval> alignmentIntervals = new ArrayList<>(contigAlignments.size());
        String contigName = null;
        byte[] contigSequence = null;
        for ( final GATKRead read : contigAlignments ) {
            contigName = read.getName();
            if ( !read.isSupplementaryAlignment() ) {
                contigSequence = read.getBasesNoCopy();
                if ( read.isReverseStrand() ) {
                    contigSequence = BaseUtils.simpleReverseComplement(contigSequence);
                }
            }
            alignmentIntervals.add(new AlignmentInterval(read));
        }
        if ( contigSequence == null ) {
            throw new UserException("No primary line for " + contigName);
        }
        final AlignedContig alignedContig = new AlignedContig(contigName, contigSequence, alignmentIntervals);
        if ( !alignedContig.hasGoodMQ() ) return;

        final List<AssemblyContigWithFineTunedAlignments> fineTunedAlignmentsList =
            alignedContig.reconstructContigFromBestConfiguration(canonicalChromosomes, SCORE_DIFF_TOLERANCE);
        for ( final AssemblyContigWithFineTunedAlignments fineTunedAlignment : fineTunedAlignmentsList ) {
            if ( fineTunedAlignment.getAlignmentSignatureBasicType() == AlignmentSignatureBasicType.SIMPLE_CHIMERA ) {
                if ( SimpleChimera.splitPairStrongEnoughEvidenceForCA(fineTunedAlignment.getHeadAlignment(),
                                                                      fineTunedAlignment.getTailAlignment()) ) {
                    final SimpleChimera simpleChimera = fineTunedAlignment.extractSimpleChimera(refDict);
                    final NovelAdjacencyAndAltHaplotype novelAdjacency =
                            new NovelAdjacencyAndAltHaplotype(simpleChimera,
                                                                fineTunedAlignment.getContigSequence(),
                                                                refDict);
                    final SimpleNovelAdjacencyAndChimericAlignmentEvidence simpleEvidence =
                            simpleEvidenceForNovelAdjacencyMap.get(novelAdjacency);
                    if ( simpleEvidence != null ) {
                        simpleEvidence.getAlignmentEvidence().add(simpleChimera);
                    } else {
                        final SimpleNovelAdjacencyAndChimericAlignmentEvidence novelAdjacencyEvidence =
                                new SimpleNovelAdjacencyAndChimericAlignmentEvidence(novelAdjacency,
                                                                Collections.singletonList(simpleChimera));
                        simpleEvidenceForNovelAdjacencyMap.put(novelAdjacency, novelAdjacencyEvidence);
                    }
                }
            } else if ( fineTunedAlignment.getAlignmentSignatureBasicType() == AlignmentSignatureBasicType.COMPLEX ) {
                complexContigs.add(fineTunedAlignment.getSourceContig());
                final Tuple2<CpxVariantCanonicalRepresentation, CpxVariantInducingAssemblyContig> entry =
                                CpxVariantInterpreter.getOneVariantFromOneContig(fineTunedAlignment, refDict);
                final List<CpxVariantInducingAssemblyContig> contigsList = complexEventContigsMap.get(entry._1);
                if ( contigsList != null ) {
                    contigsList.add(entry._2);
                } else {
                    final List<CpxVariantInducingAssemblyContig> newList = new ArrayList<>(2);
                    newList.add(entry._2);
                    complexEventContigsMap.put(entry._1, newList);
                }
            }
        }
    }
}
