package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SVTestUtils;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVType;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.*;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import scala.Tuple2;

import java.io.IOException;
import java.nio.file.Files;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants.*;

/**
 * For providing test data and utils specifically for complex sv.
 * Any attempts to reuse these resources elsewhere,
 * even for simple SV's
 * (see {@link org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVDiscoveryTestDataProvider} for that purpose)
 * should be take extreme care.
 */
public final class CpxSVInferenceTestUtils extends GATKBaseTest {

    // get their de-overlapped alignments, basic info, jumps, segments, and description, and finally, VCF records
    static final class PreprocessedAndAnalysisReadyContigWithExpectedResults {

        final CpxVariantInducingAssemblyContig expectedCpxVariantInducingAssemblyContig;

        final CpxVariantCanonicalRepresentation expectedCpxVariantCanonicalRepresentation;

        final VariantContext expectedVariantContext;
        final byte[] assumedReferenceSequence; // reference sequence used in constructed the expected results that can be used in real tests

        final Set<SimpleInterval> expectedTwoBaseBoundaries;

        PreprocessedAndAnalysisReadyContigWithExpectedResults(final CpxVariantInducingAssemblyContig expectedCpxVariantInducingAssemblyContig,
                                                              final CpxVariantCanonicalRepresentation expectedCpxVariantCanonicalRepresentation,
                                                              final VariantContext expectedVariantContext,
                                                              final byte[] assumedReferenceSequence,
                                                              final Set<SimpleInterval> expectedTwoBaseBoundaries) {
            this.expectedCpxVariantInducingAssemblyContig = expectedCpxVariantInducingAssemblyContig;
            this.expectedCpxVariantCanonicalRepresentation = expectedCpxVariantCanonicalRepresentation;
            this.expectedVariantContext = expectedVariantContext;
            this.assumedReferenceSequence = assumedReferenceSequence;
            this.expectedTwoBaseBoundaries = expectedTwoBaseBoundaries;
        }
    }

    static final List<PreprocessedAndAnalysisReadyContigWithExpectedResults> PREPROCESSED_AND_ANALYSIS_READY_CONTIGS_AND_EXPECTED_RESULTS;


    public final static Set<String> hg38CanonicalChromosomes;
    /**
     * We are having this because it is SV, especially complex ones, are rare and events on chr20 and 21 are not enough.
     */
    final static SAMSequenceDictionary bareBoneHg38SAMSeqDict;
    static {
        final List<SAMSequenceRecord> hg38Chromosomes = new ArrayList<>();
        final String hg38ChrBareBoneListFile =  GATKBaseTest.toolsTestDir + "/spark/sv/utils/hg38ChrBareBone.txt";
        try (final Stream<String> records = Files.lines(IOUtils.getPath(( hg38ChrBareBoneListFile )))) {
            records.forEach(line -> {
                final String[] fields = line.split("\t", 2);
                hg38Chromosomes.add(new SAMSequenceRecord(fields[0], Integer.valueOf(fields[1])));
            });
            bareBoneHg38SAMSeqDict = new SAMSequenceDictionary(hg38Chromosomes);
        } catch ( final IOException ioe ) {
            throw new UserException("Can't read nonCanonicalContigNamesFile file " + hg38ChrBareBoneListFile, ioe);
        }

        hg38CanonicalChromosomes =
                hg38Chromosomes.subList(0, 24).stream().map(SAMSequenceRecord::getSequenceName).collect(Collectors.toCollection(LinkedHashSet::new));


        PREPROCESSED_AND_ANALYSIS_READY_CONTIGS_AND_EXPECTED_RESULTS = new ArrayList<>();
        PREPROCESSED_AND_ANALYSIS_READY_CONTIGS_AND_EXPECTED_RESULTS.add(buildForContig_2428_0(hg38CanonicalChromosomes, bareBoneHg38SAMSeqDict));
        PREPROCESSED_AND_ANALYSIS_READY_CONTIGS_AND_EXPECTED_RESULTS.add(buildForContig_2548_1(hg38CanonicalChromosomes, bareBoneHg38SAMSeqDict));
        PREPROCESSED_AND_ANALYSIS_READY_CONTIGS_AND_EXPECTED_RESULTS.add(buildForContig_30077_1(hg38CanonicalChromosomes, bareBoneHg38SAMSeqDict));
        PREPROCESSED_AND_ANALYSIS_READY_CONTIGS_AND_EXPECTED_RESULTS.add(buildForContig_22672_4(hg38CanonicalChromosomes, bareBoneHg38SAMSeqDict));
    }

    // =================================================================================================================
    // values below are hand-picked and results manually curated.

    private static PreprocessedAndAnalysisReadyContigWithExpectedResults buildForContig_2428_0(final Set<String> canonicalChromosomes,
                                                                                               final SAMSequenceDictionary refSeqDict) {
        final AssemblyContigWithFineTunedAlignments analysisReadyContig =
                makeContigAnalysisReady("asm002428:tig00000\t16\tchr2\t4450935\t60\t1216M3I256M1030S\t*\t0\t0\tTTCCTCTGTGGAGGACACAGCATTCAAGGCGCCATCTTGGAATCAGAATCACTAACCACATCAAAGCCTAATACTGGACTTCCCAACCTCCTAAACTGTAAGCAAATAATTTTGTTCATAAGTTACCCAGCCTGTGCTATTCTTCTATACCAGCACATATGGACTAAGACACTACGTATTCTAATAAGTGTGTGAAAACCCAAGATAAAATAGGCGTTATTTCCCATTTTACAGATGAGAAAACTGAAGCTTGAGGAATTTGTGAATTTGTTTAATCCAAATCTGTTAGAAAATAGTAGACCTCACCTCTTCAATAGTTGAGCATGGATAATTATGATACTAATTTTAGACTCATATTTCTTCATTTGTACATGATTTAAATCTCAGAAAAAATTTAGAATTTCTCTTTACTGTCTGTTGGTATTATGAAATGTGATAGCAGTATGTCTAGGCACAGATTTTTTTTTTTTTTTTCCTTAGCAGGGGTAGTTATTTTTTTCTTTCTTAGCAAGGGTGGCTTTTTCAGTCTGAAAGCTGAGACATACACTTTTCTCTGAGAAGAAATATTTTCTATTATTTGAATGTTTCACCTCTCCTTTTTCCTATTCTCTCTTTCTGAGACATACATTGATTCTCCAAGTTATCATCTTTCTCTCTTATCTCTATGTTTGTCTTTGGTTCTCTTTTTTCTCTGTATTCTGGGAAATTTTCATAGCTTTATTCAAGTTCATTTTTCAAGCTCTAAGAACTCTCTTTTTGTCTGATTTTACCATTTTTATTATACTAGTTGTTACAGTTTTGATAGATAAACTATCTTCCAAAATTTCTAAACGTATGGGTTTGGGAAATTGCTTTTAAAATCCTTAGGGATATAGCCAAATTACAGCATGCCAAGAAATGACAGAAATGTATTTTTTAATTAAAAAAAGGAGTATTAGTAAAGCTTCAACTTTAGAACAGGGTCTATGTCAGAGAAACACCACAGAGCAATGAAAACTGTGTTCTCCAAATTATTTAGAACGTAGGAGGATTTCCAAAATATTTATTTATTTTGAAAGTAGAATGATTAACATCTGCTGTCTATAAGCAACACTTTTGAAAAAAGGAATTTACTCCCTCCAATTTGCTCCATAGCCTACAGTTTCTTCTTTTATCTACCTCTTCCTCCTCCTCCTCCACCACCTTCTTTTCTTGTTGTTGTTGTTGTTCTTCTTCTTCCTCCTCCTCCTCTTCCTCTTCTTCTTCTTCTCCTTCTTCTTCCCCTCCCTTCCTTCCTTCCTTCTTTCCTTCCTTCCTTCCTCTTGTTCTTCTTCTTCCTCTTATTCTCCTTCTCCTTCTTCTTCTTCTCCTGCTTCTTCCTCTTCTCCTTCTTCTTCTTCCTCCTCCTCCTTCTTCTCTTCTTCTGCTTTTCTTCTCTTCTTCTTCTCTTTTCCTTCCTTCCTTCTTTCCTTCCTCTTCTTTCTCTTCTTCTTCCTCTTATTCTCCTTCTTCCTCCTCTTCTCCTTCTTCCTTCTCCTCCTCTTCCTCCTCGTTGTTCTTGTCATTGTTGTTTTTGTTGTTCTTCCTCCTCTTCTCTCTCCTCCTTCTCCTCCTCCTCCTCCTTCTTCTTCTCCTTGTCCTCCTCCTCTGTCTCTGTCTTCTTCTTTTTCTTCTTGTTCCTCTTCTTCTTCCTCTTATTCTCCTTCTTCTTCCTCCCCTTCTCCCTCTTCCGTCTTCTCCTCCTCCTCCCTCTTCCTTCTTCTCCTCCTCCTCCTTCTTTTTCTTTCTTCTTCTTTCTTCTTGTTCTTGTTCTTCTTCTTCTTCCTCCTCCTCATCCTTCTTCTCTTCTTCTGCTTTTCTTCTGTTCTTCTTCTTCTTTCCCTCCCTCCCTCCCTCCCTCCCTTCCTTCCTTCCTTCCTTCCTCTTCTTTCTCCTCTTCTTCTTCTTGTTCTTGTTCTTCTTCTTCCTCTTATTCTCCTTCTTCCTCCTCTTCTCCTTCTTCTTCCTTTTGCTCCTCTTCCTCCTTCTTGTTGTTGTTGTTCTTCTTCTTCTTCCTCTCATTCTTCTTCTTCCTCCTCTTCTCCCTCCTCCTTCTCCTCCTCCTCCTTCTTCTTCTCCTTCTCCTCCTCCTTTTTCTTCTCCTTCTCCTCCTCTTCCATCTCTGTCTTTGCCTTCTTCTTCTTGTTCCTCATCTTCTTCTTCCTCTTATTCTCCTTCTTCTTCCTTCCCTTCTCCCTCTTCCTTTTTCTCCTCCTCCTCCTTCTTTCTTCTTTCTTCCTTCTTCTTCTTTCTCTTCTTCTTCTTCTCCTTCTTCTTCTTCCTCCTCCTCCTCCCCTTCCCCCTTCCCCTCCTCCTCCTCCTTCTCCTCCTTCTCCTCCTCCTTCTTCTTCTTCCTCTTCTTCTTCTTCTTCTCCTTCTTCTTCTTCTCCTTCTTCTTCTTCTTCCTCCTCCTCCTCCCCTTCTCCCTTCCCCTCCTCCTCCTCCTCCTTCTCCTCCTCCTTCTTCTTCTTCTTCTTCTTCTTCTTC\t*\tSA:Z:chr2,4452399,-,2250H75M3D52M15I113M,60,25,157;chr9,34840411,-,2065H67M373H,33,2,57;chr10,75714830,-,1851H51M603H,0,0,51;chr2,4452298,-,1793H60M652H,60,3,45;chr6,15853372,-,1907H43M555H,60,0,43;\tMD:Z:259C382G365G125G337\tRG:Z:GATKSVContigAlignments\tNM:i:7\tAS:i:1433\tXS:i:0", canonicalChromosomes, refSeqDict);

        final CpxVariantInducingAssemblyContig.BasicInfo manuallyCalculatedBasicInfo =
                new CpxVariantInducingAssemblyContig.BasicInfo("chr2", false,
                        new SimpleInterval("chr2", 4450935, 4450935), new SimpleInterval("chr2", 4452641, 4452641));

        //////////
        final List<CpxVariantInducingAssemblyContig.Jump> manuallyCalculatedJumps =
                Arrays.asList(
                        new CpxVariantInducingAssemblyContig.Jump(new SimpleInterval("chr2", 4452399, 4452399), new SimpleInterval("chr9", 34840477, 34840477), StrandSwitch.NO_SWITCH, 118),
                        new CpxVariantInducingAssemblyContig.Jump(new SimpleInterval("chr9", 34840411, 34840411), new SimpleInterval("chr6", 15853414, 15853414), StrandSwitch.NO_SWITCH, 115),
                        new CpxVariantInducingAssemblyContig.Jump(new SimpleInterval("chr6", 15853372, 15853372), new SimpleInterval("chr2", 4452357, 4452357), StrandSwitch.NO_SWITCH, 54),
                        new CpxVariantInducingAssemblyContig.Jump(new SimpleInterval("chr2", 4452298, 4452298), new SimpleInterval("chr2", 4452406, 4452406), StrandSwitch.NO_SWITCH, 318)
                );
        final List<SimpleInterval> manuallyCalculatedEventPrimaryChromosomeSegmentingLocations =
                Arrays.asList(
                        new SimpleInterval("chr2", 4452298, 4452298),
                        new SimpleInterval("chr2", 4452357, 4452357),
                        new SimpleInterval("chr2", 4452399, 4452399),
                        new SimpleInterval("chr2", 4452406, 4452406)
                );
        final CpxVariantInducingAssemblyContig manuallyCalculatedCpxVariantInducingAssemblyContig =
                new CpxVariantInducingAssemblyContig(analysisReadyContig,
                        manuallyCalculatedBasicInfo,
                        manuallyCalculatedJumps,
                        manuallyCalculatedEventPrimaryChromosomeSegmentingLocations);

        //////////
        final SimpleInterval manuallyCalculatedAffectedRefRegion =
                new SimpleInterval("chr2", 4452298, 4452406);
        final List<SimpleInterval> manuallyCalculatedSegments =
                Arrays.asList(
                        new SimpleInterval("chr2", 4452298, 4452357),
                        new SimpleInterval("chr2", 4452357, 4452399),
                        new SimpleInterval("chr2", 4452399, 4452406));
        final List<String> manuallyCalculatedAltArrangementDescription =
                Arrays.asList("1", "2", "3", "UINS-318", "1", "UINS-54", "chr6:15853372-15853414", "UINS-115", "chr9:34840411-34840477", "UINS-118", "3");
        final byte[] manuallyCalculatedAltSeq = Arrays.copyOfRange(analysisReadyContig.getContigSequence(), 247, 1139);
        SequenceUtil.reverseComplement(manuallyCalculatedAltSeq);
        final CpxVariantCanonicalRepresentation manuallyCalculatedCpxVariantCanonicalRepresentation =
                new CpxVariantCanonicalRepresentation(
                        manuallyCalculatedAffectedRefRegion,
                        manuallyCalculatedSegments,
                        manuallyCalculatedAltArrangementDescription,
                        manuallyCalculatedAltSeq);

        //////////
        final byte[] dummyRefSequence = "TCGA".getBytes();
        final VariantContextBuilder variantContextBuilder = new VariantContextBuilder()
                .chr("chr2").start(4452298).stop(4452406)
                .alleles(Arrays.asList(Allele.create(dummyRefSequence, true), Allele.create(SimpleSVType.createBracketedSymbAlleleString(CPX_SV_SYB_ALT_ALLELE_STR), false)))
                .id("CPX_chr2:4452298-4452406")
                .attribute(VCFConstants.END_KEY, 4452406)
                .attribute(SVTYPE, GATKSVVCFConstants.CPX_SV_SYB_ALT_ALLELE_STR)
                .attribute(SVLEN, 783)
                .attribute(SEQ_ALT_HAPLOTYPE, new String(manuallyCalculatedAltSeq));
        variantContextBuilder.attribute(CPX_EVENT_ALT_ARRANGEMENTS,
                String.join(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR, manuallyCalculatedAltArrangementDescription));
        variantContextBuilder.attribute(CPX_SV_REF_SEGMENTS,
                String.join(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR, manuallyCalculatedSegments.stream().map(SimpleInterval::toString).collect(Collectors.toList())));
        final VariantContext manuallyCalculatedVariantContext = variantContextBuilder.make();

        //////////
        final Set<SimpleInterval> manuallyCalculatedTwoBaseBoundaries = new HashSet<>();
        manuallyCalculatedTwoBaseBoundaries.add(new SimpleInterval("chr2", 4452399, 4452400));
        manuallyCalculatedTwoBaseBoundaries.add(new SimpleInterval("chr2", 4452640, 4452641));
        manuallyCalculatedTwoBaseBoundaries.add(new SimpleInterval("chr2", 4452298, 4452299));
        manuallyCalculatedTwoBaseBoundaries.add(new SimpleInterval("chr2", 4452356, 4452357));
        manuallyCalculatedTwoBaseBoundaries.add(new SimpleInterval("chr2", 4450935, 4450936));
        manuallyCalculatedTwoBaseBoundaries.add(new SimpleInterval("chr2", 4452405, 4452406));

        //////////
        return new PreprocessedAndAnalysisReadyContigWithExpectedResults(
                manuallyCalculatedCpxVariantInducingAssemblyContig,
                manuallyCalculatedCpxVariantCanonicalRepresentation,
                manuallyCalculatedVariantContext,
                dummyRefSequence,
                manuallyCalculatedTwoBaseBoundaries);
    }


    private static PreprocessedAndAnalysisReadyContigWithExpectedResults buildForContig_2548_1(final Set<String> canonicalChromosomes,
                                                                                               final SAMSequenceDictionary refSeqDict) {
        final AssemblyContigWithFineTunedAlignments
                analysisReadyContig = makeContigAnalysisReady("asm002548:tig00001\t16\tchr2\t16224630\t60\t513M634S\t*\t0\t0\tTCAGCACAATGTCTGCTATCCAGTCACTGCTCTGTATTAGCCTCGCCTTAAATTGCTCGAATAGTAAGCAAGGGTTTATAAAGCATCTGTGTGCTAATCCCATCCTAAAGTGACAGGGATCAATTTAACCTCCTCATTTAAAAGATGAGGAGACAGGCTCAAAGGGCAAAGCCGACTCACCCAGGCTCACACAACCAGAGTGGCAGAGCCTCAGGGCTTGCCCTGCACAGCGTAGCTTCTGTTCTGAGAGGAATGCACCTGGGTTTTGTGTCTTTGGTGCTGGGGTGGTGAATGGGAGGCTCTGACCGATGGCAGAGTTCCCAGGGGATCCTTCAGGCAGCCGGGGGGCCTTTCTTATGGCAGGGCTCTGAAGCAAGTCCTGGCTGACCTTTTCTTGAAGCTGAGCTCTTCTGCCAAGACCCTTGGCCCTGCCCTTGGCCCCTGTAGAAATGCACACCTGTAAGCCTGTAATACCTTCTCCAGATCACAGCAACTAGGGAACTGCTGGAAGCCACGGAGATCAGCAGGGACCGCACCATCCCCAGAAGTGTGGATCCACCCGAGGGTCACGTGGATCCTCAAGGAGTGGAGAGTAGACCCTGAGCCACGTGGGCTTCTAGCAGTTCCCTGGGCCTGCTGGCACCAAGATGGGCTCCCTGTCACAGGCCTGACGGGACCTCTGTGAGGGTCCCAGGAGATGAAGCACATGAAGAGGGTTGGTGCACCAAGCAGGCTCACCAAATATTCACTCCCTGACCCCTTTCCCCACGCAAGCTGTGGTGGGGAAGGAGTGTTTCATGGTGGGCAGCCCAGGTGCCAGTCCCCCCCTTCCCCAGTGATGGGGTACCATTACCCCCACGACAGCCCAGCTTGGTGCCAGCAGGCCCAGGGAACTGCTGGAAGCCCACGTGGCTCAGGTCTACTCCTCACTCCTTGAGGATCCACGTGGCCCTTGGGTGGATCCACACTTCTGGGGATGGTGCGGTCCCAGCTGATCTCCCTGGGGGTCCAGTTAGCCTCATCCTCCTTCACCCAGGGTGGCTGCATGTCAGCACCTGCTCTACCTTCTACCTTAATGGCGTCCCCTAGGCTATCTGTCAGGGGAGTCTATGCAAAGGGGCCCCTTTGAACTTGCCTGTGCCACTCC\t*\tSA:Z:chr2,16225125,-,887H260M,60,0,260;chr2,16226497,-,657H170M1D53M267H,60,14,141;chr2,16225125,+,518H29M1I84M515H,60,7,66;\tMD:Z:513\tRG:Z:GATKSVContigAlignments\tNM:i:0\tAS:i:513\tXS:i:0",
                canonicalChromosomes, refSeqDict);

        final CpxVariantInducingAssemblyContig.BasicInfo manuallyCalculatedBasicInfo =
                new CpxVariantInducingAssemblyContig.BasicInfo("chr2", false,
                        new SimpleInterval("chr2", 16224630, 16224630), new SimpleInterval("chr2", 16225384, 16225384));

        //////////
        final List<CpxVariantInducingAssemblyContig.Jump> manuallyCalculatedJumps =
                Arrays.asList(
                        new CpxVariantInducingAssemblyContig.Jump(new SimpleInterval("chr2", 16225125, 16225125), new SimpleInterval("chr2", 16226720, 16226720), StrandSwitch.NO_SWITCH, 7),
                        new CpxVariantInducingAssemblyContig.Jump(new SimpleInterval("chr2", 16226497, 16226497), new SimpleInterval("chr2", 16225125, 16225125), StrandSwitch.REVERSE_TO_FORWARD, 28),
                        new CpxVariantInducingAssemblyContig.Jump(new SimpleInterval("chr2", 16225237, 16225237), new SimpleInterval("chr2", 16225142, 16225142), StrandSwitch.FORWARD_TO_REVERSE, 2)
                );
        final List<SimpleInterval> manuallyCalculatedEventPrimaryChromosomeSegmentingLocations =
                Arrays.asList(
                        new SimpleInterval("chr2", 16225125, 16225125),
                        new SimpleInterval("chr2", 16225142, 16225142),
                        new SimpleInterval("chr2", 16225237, 16225237)
                );
        final CpxVariantInducingAssemblyContig manuallyCalculatedCpxVariantInducingAssemblyContig =
                new CpxVariantInducingAssemblyContig(analysisReadyContig,
                        manuallyCalculatedBasicInfo,
                        manuallyCalculatedJumps,
                        manuallyCalculatedEventPrimaryChromosomeSegmentingLocations);

        //////////
        final SimpleInterval manuallyCalculatedAffectedRefRegion =
                new SimpleInterval("chr2", 16225125, 16225237);
        final List<SimpleInterval> manuallyCalculatedSegments =
                Arrays.asList(
                        new SimpleInterval("chr2", 16225125, 16225142),
                        new SimpleInterval("chr2", 16225142, 16225237));
        final List<String> manuallyCalculatedAltArrangementDescription =
                Arrays.asList("1", "UINS-2", "-2", "-1", "UINS-28", "chr2:16226497-16226720", "UINS-7", "1", "2");
        final byte[] manuallyCalculatedAltSeq = Arrays.copyOfRange(analysisReadyContig.getContigSequence(), 147, 652);
        SequenceUtil.reverseComplement(manuallyCalculatedAltSeq);
        final CpxVariantCanonicalRepresentation manuallyCalculatedCpxVariantCanonicalRepresentation =
                new CpxVariantCanonicalRepresentation(
                        manuallyCalculatedAffectedRefRegion,
                        manuallyCalculatedSegments,
                        manuallyCalculatedAltArrangementDescription,
                        manuallyCalculatedAltSeq);

        //////////
        final byte[] dummyRefSequence = "TCGA".getBytes();
        final VariantContextBuilder variantContextBuilder = new VariantContextBuilder()
                .chr("chr2").start(16225125).stop(16225237)
                .alleles(Arrays.asList(Allele.create(dummyRefSequence, true), Allele.create(SimpleSVType.createBracketedSymbAlleleString(CPX_SV_SYB_ALT_ALLELE_STR), false)))
                .id("CPX_chr2:16225125-16225237")
                .attribute(VCFConstants.END_KEY, 16225237)
                .attribute(SVTYPE, GATKSVVCFConstants.CPX_SV_SYB_ALT_ALLELE_STR)
                .attribute(SVLEN, 392)
                .attribute(SEQ_ALT_HAPLOTYPE, new String(manuallyCalculatedAltSeq));
        variantContextBuilder.attribute(CPX_EVENT_ALT_ARRANGEMENTS,
                String.join(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR, manuallyCalculatedAltArrangementDescription));
        variantContextBuilder.attribute(CPX_SV_REF_SEGMENTS,
                String.join(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR, manuallyCalculatedSegments.stream().map(SimpleInterval::toString).collect(Collectors.toList())));
        final VariantContext manuallyCalculatedVariantContext = variantContextBuilder.make();

        //////////
        final Set<SimpleInterval> manuallyCalculatedTwoBaseBoundaries = new HashSet<>();
        manuallyCalculatedTwoBaseBoundaries.add(new SimpleInterval("chr2", 16225125, 16225126));
        manuallyCalculatedTwoBaseBoundaries.add(new SimpleInterval("chr2", 16225383, 16225384));
        manuallyCalculatedTwoBaseBoundaries.add(new SimpleInterval("chr2", 16225236, 16225237));
        manuallyCalculatedTwoBaseBoundaries.add(new SimpleInterval("chr2", 16224630, 16224631));
        manuallyCalculatedTwoBaseBoundaries.add(new SimpleInterval("chr2", 16225141, 16225142));

        //////////
        return new PreprocessedAndAnalysisReadyContigWithExpectedResults(
                manuallyCalculatedCpxVariantInducingAssemblyContig,
                manuallyCalculatedCpxVariantCanonicalRepresentation,
                manuallyCalculatedVariantContext,
                dummyRefSequence,
                manuallyCalculatedTwoBaseBoundaries);
    }

    private static PreprocessedAndAnalysisReadyContigWithExpectedResults buildForContig_30077_1(final Set<String> canonicalChromosomes,
                                                                                                final SAMSequenceDictionary refSeqDict) {
        final AssemblyContigWithFineTunedAlignments
                analysisReadyContig = makeContigAnalysisReady("asm030077:tig00001\t0\tchrX\t30792651\t60\t791M78I1312M96I610M51D129M1802S\t*\t0\t0\tCTTGCCGCTTGCCTTGCTTTGCCTTGCCTTGACTAGGGCCTTGCCTTGCCTATTGCCTTGACTTGCCTAGCTTTGTGCCTTGCCTTGCCTGGTGCCTTGACTTGCCTTGTGCCTTGCCTTCCCTTCGCTTGCCTTTTGCTTGCCTTGCCTTGCCTTGTGCCTTCCCTTGCCTTGCCTTGCCTTGTGCCTTGATTTGCCTTGCCTTGGCTTGTGCCTTGCCTTGCTTTGCCTTGTACTTGCCTTGCCTTGTGCCTTGCCTTGTCTTGCTTATGCCATGCCTTGCTTATGCCTTGCCTTGCCTTGTGCCTTGCCTTGTGCCCTCCCTTTCCTTGTGCCTTGCGTTTTTCTGTGCCTTACCTTGCCTTGCCTCCTGCCTTGCTTTGGCTTGCCATGACTAGGGCCTTGCCTTGTGCCTTGCCTTGCCTTGCGTCCTCTCTTGCCTTGTGCCTTGCATTGTTCTGTGCGTTACCTTGCCTTGCCTCCTGCCTTGCTTTGGCTTGCCATGACTAGGGCCTTGCCTTGTTGTGCCTTGCTTTCTGCCTTGCCTTGACTTGACTTGCCTTGCCTTGCTTTTTGCTTTGCCTTGCCATTCCTTTGTGCCTTGCCTTGACTTGCCTTATGCCTTGCCTTTTGCCTTCACTTGCCTTGCTTTGTGCCTCTCCTTGCCTGGTGCCTTGCCTTGCTCTGCCTTGCCTTGCTTTGCCTTGTGCCTTGCCTTGCATTGCTTTGCCTTGCCTTGTGCCTAGCATTGCCTTGCATTTTGCTTGCCTTTCCTTGACTTGTGTGTCAATTTGCCTTGCCTTGTTACTTGCCTTCCCTTGGCTTGTGCCTTGCCTTGCTTTTTGCCTTGCCTTGAATTGCGCATCGAGTTGCCTTGCCTTGTTACTTGTCTTCCCTTGGCTTGTGCCTTGCCTTGTGCCTTGCCTTACCTTGCCTTGTGTCTTGATTTGCCTTGCCTTGTGCCTTGCCTTGCCTTGCCTTGCCTTGCCTTGTGTTTTTTCTTGCCTTGTTCCATGCTTTACCTTGCCTTGCCTCACCTCGTGCCTTGCTTTGGCTTGCCTTGACTTGGACCTTGCCTTGCTGTGCCTTACTTTCTGCCTAGCCTAGCCTAGGCTTGCTTTTTGCTTTGCCTTGCCATTCCCTTGTGCTTTGCCTTGACTTGCCTTGCACCTTGCCTTGCCTTGCCTTGCCTTTTGCCTTCACTTGTTTTGCTTTGTGCATTGCATTGCCTTGCCTTGCCTGGTCTCTTGACTTGCCTTGTGCCTTGCCTTCCCTTTGCTTGCCTTTTGCTTGCCTTGCCTTGCCTTGTGCCTGGCCTTGGCCTTGCCTTGCTTGGCTTTGCCTTGCCTTGTGCCTTGCCTTTTGCCTTGCCTTCTGCCTTGTGCCTTGCACCTTGTCTTGCCTTGCTTGTGCCTTACCTTGCCTTGTGTTTTACCTTGCGTTGTCTTGCCTTGGCCTTGTCTTGCCTAGCTGTGCCTTGCCTTGTGCCTTGTTTTGCCTTGTGCCTTGCCTTGTGCCTTGTTTTGCCTTGTGCCTTGCCTTGTGCCTTGTTTTGCCTTGTGCCTTGCCTTGCTTTGGCCTTGCCTTGTCTGGCTTTGCCTTGATTTGTGCCTTGACTTGTCGTGCCTTGCCTTGTACCTTGCCTTGTGCTTTCCCTTGCTTTGCCTTGTGCCTTGATTTGCCTTGTGTTGTGCCTTGCCTTGCCATCCCTTGTGCCTTGCCTTGCTTTGCCTTGCCTTGTAACTTGCCTTGCCTTGCCTTATACCTTGCCTTACTTTGCCTTGCCATGCTTTCTGCCTTGCCTTGTCTTGTCTTGTGCCTTGCCTTGTATTGCTTATCCTTGGCTTATGCCTTGCATTGCCTTGCCTTGTATTGCTTATCCTTGGCTTATGCCTTGCATTGCCTTGCCTTGTGACTTACCTTGCCTTGCCCTGTGTGTTGATTTTTTTGCCTTGCCTTGTTACTTGCCTTCCCTTGGCTTGTGCCTTACCTTGCCTTGTTTGCGCCTTGCCTTGCCTGTGCCTTGCCTTGCCTTGCCTTGTGCCTTGCCTTGATTTGGCCTTGCTTTGCCTGGCTTTGCCTTGCCTTGTGCCTTGCCTTGCCTTGCCTTTGCTTTGCCTTGGCTTGCCTTGCCTTGTGCCTTGCATTGCCTTGCCTTACGCCTTGCATTGCCTTGCCTTTTGCCTTGCCTTGCCTTGTAACTTGCCTTGCCTTGCCATTCCTTGTCTTCCTGTGCCTTGATTTGCCTTACCTTGTGCCTTGTCTTGCCTTGGCCTGTGCCTTGCCTTGCCTTGTAACTTGCCTTGCCTTGCCATGCCTTGCCTTGCCATTCCTTGTCTTTCCTTGTGTTTCCTTGCCTTGTGCCTTGCCTTGTTTTGTGCTTTGCCTTGCCTTGCTTGTGCCTTGAATTGCCTTGTGCCTTGACTTGCCTTGCCTTGCCTTGTGCCTTGTCTTGACTTCTTCTGTGCCTTACCTTGCCTTACCTTGTGCCTCATTTGGCTTGCCTTGCCTTGGGCTCTGCCTTTCTGTGTCTTGCTCTTTGCTTTGTTTTGACGTGCCTTGCCTTGTGCCTTGCCTTGCCTTGCTCTTTGCTTTGCCTTGCCATTCCCTTGCTCCATGCCTTGCCTTTTGCCTTCACTTGCCTTGCTTTGTGCCTTGCGTTGCCTTGCCTTGCTTGGTGCCTTGACTTGCCTTGTGTCTTGAGTTCCCTTGGCTTGCCTTTTGTTTGCCTTGCCTTGTCCTTGCCTTGCCTGGCTTTGGCTTGCCTTGTGCCTTGCCTTGCCTTCCCTTTTGCTTTGACTTGCCTTGCCTTGAGCCTTGCCTTGCCTTGCCTTGCCTTGCGCCTTGCCTTGCTTGCCTTGTGCCTTGCCTTGCCTTGCCTGGTGCCTTGCCTTGTGCCTTGGCTTGCCTTGCCTTGCTTTGCCCTGCACTTTGCATTGCCTTGCTTTGCCTTGTGCCTTGCCTTGCCTTGCTTTGTGCATTGCCTTGACTTGCTTTGTGCCTTGCCTTGCCTTGCCTTGTGCCTTGCCTTGCCTTCCACCTTGTCTTGCCTTACTTTGCCTTGTGCCTTGCCTTGCTTTGCCTTGCCATGCCTTGTGCCTTGTGCCTTGCCTTGCATTGCTTTGCCTTGTCTTGTGCCTTGCCCTGCTTTGCCTTGACTTGTTTGTTGATTTGCCTTGCCTTGATACTTGCCTTCCCTTAGCTTGTGCCTTGCCTTCGCTTGTGCCTTGGGTTGCCTTGCCTTGTGCCTTGGTGTCCCTTACCTTGTGCCTTACCTTGTCTTGCCTTGCCCCTTGCATTGCCTTGCCTTGTCCCTTGCCTTGCCTTGTAACTTGCCTTGCCTTGTGCCTTGCCATGGCTTGTCTTGCTTTGCCTTTTGCCTTGTCTTACACTGTGCCTTGCCTTGCCTTGTGCCTTGCTTTGCCTTTGTTCCTTGGTTTGTTCTGTGCCTTACTTTGCCTTGCCTCGTGCCTTGCTTTGGCTTGCCTTGACTTGGGCCTTGCCTTGCTATACCTTGATTTCTGTCTTGCCTTGACTTTCCCTGCCTTGTGCCTTGCCTTGCCTTGATTTTTGCTTTGCCTTGCCATTTCCTTATGGCTTGCCCTGACTTGCCTTGCCCCTTGCCCCTTGCCTTGCTTTTTGCCTTCCCTTGCCTTGCTTTGTGCCTTGTCTTGCCTGATGCCTTTACTTGCCTTGTGCCTTGCGTTGCCTTGTGCCTTGGCTTGCCTTGCCTTGCTTTACCCTGCACTTTGCATTGCCTTGCTTTGCCTTGTGCCTTGCCTTGCTTTGCCTTGTGCCTTGCCTTGCCTTGCTTTGTGCCTTGCCTTGCCTTGCTTTGCTTTGTGCCTTGCCTTGCCTTGTGCCTTGCCTTGCCTTGTGCCATGCCTTGCCTTGTGCCGTGTCTTGCCTTGCTTTCTGCCTTGCCTTGCCTTGCCTTGCCTTCTGCCTTTTTTTGCTTTTCCTTGCTTTTCCTTGCCATGCCTTGTGCCTTGCCTTGCCTTGCTTTGCCTTGCCTTGTGCCTTGCCTTGCCTTGACTTGTGCGTTGATTTGCCTTGCCTTGTTACTTGCCTTCCCTTGGCTTGTGGCTTGTCTTGCCTTGCCTTGTGCCTTGCTTGTGCCTTGGTCTCCCTTACCTTGTGCCTTATCTTGTCTTGCCTTGTCTTGCCTTGCCTTGTATCTTGCATTGCCTTGCCTTGTGCCTTACCTTGTGCCTTCCCTTCCCTTGCCTTGCCTTGTGCCTTCCCTTGGCTTGCCTTGTAACTTGCCTTGCCTTATGCCTTGCCTTGCCTTATTCCTTGTCTTGCCATGCCTTGCCTTGCTTTGCCTTTTGCCTTGCCTTGCCTTGCATTGTGCCTTGCCTTGCCTTGCCTTGCCTTGTGCCTTGCCTTGCCTTGCCTTGCCTCGTGCCTTGCTTTGGCTTGCCTTGACTTGGGTCTTGCCTTGCCTTGCCTTGATTTTTGCTTTGCCTTGCCATTTCCTTGCAGCTTGCCTGGACTTGCCTTGTCCCTTGTCCCTTGCCTTGCCTTTTGCCTTCACTTGCCTTGCTTTGTGCTTTGCTTTGCCTGGTGCCTTGACTTGACTTGTGCCTTGCCTTCCCTTCACTTGCCTTTTGCTTGCGTTGCCTTGTGCCTTGCCTTGCCTTGGCCTTGCCTTGCCTGGCTTTGCCTTGCCTTGTGGCTTGCCTTGTCTTGCCTTGCCTTGTGCCTTGCCATGTGCTTTCCCTTGCCTTGCCTTGCCATGTGCCTTGCGTTACCTTTGCTTGTGCCTTGCCTTGCCTTCTGCCTTGCCTTGCCATGCTTTGTCTTGCCTTTTCTAGT\t*\tSA:Z:chrX,30795360,+,3650H1168M,60,0,1168;chrX,30796147,+,3549H155M1114H,60,15,80;chrX,30789491,+,3426H83M1309H,22,5,58;chrX,30794688,+,3245H100M1473H,60,10,50;\tMD:Z:123A1214G706T396G270^TTGCTTTGTGCCTTGTCTTGCCTGATGCCTTTACTTGCCTTGTGCCTTGCG34A45T4C6C8C27\tRG:Z:GATKSVContigAlignments\tNM:i:234\tAS:i:2524\tXS:i:50",
                canonicalChromosomes, refSeqDict);

        final CpxVariantInducingAssemblyContig.BasicInfo manuallyCalculatedBasicInfo =
                new CpxVariantInducingAssemblyContig.BasicInfo("chrX", true,
                        new SimpleInterval("chrX", 30792651, 30792651), new SimpleInterval("chrX", 30796527, 30796527));

        //////////
        final List<CpxVariantInducingAssemblyContig.Jump> manuallyCalculatedJumps =
                Arrays.asList(
                        new CpxVariantInducingAssemblyContig.Jump(new SimpleInterval("chrX", 30793441, 30793441), new SimpleInterval("chrX", 30793442, 30793442), StrandSwitch.NO_SWITCH, 78),
                        new CpxVariantInducingAssemblyContig.Jump(new SimpleInterval("chrX", 30794753, 30794753), new SimpleInterval("chrX", 30794754, 30794754), StrandSwitch.NO_SWITCH, 96),
                        new CpxVariantInducingAssemblyContig.Jump(new SimpleInterval("chrX", 30795363, 30795363), new SimpleInterval("chrX", 30795415, 30795415), StrandSwitch.NO_SWITCH, 0),
                        new CpxVariantInducingAssemblyContig.Jump(new SimpleInterval("chrX", 30795543, 30795543), new SimpleInterval("chrX", 30794688, 30794688), StrandSwitch.NO_SWITCH, 229),
                        new CpxVariantInducingAssemblyContig.Jump(new SimpleInterval("chrX", 30794787, 30794787), new SimpleInterval("chrX", 30789491, 30789491), StrandSwitch.NO_SWITCH, 81),
                        new CpxVariantInducingAssemblyContig.Jump(new SimpleInterval("chrX", 30789573, 30789573), new SimpleInterval("chrX", 30796147, 30796147), StrandSwitch.NO_SWITCH, 40),
                        new CpxVariantInducingAssemblyContig.Jump(new SimpleInterval("chrX", 30796247, 30796247), new SimpleInterval("chrX", 30795360, 30795360), StrandSwitch.NO_SWITCH, 0)
                );
        final List<SimpleInterval> manuallyCalculatedEventPrimaryChromosomeSegmentingLocations =
                Arrays.asList(
                        new SimpleInterval("chrX", 30793441, 30793441),
                        new SimpleInterval("chrX", 30793442, 30793442),
                        new SimpleInterval("chrX", 30794688, 30794688),
                        new SimpleInterval("chrX", 30794753, 30794753),
                        new SimpleInterval("chrX", 30794754, 30794754),
                        new SimpleInterval("chrX", 30794787, 30794787),
                        new SimpleInterval("chrX", 30795360, 30795360),
                        new SimpleInterval("chrX", 30795363, 30795363),
                        new SimpleInterval("chrX", 30795415, 30795415),
                        new SimpleInterval("chrX", 30795543, 30795543),
                        new SimpleInterval("chrX", 30796147, 30796147),
                        new SimpleInterval("chrX", 30796247, 30796247)
                );
        final CpxVariantInducingAssemblyContig manuallyCalculatedCpxVariantInducingAssemblyContig =
                new CpxVariantInducingAssemblyContig(analysisReadyContig,
                        manuallyCalculatedBasicInfo,
                        manuallyCalculatedJumps,
                        manuallyCalculatedEventPrimaryChromosomeSegmentingLocations);

        //////////
        final SimpleInterval manuallyCalculatedAffectedRefRegion =
                new SimpleInterval("chrX", 30793441, 30796247);
        final List<SimpleInterval> manuallyCalculatedSegments =
                Arrays.asList(
                        new SimpleInterval("chrX", 30793442, 30794688),
                        new SimpleInterval("chrX", 30794688, 30794753),
                        new SimpleInterval("chrX", 30794754, 30794787),
                        new SimpleInterval("chrX", 30794787, 30795360),
                        new SimpleInterval("chrX", 30795360, 30795363),
                        new SimpleInterval("chrX", 30795363, 30795415),
                        new SimpleInterval("chrX", 30795415, 30795543),
                        new SimpleInterval("chrX", 30795543, 30796147),
                        new SimpleInterval("chrX", 30796147, 30796247));
        final List<String> manuallyCalculatedAltArrangementDescription =
                Arrays.asList("UINS-78", "1", "2", "UINS-96", "3", "4", "5", "7", "UINS-229", "2", "3", "UINS-81", "chrX:30789491-30789573", "UINS-40", "9", "5", "6", "7", "8", "9");
        final byte[] manuallyCalculatedAltSeq = Arrays.copyOfRange(analysisReadyContig.getContigSequence(), 790, 4538);
        final CpxVariantCanonicalRepresentation manuallyCalculatedCpxVariantCanonicalRepresentation =
                new CpxVariantCanonicalRepresentation(
                        manuallyCalculatedAffectedRefRegion,
                        manuallyCalculatedSegments,
                        manuallyCalculatedAltArrangementDescription,
                        manuallyCalculatedAltSeq);

        //////////
        final byte[] dummyRefSequence = "TCGA".getBytes();
        final VariantContextBuilder variantContextBuilder = new VariantContextBuilder()
                .chr("chrX").start(30793441).stop(30796247)
                .alleles(Arrays.asList(Allele.create(dummyRefSequence, true), Allele.create(SimpleSVType.createBracketedSymbAlleleString(CPX_SV_SYB_ALT_ALLELE_STR), false)))
                .id("CPX_chrX:30793441-30796247")
                .attribute(VCFConstants.END_KEY, 30796247)
                .attribute(SVTYPE, GATKSVVCFConstants.CPX_SV_SYB_ALT_ALLELE_STR)
                .attribute(SVLEN, 941)
                .attribute(SEQ_ALT_HAPLOTYPE, new String(manuallyCalculatedAltSeq));
        variantContextBuilder.attribute(CPX_EVENT_ALT_ARRANGEMENTS,
                String.join(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR, manuallyCalculatedAltArrangementDescription));
        variantContextBuilder.attribute(CPX_SV_REF_SEGMENTS,
                String.join(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR, manuallyCalculatedSegments.stream().map(SimpleInterval::toString).collect(Collectors.toList())));
        final VariantContext manuallyCalculatedVariantContext = variantContextBuilder.make();

        //////////
        final Set<SimpleInterval> manuallyCalculatedTwoBaseBoundaries = new HashSet<>();
        manuallyCalculatedTwoBaseBoundaries.add(new SimpleInterval("chrX", 30792651, 30792652));
        manuallyCalculatedTwoBaseBoundaries.add(new SimpleInterval("chrX", 30793440, 30793441));
        manuallyCalculatedTwoBaseBoundaries.add(new SimpleInterval("chrX", 30793442, 30793443));
        manuallyCalculatedTwoBaseBoundaries.add(new SimpleInterval("chrX", 30794752, 30794753));
        manuallyCalculatedTwoBaseBoundaries.add(new SimpleInterval("chrX", 30794754, 30794755));
        manuallyCalculatedTwoBaseBoundaries.add(new SimpleInterval("chrX", 30795362, 30795363));
        manuallyCalculatedTwoBaseBoundaries.add(new SimpleInterval("chrX", 30795415, 30795416));
        manuallyCalculatedTwoBaseBoundaries.add(new SimpleInterval("chrX", 30795542, 30795543));
        manuallyCalculatedTwoBaseBoundaries.add(new SimpleInterval("chrX", 30794688, 30794689));
        manuallyCalculatedTwoBaseBoundaries.add(new SimpleInterval("chrX", 30794786, 30794787));
        manuallyCalculatedTwoBaseBoundaries.add(new SimpleInterval("chrX", 30796147, 30796148));
        manuallyCalculatedTwoBaseBoundaries.add(new SimpleInterval("chrX", 30796246, 30796247));
        manuallyCalculatedTwoBaseBoundaries.add(new SimpleInterval("chrX", 30795360, 30795361));
        manuallyCalculatedTwoBaseBoundaries.add(new SimpleInterval("chrX", 30796526, 30796527));

        //////////
        return new PreprocessedAndAnalysisReadyContigWithExpectedResults(
                manuallyCalculatedCpxVariantInducingAssemblyContig,
                manuallyCalculatedCpxVariantCanonicalRepresentation,
                manuallyCalculatedVariantContext,
                dummyRefSequence,
                manuallyCalculatedTwoBaseBoundaries);
    }

    private static PreprocessedAndAnalysisReadyContigWithExpectedResults buildForContig_22672_4(final Set<String> canonicalChromosomes,
                                                                                                final SAMSequenceDictionary refSeqDict) {
        final AssemblyContigWithFineTunedAlignments
                analysisReadyContig = makeContigAnalysisReady("asm022672:tig00004\t16\tchr9\t130954964\t60\t1631S192M60I99M24I54M156I1430M\t*\t0\t0\tTCATCCAGGTTGGAGTGCAGTGGTGCTATCTCAGCTCACTGCAGCCTCCACCCCTGGATTGAATCGATTCTCCTGCCTCAGCCTCCTGAGTAGCTGGGATTACAGGAGTGCACCACCACGCCCGGCTCATGTTTGTGTTTTTAGTAGAGACAGGGTTTCACCACATTGGCCAGGCTGGTCTCAAACTCCTGACCTCAAGTGATCTGGCTGCCTTGGCCTCCCAAAGTGCTGGGATTATAGGTGTGAGCCACCACACCCAGCCAGAGTGGGTAGTTTTTAAAACCACCACAATTGCCCGCCAGGGGCACAAAGAGGACTCATGGGGATGGGTCTGATAAGTGCTGAGCCCAGTGCTGGGCACCTGCAGGCACTCAGTAGGTGGTGGACACTTGTTTTATTATGATTATCCAGCACCTAGCACCCAGTCTACCCCAAATAGCACCATCAACATATCTTCCTGAATCCTGCCTCCCTCATTTCAGCCTTTGCTCTCAGGCAGCCAGGGGTTCCCCATTTTCAACAAGAAAAAGCCCAAACCCTTTTGCTTGGCAGTCAACCCCACCCCATCCAAACCAATCCCTCCTTTCTCTGCTGAATGTTTCCCCTGTGCCAGCCACTTCCTTCCCACCACCACACCCTTAGCCAGGCTATCACCCACCAAGAATCCTGCCCCATCTTGAGCCTCTGCCTTCATTTGAAGCCTATCATCCCCTGTGCCTCAGTATTTTCCTTCAGCTCTAGATTCTTGCAGTGCTCTAGCCATTTACACAGCCATCAAATTCCAGGCCTCTTCTGTCACTCACCGGTTCACCTCTGTCATCTTGTCTCCTGCCCAAGACTGCATCTCCTCTTGTCCCTTCTCTGAGAGGCCTGTGGTTGAGCACAGGGATGAATGAGAAAGAAACCCAGCCTGTCCTCAAGAAGTTGACAATGTGTCCCTCAAAGACTGGGCTCATTGCTGCTACCTCTGGGAGTCCCAATGTGGAGGAGAGCTCTCAGTGGGTGTGCAGAAAATCTTGTCAGAAGTTACTGTCCCTGGGCAGGGCTAATACCACCAGTATCACCATTATCATCATCACCACTATCGTCATCATCACCACCATCACCATCATTACCACCATCATCACCATCAATATCGACACCATCCTCATCATCATCATCACCACCATCATCAGCATCATCATCAGCAGCAGCACCACCATTACCATCATTATCCCTACCATCTTCATCACCACCATCACCATCACCATCACCATCATCATCACCACCATCACCATCATCACCACCACCATCATCATCACCACCATCACCACCATCATCACCATCATCATCACCATCATCACCACCATCATCACCATCATCATCATCAGCACCACCACCACTATTACCATCATTATCACTACCATCCTCATCACCACCATCACTATCACCATCACCATCATCATCACCATCATCACCATCATCATCACCATCACCATCATCACCATCACCATCACCATCATCACCATCATCATCATCACCATCACCATCATCATCACCATCATCACCACCATCATCACCATCATCATCATCAGCACCACCACCACTATTACCATCATTATCCCTACCATCTTCATCACCACCATCACCATCACCATCACCATCATCATCACCACCATCACCATCACCATCACCATCATCATCACCATCATCACCACCATCATCACCATCATCATCATCAGCACCACCACCACTATTACCATCATTATCCCTACCATCCTCATCACCACCATCACCATCACCATCACCATCATCATCACCATCACCATCATCACCATCACCATCATCAACATCACCATCACCATCATCACCATCATCATCATCACCATCACCATCATCATCACCATCATCACCACCATCATCACCATCATCAGCAGCAGCAGCACCACCACCACTATTACCATCATTATCCCTACCATCCTCATCACCACCATCACCATCACCATCATCACCATCACCATCACCATCATCACCATCACCATCATCATCATCACCACCATCATCACCATCATCATCAGCAGCAGCACCACCACCACTATTACCATCATTATTCCTACCATCCTCATCACCACCCTCACCATCACCATCATCATCATGACCATCATCATCAGCACCACCGTTATCATCATTATCCCCACCGTCCTCATCACCATCATCACCATCACCATCATCATTGTCACCATCACCATGATCATCGTCACCATCTTTATCCCCACCATCCTCATCACCATCATCACTACCATCAGCACCATCACCATCATCACCATTACCACCATTATCACCACCATCATCATCACCAACACCATCCTCAGCACTACCACCACCATCACATTCTCATTCATCCCTGTGCTCAACCACAGGCCTCTCAGAGAAGGGACAAGAGGAGATGCAGTCTTGGGCAGGAGACCACTTTTGGGGGCAAAAGCAGCTCCTGAAACCACACCCCAAAGCAGGCTCCACTCCATTCCATCAGACCCTGCAGTCAGGAAGGGCCGTGGGGTGTCCTGGCCTTCATCTGTGAGAACTGCCTTACCCATGCTGATTTCCACCCACATGGCATCAGAGGACTCCGTGCCCTGACACTACAATCCTTGTCCCCTTGTGTAGCCACCTCAGGGTCCAGCACCAACTGTCCTGTCTACCTGGAGCATAACACCATGAGGTCCCCAGCTCCTGGCTCTCCCCTGGGGGCTTGCCATGACCCGTTGGCCTGAGCTTTGGGGTGTCAGAAGCCCCCATGGAATGCACTGCTAAGACAAGCCCATCCTGCCAGCCTTCCTACCTCTCATAGCTGACCCTGCTGGGTTCTATGTTACAAATTCCCCTGCCCCAGCACCACTCTGTGACCTGCAGTGAGCGCAGCTGTCACTCCCCTCAAGGGTGCCCAGGCAATAGGGAGATGTGAGGACTGTGGGGGCATAAGAGTTGTCTCAGGGCTGTGCAGAGGAAGCCAGGGCCCCCCTAAATATCTCCAGCTCATCGCACCAGCTCCATTCTGCCCGGAGCCCCATGTCTGAGCCCTCAGCTCAGAGCAGACAACCGTGCTAATCCTGTCCCAGGCTGGTGGCCGCCCCCACTGAGGAGGGGGAGCAAGCTGCCCAGAAAGCTGGTGTGGAGTGCCCGTTAGCTGTGGATGGCATCGTGGTGCCACCGGGAGATTATCCACACGCAGCTAGTTCCCAGCAGCCGACCCCAGTGCCCCAAAAGAAGGCAGCTAGCTATTAAGGAGATTCATGGGCAGGGGTGAGGCGAGGAGGAAGGCTAATCAAGCTGTCCTGATTGTAATTGTCCGAAGGTGCTGGCCTCTCTCGTAAACATGCCAACTGCAGGCCCTGCTTCTGCGTCTCAGCAGGACTTCTCCCTGCTGGGACCGTGCTGGGGAATGTTCTTTCAACATAACTCTATTTAAATTCACATTTCCATCATCCCCCAGAGGAGCCGAGAGAGCTGAGCTGCGTGGTATAAGCCGGGAAAGGATTAGATGGGGTGGGTGTTATTTTTTTTCCTTTTATTTTCCCTTAGTGATGGAGATGGGGTTGTGGGGGGTGGGCCATTCTAGAATTCTGCAGTATTGGAACTGGAAGAGCTGTTACAAACCATCCAGCTC\t*\tSA:Z:chr9,130953867,-,1274M42I167M2163H,60,55,1318;chr9,130955093,-,1469H33M3I179M1962H,19,14,138;\tMD:Z:15T1C5T11A0T3G0A2C4C3T8T5C5T2G8G5G2G6C24T45T2C23T5C14T46T33T5C32T1433\tRG:Z:GATKSVContigAlignments\tNM:i:268\tAS:i:1347\tXS:i:176",
                canonicalChromosomes, refSeqDict);

        final CpxVariantInducingAssemblyContig.BasicInfo manuallyCalculatedBasicInfo =
                new CpxVariantInducingAssemblyContig.BasicInfo("chr9", false,
                        new SimpleInterval("chr9", 130953867, 130953867), new SimpleInterval("chr9", 130956738, 130956738));

        //////////
        final List<CpxVariantInducingAssemblyContig.Jump> manuallyCalculatedJumps =
                Arrays.asList(
                        new CpxVariantInducingAssemblyContig.Jump(new SimpleInterval("chr9", 130955309, 130955309), new SimpleInterval("chr9", 130955308, 130955308), StrandSwitch.NO_SWITCH, 156),
                        new CpxVariantInducingAssemblyContig.Jump(new SimpleInterval("chr9", 130955156, 130955156), new SimpleInterval("chr9", 130955155, 130955155), StrandSwitch.NO_SWITCH, 60),
                        new CpxVariantInducingAssemblyContig.Jump(new SimpleInterval("chr9", 130954964, 130954964), new SimpleInterval("chr9", 130955307, 130955307), StrandSwitch.NO_SWITCH, 148)
                );
        final List<SimpleInterval> manuallyCalculatedEventPrimaryChromosomeSegmentingLocations =
                Arrays.asList(
                        new SimpleInterval("chr9", 130954964, 130954964),
                        new SimpleInterval("chr9", 130955155, 130955155),
                        new SimpleInterval("chr9", 130955156, 130955156),
                        new SimpleInterval("chr9", 130955307, 130955307),
                        new SimpleInterval("chr9", 130955308, 130955308),
                        new SimpleInterval("chr9", 130955309, 130955309)
                );
        final CpxVariantInducingAssemblyContig manuallyCalculatedCpxVariantInducingAssemblyContig =
                new CpxVariantInducingAssemblyContig(analysisReadyContig,
                        manuallyCalculatedBasicInfo,
                        manuallyCalculatedJumps,
                        manuallyCalculatedEventPrimaryChromosomeSegmentingLocations);

        //////////
        final SimpleInterval manuallyCalculatedAffectedRefRegion =
                new SimpleInterval("chr9", 130954964, 130955309);
        final List<SimpleInterval> manuallyCalculatedSegments =
                Arrays.asList(
                        new SimpleInterval("chr9", 130954964, 130955155),
                        new SimpleInterval("chr9", 130955156, 130955307),
                        new SimpleInterval("chr9", 130955307, 130955308));
        final List<String> manuallyCalculatedAltArrangementDescription =
                Arrays.asList("1", "2", "UINS-148", "1", "UINS-60", "2", "3", "UINS-156");
        final byte[] manuallyCalculatedAltSeq = Arrays.copyOfRange(analysisReadyContig.getContigSequence(), 1429, 2549);
        SequenceUtil.reverseComplement(manuallyCalculatedAltSeq);
        final CpxVariantCanonicalRepresentation manuallyCalculatedCpxVariantCanonicalRepresentation =
                new CpxVariantCanonicalRepresentation(
                        manuallyCalculatedAffectedRefRegion,
                        manuallyCalculatedSegments,
                        manuallyCalculatedAltArrangementDescription,
                        manuallyCalculatedAltSeq);

        //////////
        final byte[] dummyRefSequence = "TCGA".getBytes();
        final VariantContextBuilder variantContextBuilder = new VariantContextBuilder()
                .chr("chr9").start(130954964).stop(130955309)
                .alleles(Arrays.asList(Allele.create(dummyRefSequence, true), Allele.create(SimpleSVType.createBracketedSymbAlleleString(CPX_SV_SYB_ALT_ALLELE_STR), false)))
                .id("CPX_chr9:130954964-130955309")
                .attribute(VCFConstants.END_KEY, 130955309)
                .attribute(SVTYPE, GATKSVVCFConstants.CPX_SV_SYB_ALT_ALLELE_STR)
                .attribute(SVLEN, 774)
                .attribute(SEQ_ALT_HAPLOTYPE, new String(manuallyCalculatedAltSeq));
        variantContextBuilder.attribute(CPX_EVENT_ALT_ARRANGEMENTS,
                String.join(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR, manuallyCalculatedAltArrangementDescription));
        variantContextBuilder.attribute(CPX_SV_REF_SEGMENTS,
                String.join(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR, manuallyCalculatedSegments.stream().map(SimpleInterval::toString).collect(Collectors.toList())));
        final VariantContext manuallyCalculatedVariantContext = variantContextBuilder.make();

        //////////
        final Set<SimpleInterval> manuallyCalculatedTwoBaseBoundaries = new HashSet<>();
        manuallyCalculatedTwoBaseBoundaries.add(new SimpleInterval("chr9", 130955309, 130955310));
        manuallyCalculatedTwoBaseBoundaries.add(new SimpleInterval("chr9", 130956737, 130956738));
        manuallyCalculatedTwoBaseBoundaries.add(new SimpleInterval("chr9", 130955156, 130955157));
        manuallyCalculatedTwoBaseBoundaries.add(new SimpleInterval("chr9", 130955307, 130955308));
        manuallyCalculatedTwoBaseBoundaries.add(new SimpleInterval("chr9", 130954964, 130954965));
        manuallyCalculatedTwoBaseBoundaries.add(new SimpleInterval("chr9", 130955154, 130955155));
        manuallyCalculatedTwoBaseBoundaries.add(new SimpleInterval("chr9", 130953867, 130953868));
        manuallyCalculatedTwoBaseBoundaries.add(new SimpleInterval("chr9", 130955306, 130955307));

        //////////
        return new PreprocessedAndAnalysisReadyContigWithExpectedResults(
                manuallyCalculatedCpxVariantInducingAssemblyContig,
                manuallyCalculatedCpxVariantCanonicalRepresentation,
                manuallyCalculatedVariantContext,
                dummyRefSequence,
                manuallyCalculatedTwoBaseBoundaries);
    }

    private static AssemblyContigWithFineTunedAlignments makeContigAnalysisReady(final String primarySAMRecord,
                                                                                 final Set<String> canonicalChromosomes,
                                                                                 final SAMSequenceDictionary refSeqDict) {
        final AlignedContig alignedContig =
                SVTestUtils.fromPrimarySAMRecordString(primarySAMRecord, true);
        final AssemblyContigWithFineTunedAlignments intermediate =
                AssemblyContigAlignmentsConfigPicker.reConstructContigFromPickedConfiguration(
                        new Tuple2<>(new Tuple2<>(alignedContig.getContigName(), alignedContig.getContigSequence()),
                                AssemblyContigAlignmentsConfigPicker.pickBestConfigurations(alignedContig, canonicalChromosomes,
                                                                            0.0)))
                        .next();

        final AssemblyContigAlignmentsConfigPicker.GoodAndBadMappings goodAndBadMappings =
                AssemblyContigAlignmentsConfigPicker.removeNonUniqueMappings(
                        intermediate.getAlignments(),
                        AssemblyContigAlignmentsConfigPicker.LAST_ROUND_TUNING_ALIGNMENT_MAPQUAL_THREHOLD,
                        AssemblyContigAlignmentsConfigPicker.LAST_ROUND_TUNING_ALIGNMENT_READSPAN_THRESHOLD);
        final AssemblyContigWithFineTunedAlignments result = new AssemblyContigWithFineTunedAlignments(
                new AlignedContig(intermediate.getContigName(), intermediate.getContigSequence(), goodAndBadMappings.getGoodMappings()),
                goodAndBadMappings.getBadMappingsAsCompactStrings(), false, (AlignmentInterval) null);
        return CpxVariantInterpreter.furtherPreprocess(result, refSeqDict);
    }
}
