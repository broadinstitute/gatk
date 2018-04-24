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
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignedContig;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AssemblyContigAlignmentsConfigPicker;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AssemblyContigWithFineTunedAlignments;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.StrandSwitch;
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

        private PreprocessedAndAnalysisReadyContigWithExpectedResults(final CpxVariantInducingAssemblyContig expectedCpxVariantInducingAssemblyContig,
                                                                      final CpxVariantCanonicalRepresentation expectedCpxVariantCanonicalRepresentation,
                                                                      final VariantContext expectedVariantContext,
                                                                      final byte[] assumedReferenceSequence) {
            this.expectedCpxVariantInducingAssemblyContig = expectedCpxVariantInducingAssemblyContig;
            this.expectedCpxVariantCanonicalRepresentation = expectedCpxVariantCanonicalRepresentation;
            this.expectedVariantContext = expectedVariantContext;
            this.assumedReferenceSequence = assumedReferenceSequence;
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
                        new CpxVariantInducingAssemblyContig.Jump(new SimpleInterval("chr2", 4452399, 4452399), new SimpleInterval("chr9", 34840477, 34840477), StrandSwitch.NO_SWITCH, 117),
                        new CpxVariantInducingAssemblyContig.Jump(new SimpleInterval("chr9", 34840411, 34840411), new SimpleInterval("chr6", 15853414, 15853414), StrandSwitch.NO_SWITCH, 114),
                        new CpxVariantInducingAssemblyContig.Jump(new SimpleInterval("chr6", 15853372, 15853372), new SimpleInterval("chr2", 4452357, 4452357), StrandSwitch.NO_SWITCH, 53),
                        new CpxVariantInducingAssemblyContig.Jump(new SimpleInterval("chr2", 4452298, 4452298), new SimpleInterval("chr2", 4452406, 4452406), StrandSwitch.NO_SWITCH, 317)
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
                Arrays.asList("1", "2", "3", "UINS-317", "1", "UINS-53", "chr6:15853372-15853414", "UINS-114", "chr9:34840411-34840477", "UINS-117", "3");
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
                .attribute(SVLEN, 109)
                .attribute(SEQ_ALT_HAPLOTYPE, new String(manuallyCalculatedAltSeq));
        variantContextBuilder.attribute(CPX_EVENT_ALT_ARRANGEMENTS,
                String.join(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR, manuallyCalculatedAltArrangementDescription));
        variantContextBuilder.attribute(CPX_SV_REF_SEGMENTS,
                String.join(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR, manuallyCalculatedSegments.stream().map(SimpleInterval::toString).collect(Collectors.toList())));
        final VariantContext manuallyCalculatedVariantContext = variantContextBuilder.make();

        //////////
        return new PreprocessedAndAnalysisReadyContigWithExpectedResults(
                manuallyCalculatedCpxVariantInducingAssemblyContig,
                manuallyCalculatedCpxVariantCanonicalRepresentation,
                manuallyCalculatedVariantContext,
                dummyRefSequence);
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
                        new CpxVariantInducingAssemblyContig.Jump(new SimpleInterval("chr2", 16225125, 16225125), new SimpleInterval("chr2", 16226720, 16226720), StrandSwitch.NO_SWITCH, 6),
                        new CpxVariantInducingAssemblyContig.Jump(new SimpleInterval("chr2", 16226497, 16226497), new SimpleInterval("chr2", 16225125, 16225125), StrandSwitch.REVERSE_TO_FORWARD, 27),
                        new CpxVariantInducingAssemblyContig.Jump(new SimpleInterval("chr2", 16225237, 16225237), new SimpleInterval("chr2", 16225142, 16225142), StrandSwitch.FORWARD_TO_REVERSE, 1)
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
                Arrays.asList("1", "UINS-1", "-2", "-1", "UINS-27", "chr2:16226497-16226720", "UINS-6", "1", "2");
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
                .attribute(SVLEN, 113)
                .attribute(SEQ_ALT_HAPLOTYPE, new String(manuallyCalculatedAltSeq));
        variantContextBuilder.attribute(CPX_EVENT_ALT_ARRANGEMENTS,
                String.join(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR, manuallyCalculatedAltArrangementDescription));
        variantContextBuilder.attribute(CPX_SV_REF_SEGMENTS,
                String.join(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR, manuallyCalculatedSegments.stream().map(SimpleInterval::toString).collect(Collectors.toList())));
        final VariantContext manuallyCalculatedVariantContext = variantContextBuilder.make();

        //////////
        return new PreprocessedAndAnalysisReadyContigWithExpectedResults(
                manuallyCalculatedCpxVariantInducingAssemblyContig,
                manuallyCalculatedCpxVariantCanonicalRepresentation,
                manuallyCalculatedVariantContext,
                dummyRefSequence);
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
                        new CpxVariantInducingAssemblyContig.Jump(new SimpleInterval("chrX", 30793441, 30793441), new SimpleInterval("chrX", 30793442, 30793442), StrandSwitch.NO_SWITCH, 77),
                        new CpxVariantInducingAssemblyContig.Jump(new SimpleInterval("chrX", 30794753, 30794753), new SimpleInterval("chrX", 30794754, 30794754), StrandSwitch.NO_SWITCH, 95),
                        new CpxVariantInducingAssemblyContig.Jump(new SimpleInterval("chrX", 30795363, 30795363), new SimpleInterval("chrX", 30795415, 30795415), StrandSwitch.NO_SWITCH, 0),
                        new CpxVariantInducingAssemblyContig.Jump(new SimpleInterval("chrX", 30795543, 30795543), new SimpleInterval("chrX", 30794688, 30794688), StrandSwitch.NO_SWITCH, 228),
                        new CpxVariantInducingAssemblyContig.Jump(new SimpleInterval("chrX", 30794787, 30794787), new SimpleInterval("chrX", 30789491, 30789491), StrandSwitch.NO_SWITCH, 80),
                        new CpxVariantInducingAssemblyContig.Jump(new SimpleInterval("chrX", 30789573, 30789573), new SimpleInterval("chrX", 30796147, 30796147), StrandSwitch.NO_SWITCH, 39),
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
                Arrays.asList("UINS-77", "1", "2", "UINS-95", "3", "4", "5", "7", "UINS-228", "2", "3", "UINS-80", "chrX:30789491-30789573", "UINS-39", "9", "5", "6", "7", "8", "9");
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
                .attribute(SVLEN, 2807)
                .attribute(SEQ_ALT_HAPLOTYPE, new String(manuallyCalculatedAltSeq));
        variantContextBuilder.attribute(CPX_EVENT_ALT_ARRANGEMENTS,
                String.join(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR, manuallyCalculatedAltArrangementDescription));
        variantContextBuilder.attribute(CPX_SV_REF_SEGMENTS,
                String.join(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR, manuallyCalculatedSegments.stream().map(SimpleInterval::toString).collect(Collectors.toList())));
        final VariantContext manuallyCalculatedVariantContext = variantContextBuilder.make();

        //////////
        return new PreprocessedAndAnalysisReadyContigWithExpectedResults(
                manuallyCalculatedCpxVariantInducingAssemblyContig,
                manuallyCalculatedCpxVariantCanonicalRepresentation,
                manuallyCalculatedVariantContext,
                dummyRefSequence);
    }


    private static AssemblyContigWithFineTunedAlignments makeContigAnalysisReady(final String primarySAMRecord,
                                                                                 final Set<String> canonicalChromosomes,
                                                                                 final SAMSequenceDictionary refSeqDict) {
        final AlignedContig alignedContig =
                SVTestUtils.fromPrimarySAMRecordString(primarySAMRecord, true);
        final AssemblyContigWithFineTunedAlignments toBeDeOverlapped =
                AssemblyContigAlignmentsConfigPicker.reConstructContigFromPickedConfiguration(
                        new Tuple2<>(new Tuple2<>(alignedContig.getContigName(), alignedContig.getContigSequence()),
                                AssemblyContigAlignmentsConfigPicker.pickBestConfigurations(alignedContig, canonicalChromosomes,
                                                                            0.0)))
                        .next();
        return CpxVariantInterpreter.furtherPreprocess(toBeDeOverlapped, refSeqDict);
    }
}
