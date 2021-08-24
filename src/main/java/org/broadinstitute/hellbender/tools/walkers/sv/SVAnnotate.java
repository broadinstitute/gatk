package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.samtools.util.OverlapDetector;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.codecs.gtf.*;

import java.io.File;
import java.util.*;

/**
 * Adds gene overlap and variant consequence annotations to SV VCF from GATK-SV pipeline
 * Input files are an SV VCF and a GTF file containing primary or canonical transcripts
 * Output file is an annotated SV VCF
 */
@CommandLineProgramProperties(
        summary = "Adds gene overlap and variant consequence annotations to SV VCF from GATK-SV pipeline." +
                "Input files are an SV VCF and a GTF file containing primary or canonical transcripts." +
                "Output file is an annotated SV VCF.",
        oneLineSummary = "Adds gene overlap and variant consequence annotations to SV VCF from GATK-SV pipeline",
        programGroup = StructuralVariantDiscoveryProgramGroup.class
)
public final class SVAnnotate extends VariantWalker {

    @Argument(
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "Output file (if not provided, defaults to STDOUT)",
            common = false,
            optional = true
    )
    private GATKPath outputFile = null;

    @Argument(
            fullName="proteinCodingGTF",
            shortName="P",
            doc="protein-coding GTF file (canonical only)",
            optional=true
    )
    private File proteinCodingGTFFile;

    private VariantContextWriter vcfWriter = null;

    private OverlapDetector<GencodeGtfGeneFeature> gtfOverlapDetector;

    private final Set<String> MSVExonOverlapClassifications = new HashSet<>(Arrays.asList(GATKSVVCFConstants.LOF, GATKSVVCFConstants.INT_EXON_DUP, GATKSVVCFConstants.DUP_PARTIAL, GATKSVVCFConstants.PARTIAL_EXON_DUP, GATKSVVCFConstants.COPY_GAIN));

    @Override
    public void onTraversalStart() {
        FeatureDataSource<GencodeGtfGeneFeature> proteinCodingGTFSource = new FeatureDataSource<>(proteinCodingGTFFile);
        List<GencodeGtfGeneFeature> gtfFeaturesList = new ArrayList<>();
        proteinCodingGTFSource.forEach(gtfFeaturesList::add);  // TODO: faster method?
        gtfOverlapDetector = OverlapDetector.create(gtfFeaturesList);

        vcfWriter = createVCFWriter(outputFile);
        updateAndWriteHeader();
    }

    private void addAnnotationInfoKeysToHeader(VCFHeader header) {
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.LOF, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Gene(s) (comma-separated) on which the SV is predicted to have a loss-of-function effect."));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.INT_EXON_DUP, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Gene(s) (comma-separated) on which the SV is predicted to result in intragenic exonic duplication."));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.COPY_GAIN, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Gene(s) (comma-separated) on which the SV is predicted to have a copy-gain effect."));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.DUP_PARTIAL, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Gene(s) (comma-separated) which are partially overlapped by an SV's duplication."));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.INTRONIC, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Gene(s) (comma-separated) where the SV was found to lie entirely within an intron."));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.PARTIAL_EXON_DUP, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Gene(s) (comma-separated) where the SV was found to partially overlap a single exon."));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.INV_SPAN, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Gene(s) (comma-separated) which are entirely spanned by an SV's inversion."));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.UTR, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Gene(s) (comma-separated) for which the SV is predicted to disrupt a UTR."));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.MSV_EXON_OVERLAP, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Gene(s) (comma-separated) on which the multiallelic SV would be predicted to have a LOF, INTRAGENIC_EXON_DUP, COPY_GAIN, DUP_PARTIAL, or PARTIAL_EXON_DUP annotation if the SV were biallelic."));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.PROMOTER, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Gene(s) (comma-separated) for which the SV is predicted to overlap its promoter."));

    }

    private void updateAndWriteHeader() {
        final VCFHeader header = getHeaderForVariants();
        addAnnotationInfoKeysToHeader(header);
        vcfWriter.writeHeader(header);
    }

    protected final static boolean variantSpansFeature(SimpleInterval variantInterval, SimpleInterval featureInterval) {
        if (!variantInterval.getContig().equals(featureInterval.getContig())) {
            return false;
        } else {
            return variantInterval.getStart() <= featureInterval.getStart() && variantInterval.getEnd() >= featureInterval.getEnd();
        }
    }

    protected final static int countBreakpointsInsideFeature(SimpleInterval variantInterval, SimpleInterval featureInterval) {
        int count = 0;
        if (variantInterval.getContig().equals(featureInterval.getContig())) {
            if (variantInterval.getStart() >= featureInterval.getStart() && variantInterval.getStart() <= featureInterval.getEnd()) {
                count++;
            }
            if (variantInterval.getEnd() >= featureInterval.getStart() && variantInterval.getEnd() <= featureInterval.getEnd()) {
                count++;
            }
        }
        return count;
    }

    protected final static boolean variantOverlapsFeature(SimpleInterval variantInterval, SimpleInterval featureInterval) {
        return IntervalUtils.overlaps(variantInterval, featureInterval);
    }

    private void updateVariantConsequenceDict(Map<String, Set<String>> variantConsequenceDict, String key, String value) {
        if (!variantConsequenceDict.containsKey(key)) {
            variantConsequenceDict.put(key, new HashSet<>());
        }
        variantConsequenceDict.get(key).add(value);
    }

    private String annotateDEL(SimpleInterval variantInterval, GencodeGtfTranscriptFeature gtfTranscript) {
        List<GencodeGtfFeature> gtfFeaturesForTranscript = gtfTranscript.getAllFeatures();
        String consequence = GATKSVVCFConstants.INTRONIC;
        for (GencodeGtfFeature gtfFeature : gtfFeaturesForTranscript) {
            SimpleInterval featureInterval = new SimpleInterval(gtfFeature);
            if (!variantOverlapsFeature(variantInterval, featureInterval)) {
                continue;
            }
            if (gtfFeature.getFeatureType() == GencodeGtfFeature.FeatureType.CDS) {
                consequence = GATKSVVCFConstants.LOF;
                break;
            } else if (gtfFeature.getFeatureType() == GencodeGtfFeature.FeatureType.UTR) {
                consequence = GATKSVVCFConstants.UTR;
            }
        }
        return consequence;
    }

    private String annotateINS(SimpleInterval variantInterval, GencodeGtfTranscriptFeature gtfTranscript) {
        return annotateDEL(variantInterval, gtfTranscript);
    }

    private String annotateDUP(SimpleInterval variantInterval, GencodeGtfTranscriptFeature gtfTranscript) {
        String consequence = GATKSVVCFConstants.INTRONIC;
        SimpleInterval transcriptInterval = new SimpleInterval(gtfTranscript);
        if (variantSpansFeature(variantInterval, transcriptInterval)) {
            consequence = GATKSVVCFConstants.COPY_GAIN;
        } else if (countBreakpointsInsideFeature(variantInterval, transcriptInterval) == 1) {
            consequence = GATKSVVCFConstants.DUP_PARTIAL;
        } else {
            // both breakpoints inside transcript
            List<GencodeGtfFeature> gtfFeaturesForTranscript = gtfTranscript.getAllFeatures();
            int numBreakpointsInExon = 0;  // TODO: CDS or exon?
            int numBreakpointsInUTR = 0;
            int numExonsSpanned = 0;
            for (GencodeGtfFeature gtfFeature : gtfFeaturesForTranscript) {
                SimpleInterval featureInterval = new SimpleInterval(gtfFeature);
                if (!variantOverlapsFeature(variantInterval, featureInterval)) {
                    continue;
                }
                if (gtfFeature.getFeatureType() == GencodeGtfFeature.FeatureType.EXON) {
                    if (variantSpansFeature(variantInterval, featureInterval)) {
                        numExonsSpanned++;  // TODO: CDS or exon? may differ from breakpoints
                    } else {
                        numBreakpointsInExon = numBreakpointsInExon + countBreakpointsInsideFeature(variantInterval, featureInterval);
                    }
                } else if (gtfFeature.getFeatureType() == GencodeGtfFeature.FeatureType.UTR) {
                    numBreakpointsInUTR = numBreakpointsInUTR + countBreakpointsInsideFeature(variantInterval, featureInterval);
                }
            }
            if (numBreakpointsInExon == 2) {
                consequence = GATKSVVCFConstants.LOF;
            } else if (numExonsSpanned > 0) {
                consequence = GATKSVVCFConstants.INT_EXON_DUP;  // formerly DUP_LOF - consider INTERNAL
            } else if (numBreakpointsInExon == 1) {
                consequence = GATKSVVCFConstants.PARTIAL_EXON_DUP;  // new category - could collapse with DUP_PARTIAL
            } else if (numBreakpointsInUTR > 0) {
                consequence = GATKSVVCFConstants.UTR;
            }
        }
        return consequence;
    }

    private String annotateCNV(SimpleInterval variantInterval, GencodeGtfTranscriptFeature gtfTranscript) {
        String consequence = annotateDUP(variantInterval, gtfTranscript);
        if (MSVExonOverlapClassifications.contains(consequence)) {
            return GATKSVVCFConstants.MSV_EXON_OVERLAP;
        } else {
            return consequence;  // TODO: MCNV classifications ???
        }
    }

    private String annotateINV(SimpleInterval variantInterval, GencodeGtfTranscriptFeature gtfTranscript) {
        String consequence = GATKSVVCFConstants.INTRONIC;
        SimpleInterval transcriptInterval = new SimpleInterval(gtfTranscript);
        if (variantSpansFeature(variantInterval, transcriptInterval)) {
            consequence = GATKSVVCFConstants.INV_SPAN;
        } else if (countBreakpointsInsideFeature(variantInterval, transcriptInterval) == 1) {
            consequence = GATKSVVCFConstants.LOF;
        } else {
            // both breakpoints inside transcript
            List<GencodeGtfFeature> gtfFeaturesForTranscript = gtfTranscript.getAllFeatures();
            for (GencodeGtfFeature gtfFeature : gtfFeaturesForTranscript) {
                SimpleInterval featureInterval = new SimpleInterval(gtfFeature);
                if (!variantOverlapsFeature(variantInterval, featureInterval)) {
                    continue;
                }
                // TODO: if overlaps exon, it's LOF unless both breakpoints are in the same UTR?
                if (gtfFeature.getFeatureType() == GencodeGtfFeature.FeatureType.EXON) {
                    consequence = GATKSVVCFConstants.LOF;  // TODO: CDS or exon here?
                } else if (gtfFeature.getFeatureType() == GencodeGtfFeature.FeatureType.UTR) {
                    if (countBreakpointsInsideFeature(variantInterval, featureInterval) == 2) {
                        consequence = GATKSVVCFConstants.UTR;
                    }
                }
            }
        }
        return consequence;
    }

    private void annotateTranscript(SimpleInterval variantInterval, String svType, GencodeGtfTranscriptFeature transcript, Map<String, Set<String>> variantConsequenceDict) {
        if (!variantOverlapsFeature(variantInterval, new SimpleInterval(transcript))) {
            return;
        }
        String consequence = null;
        if (svType.equals("DEL")) {
            consequence = annotateDEL(variantInterval, transcript);
        } else if (svType.contains("INS")) {
            consequence = annotateINS(variantInterval, transcript);
        } else if (svType.equals("DUP")) {
            consequence = annotateDUP(variantInterval, transcript);
        } else if (svType.equals("CNV")) {
            consequence = annotateCNV(variantInterval,transcript);
        } else if (svType.equals("INV")) {
            consequence = annotateINV(variantInterval, transcript);
        }

        if (consequence != null) {
            updateVariantConsequenceDict(variantConsequenceDict, consequence, transcript.getGeneName());
        }
    }

    private String getSVType(VariantContext variant) {
        // TODO: haha majorly clean this up
        // return variant.getStructuralVariantType().name();
        Allele alt = variant.getAlternateAllele(0); // TODO: any chance of multiallelic alt field for SV?
        if (alt.isBreakpoint()) {
            return "BND";
        } else if (alt.isSingleBreakend()) {
            System.out.println("what even is single breakend???");
            System.out.println(alt.toString());
            return alt.toString();
        } else if (alt.isSymbolic()) {
            return alt.toString().substring(1, alt.toString().length()-1);  // assume <SVTYPE>
        } else {
            System.out.println("unexpected alt allele type");
            System.out.println(alt.toString());
            return alt.toString();
        }
    }

    private void annotateSVInterval(SimpleInterval variantInterval, String svType, Map<String, Set<String>> variantConsequenceDict) {
        Set<GencodeGtfGeneFeature> gtfGenesForVariant = gtfOverlapDetector.getOverlaps(variantInterval);
        for (GencodeGtfGeneFeature geneOverlapped : gtfGenesForVariant) {
            List<GencodeGtfTranscriptFeature> transcriptsForGene = geneOverlapped.getTranscripts();
            for (GencodeGtfTranscriptFeature transcript : transcriptsForGene) {
                annotateTranscript(variantInterval, svType, transcript, variantConsequenceDict);
            }
        }
    }

    private List<Pair<String, SimpleInterval>> getCPXIntervals(VariantContext variant) {
        // TODO: assert CPX only? checked before calling function so maybe unnecessary
        List<Pair<String, SimpleInterval>> intervals = new ArrayList<>();
        List<String> cpxIntervalsString = variant.getAttributeAsStringList("CPX_INTERVALS", "NONE");
        for (String cpxInterval : cpxIntervalsString) {
            String[] parsed = cpxInterval.split("_");
            String svTypeForInterval = parsed[0];
            SimpleInterval interval = new SimpleInterval(parsed[1]);
            intervals.add(Pair.of(svTypeForInterval, interval));
        }
        return intervals;
    }

    @Override
    public void apply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        if (!variant.getContig().equals("chr19")) {
            return;
        }
        Map<String, Set<String>> variantConsequenceDict = new HashMap<>();
        String overallSVType = getSVType(variant);
        if (overallSVType.equals("CPX")) {
            List<Pair<String, SimpleInterval>> svIntervals = getCPXIntervals(variant);
            for (Pair<String, SimpleInterval> typeIntervalPair : svIntervals) {
                annotateSVInterval(typeIntervalPair.getRight(), typeIntervalPair.getLeft(), variantConsequenceDict);
            }
        } else {
            annotateSVInterval(new SimpleInterval(variant), overallSVType, variantConsequenceDict);
        }

        if (overallSVType.equals("CPX")) {
            System.out.println(variant);
            System.out.println(variantConsequenceDict);
        }

        // if !variantConsequenceDict.isEmpty() then write annotations

    }

    @Override
    public void closeTool() {
        if ( vcfWriter != null ) {
            vcfWriter.close();
        }
    }
}
