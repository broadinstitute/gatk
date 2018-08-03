package org.broadinstitute.hellbender.tools.walkers.variantutils;

import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import picard.cmdline.programgroups.VariantManipulationProgramGroup;

import java.io.File;
import java.io.PrintWriter;
import java.util.*;
import java.util.stream.Collectors;

@CommandLineProgramProperties(
        summary = "This tool applies the phasing from one VCF to the genotypes in a second, producing reports on cases where genotypes are incompatible or sites are missing",
        oneLineSummary = "Apply Phasing to VCF",
        programGroup = VariantManipulationProgramGroup.class
)
@DocumentedFeature
public class ApplyPhasing extends VariantWalker {

    public static final String DISCORDANT_GENTOYPES_INFO_KEY = "PhasedDiscordantGTs";
    public static final String CONCORDANT_GENOTYPE_FILTER_STRINGS_INFO_KEY = "PhasedConcordantFiltered";
    public static final String NOPHASEMATCH_INFO_KEY = "NOPHASEMATCH";
    public static final String UNPHASED_SAMPLES_INFO_KEY = "UnphasedSamples";

    @Argument(fullName="phased-variants-file", doc="the phased VCF file", optional=false)
    private List<FeatureInput<VariantContext>> phased;

    @Argument(doc="File to which variants should be written", fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, optional = false)
    public File out = null;

    @Argument(doc="Additional genotype format tags to apply", fullName = "additional-tags", optional = true)
    public List<String> additionalTags = null;

    @Argument(doc="Missed alleles report", fullName = "missing-alleles-report", optional = true)
    public String missingAllelesReport = null;

    @Argument(doc="Concordance summary report", fullName = "concordance-summary-report", optional = true)
    public String concordanceSummaryReport = null;

    @Argument(doc="Discordant genotypes report", fullName = "discordant-genotypes-report", optional = true)
    public String discordantGenotypesReport = null;

    @Argument(doc="Unmatched variant report", fullName = "unmatched-variant-report", optional = true)
    public String unmatchedVariantReport = null;

    private VariantContextWriter vcfWriter;
    private PrintWriter missingAllelesReportWriter;
    private PrintWriter concordanceSummaryReportWriter;
    private PrintWriter discordantGenotypesReportWriter;
    private PrintWriter unmatchedVariantReportWriter;
    private ConcordanceSummary concordanceSummary;

    @Override
    public void onTraversalStart() {
        vcfWriter = createVCFWriter(out);

        final VCFHeader header = getHeaderForVariants();

        if ( ! header.hasGenotypingData() ) {
            throw new UserException.BadInput("VCF has no genotypes");
        }

        final Set<VCFHeaderLine> headerLines = new TreeSet<>(header.getMetaDataInSortedOrder());
        final Set<String> tagsToFind = new HashSet<>();
        tagsToFind.add(VCFConstants.PHASE_SET_KEY);
        tagsToFind.add(VCFConstants.PHASE_QUALITY_KEY);
        tagsToFind.addAll(additionalTags);

        final List<String> genotypeSamples = header.getGenotypeSamples();
        final HashSet<String> genotypeSampleSet = new HashSet<>(genotypeSamples);
        for (final FeatureInput<VariantContext> phasedVariantInput : phased) {
            final VCFHeader headerForFeatures = (VCFHeader) getHeaderForFeatures(phasedVariantInput);
            for (final String phasedSample : headerForFeatures.getGenotypeSamples()) {
                if (! genotypeSampleSet.contains(phasedSample)) {
                    logger.warn("input VCF does not contain phased sample " + phasedSample + ", ignoring that sample");
                }
            }
            for (Iterator<String> tagIterator = tagsToFind.iterator(); tagIterator.hasNext(); ) {
                String tag = tagIterator.next();
                if (headerForFeatures.hasFormatLine(tag)) {
                    headerLines.add(headerForFeatures.getFormatHeaderLine(tag));
                    tagIterator.remove();
                }
            }
            if (tagsToFind.isEmpty()) break;
        }

        if (! tagsToFind.isEmpty()) {
            throw new GATKException("Phased variant file does not contain a header line for requested genotype tags " + tagsToFind);
        }

        headerLines.add(new VCFInfoHeaderLine(DISCORDANT_GENTOYPES_INFO_KEY, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Samples in which the phased VCF had a discordant genotype"));
        headerLines.add(new VCFInfoHeaderLine(CONCORDANT_GENOTYPE_FILTER_STRINGS_INFO_KEY, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Filter strings for samples that had concordant genotypes with filters"));
        headerLines.add(new VCFInfoHeaderLine(NOPHASEMATCH_INFO_KEY, 0, VCFHeaderLineType.Flag, "Site had no match in phased input vcf"));
        headerLines.add(new VCFInfoHeaderLine(UNPHASED_SAMPLES_INFO_KEY, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Samples with matching heterozygous genotypes that were unphased in the phased VCF"));
        vcfWriter.writeHeader(new VCFHeader(headerLines, genotypeSamples));

        if (unmatchedVariantReport != null) {
            unmatchedVariantReportWriter = new PrintWriter(BucketUtils.createFile(unmatchedVariantReport));
            unmatchedVariantReportWriter.print("CHROM\tPOS\tTYPE\tFILTER\tREF\tALT\tAC\tQS\tMinAB\n");
        }

        if (missingAllelesReport != null) {
            missingAllelesReportWriter = new PrintWriter(BucketUtils.createFile(missingAllelesReport));
            missingAllelesReportWriter.print("CHROM\tPOS\tJOINT_ALLELES\tSAMPLE\tSAMPLE_ALLELES\n");
        }

        if (concordanceSummaryReport != null) {
            concordanceSummaryReportWriter = new PrintWriter(BucketUtils.createFile(concordanceSummaryReport));
        }

        if (discordantGenotypesReport != null) {
            discordantGenotypesReportWriter = new PrintWriter(BucketUtils.createFile(discordantGenotypesReport));
            discordantGenotypesReportWriter.print("CHROM\tPOS\tJOINT_TYPE\tJOINT_QUAL\tPHASED_QUAL\tHAPLOCALLED\tSAMPLE\tJOINT_GT\t\tJOINT_GQ\tJOINT_DP\tPHASED_GT\tPHASED_GQ\tPHASED_DP\n");
        }

        concordanceSummary = new ConcordanceSummary();
    }

    @Override
    public void apply(final VariantContext variant,
                      final ReadsContext readsContext,
                      final ReferenceContext referenceContext,
                      final FeatureContext featureContext) {
        final List<VariantContext> phasedOverlappingVariants = featureContext.getValues(phased);
        final List<VariantContext> phasedVariantsWithSameStart =
                phasedOverlappingVariants.stream().filter(v -> v.getStart() == variant.getStart()).collect(Collectors.toList());

        concordanceSummary.sawVariant(variant);

        VariantContextBuilder newVariantBuilder = new VariantContextBuilder(variant);
        final GenotypesContext genotypesContext = GenotypesContext.copy(variant.getGenotypes());

        final Map<String, String> sampleFilterStrings = new HashMap<>();
        final Set<String> discordantSamples = new HashSet<>();
        final Set<String> unphasedSamples = new HashSet<>();

        boolean foundMatch = false;
        for (final VariantContext phasedVariant : phasedVariantsWithSameStart) {
//            if (phasedVariant.isFiltered()) {
//                continue;
//            }

            final boolean phasedVariantRefLonger = phasedVariant.getReference().length() > variant.getReference().length();
            final GATKVariantContextUtils.AlleleMapper alleleMapping;

            if (phasedVariantRefLonger) {
                alleleMapping = GATKVariantContextUtils.resolveIncompatibleAlleles(phasedVariant.getReference(), variant, new LinkedHashSet<>());
            } else {

                alleleMapping = GATKVariantContextUtils.resolveIncompatibleAlleles(variant.getReference(), phasedVariant, new LinkedHashSet<>());
            }

//            if ((phasedVariantRefLonger &&  ! variant.getAlleles().stream().map(alleleMapping::remap).collect(Collectors.toList()).containsAll(phasedVariant.getAlleles())) ||
//                    ! variant.getAlleles().containsAll(phasedVariant.getAlleles().stream().map(alleleMapping::remap).collect(Collectors.toList()))) {
//                concordanceSummary.mismatchedAlleles(variant);
//                if (missingAllelesReportWriter != null) {
//                    missingAllelesReportWriter.print(variant.getContig() + "\t" + variant.getStart() + "\t" + variant.getAlleles() + "\t" + phasedVariant.getAlleles() + "\n");
//                }
//                continue;
//            }

            for (Genotype phasedVariantGenotype : phasedVariant.getGenotypes()) {

                final String sampleName = phasedVariantGenotype.getSampleName();
                if (! genotypesContext.containsSample(sampleName)) {
                    continue;
                }

                if (!phasedVariantGenotype.isCalled()) {
                    continue;
                }

                if (! foundMatch) {
                    concordanceSummary.variantSiteMatch(variant);
                    foundMatch = true;
                }

                if (!variantContainsPhasedAlleles(variant, alleleMapping, phasedVariantGenotype, phasedVariantRefLonger)) {
                    concordanceSummary.mismatchedAlleles(variant);
                    if (missingAllelesReportWriter != null) {
                        missingAllelesReportWriter.println(variant.getContig() + "\t" +
                                variant.getStart() + "\t" +
                                (phasedVariantRefLonger ? variant.getAlleles().stream().map(alleleMapping::remap).collect(Collectors.toList()) : variant.getAlleles()) + "\t" +
                                phasedVariantGenotype.getSampleName() + "\t" +
                                (phasedVariantRefLonger ? phasedVariant.getAlleles() : phasedVariantGenotype.getAlleles().stream().map(alleleMapping::remap).collect(Collectors.toList())));
                    }
                    discordantSamples.add(sampleName);
                    continue;
                }

                final Genotype variantGenotype = genotypesContext.get(sampleName);
                if (!isConcordant(variantGenotype, phasedVariantGenotype, alleleMapping, phasedVariantRefLonger)) {
                    concordanceSummary.discordantGenotype(variant);
                    if (discordantGenotypesReportWriter != null) {
                        discordantGenotypesReportWriter.println(variant.getContig() +
                                "\t" + variant.getStart() +
                                "\t" + variant.getType() +
                                "\t" + variant.getPhredScaledQual() +
                                "\t" + phasedVariant.getPhredScaledQual() +
                                "\t" + phasedVariant.getAttributeAsString("HAPLOCALLED", "") +
                                "\t" + sampleName +
                                "\t" + variantGenotype.getGenotypeString(false) +
                                "\t" + variantGenotype.getGQ() +
                                "\t" + variantGenotype.getDP() +
                                "\t" + phasedVariantGenotype.getGenotypeString(false) +
                                "\t" + phasedVariantGenotype.getGQ() +
                                "\t" + phasedVariantGenotype.getDP());
                    }
                    discordantSamples.add(sampleName);
                    continue;
                }

                concordanceSummary.concordantGenotype(variant);

                if (phasedVariant.isFiltered()) {
                    sampleFilterStrings.put(sampleName, String.join("-", phasedVariant.getFilters()));
                    continue;
                }

                if (phasedVariantGenotype.isHet() && ! phasedVariantGenotype.isPhased()) {
                    unphasedSamples.add(sampleName);
                    continue;
                }

                GenotypeBuilder genotypeBuilder = new GenotypeBuilder(variantGenotype)
                        .phased(phasedVariantGenotype.isPhased())
                        .alleles(phasedVariantGenotype.getAlleles().stream().map(alleleMapping::remap).collect(Collectors.toList()))
                        .attribute(VCFConstants.PHASE_SET_KEY, phasedVariantGenotype.getAnyAttribute(VCFConstants.PHASE_SET_KEY))
                        .attribute(VCFConstants.PHASE_QUALITY_KEY, phasedVariantGenotype.getAnyAttribute(VCFConstants.PHASE_QUALITY_KEY));

                for (String tag : additionalTags) {
                    genotypeBuilder = genotypeBuilder.attribute(tag, phasedVariantGenotype.getAnyAttribute(tag));
                }
                genotypesContext.replace(genotypeBuilder.make());
            }
        }


        if (!foundMatch) {
            if (unmatchedVariantReportWriter != null) {
                final Double minAB = getMinAB(variant);
                unmatchedVariantReportWriter.println(variant.getContig() +
                        "\t" + variant.getStart() +
                        "\t" + variant.getType() +
                        "\t" + variant.getFilters() +
                        "\t" + variant.getReference() +
                        "\t" + variant.getAlternateAlleles() +
                        "\t" + variant.getAttributeAsString(VCFConstants.ALLELE_COUNT_KEY, ".") +
                        "\t" + variant.getAttributeAsString("QD", ".") +
                        "\t" + ((minAB != null) ? minAB : "."));
            }

            newVariantBuilder = newVariantBuilder.attribute(NOPHASEMATCH_INFO_KEY,true);
        }

        newVariantBuilder = newVariantBuilder.genotypes(genotypesContext);
        if (! discordantSamples.isEmpty()) {
            newVariantBuilder = newVariantBuilder.attribute(DISCORDANT_GENTOYPES_INFO_KEY,
                    String.join(",", discordantSamples));
        }

        if (! sampleFilterStrings.isEmpty()) {
            newVariantBuilder = newVariantBuilder.attribute(CONCORDANT_GENOTYPE_FILTER_STRINGS_INFO_KEY,
                    sampleFilterStrings.entrySet().stream().map(e -> e.getKey() + ":" + e.getValue()).collect(Collectors.joining(",")));
        }

        if (! unphasedSamples.isEmpty()) {
            newVariantBuilder = newVariantBuilder.attribute(UNPHASED_SAMPLES_INFO_KEY,
                    String.join(",", unphasedSamples));
        }

        vcfWriter.add(newVariantBuilder.make());

    }


    private boolean variantContainsPhasedAlleles(final VariantContext variant, final GATKVariantContextUtils.AlleleMapper alleleMapping, final Genotype phasedVariantGenotype, final boolean phasedVariantRefLonger) {
        if (! phasedVariantRefLonger) {
            return variant.getAlleles().containsAll(phasedVariantGenotype.getAlleles().stream().map(alleleMapping::remap).collect(Collectors.toList()));
        } else {
            return variant.getAlleles().stream().map(alleleMapping::remap).collect(Collectors.toSet()).containsAll(phasedVariantGenotype.getAlleles());
        }
    }

    private boolean isConcordant(final Genotype variantGenotype,
                                 final Genotype phasedVariantGenotype,
                                 final GATKVariantContextUtils.AlleleMapper alleleMapping,
                                 final boolean phasedVariantRefLonger) {
        if (variantGenotype.getPloidy() != phasedVariantGenotype.getPloidy()) {
            return false;
        }

        if (phasedVariantRefLonger) {
            final List<Allele> remappedVariantAlleles = variantGenotype.getAlleles().stream().map(alleleMapping::remap).collect(Collectors.toList());
            for (final Allele phasedVariantAllele : phasedVariantGenotype.getAlleles()) {
                if (! remappedVariantAlleles.contains(phasedVariantAllele)) {
                    return false;
                }

                final Allele variantAllele = remappedVariantAlleles.get(remappedVariantAlleles.indexOf(phasedVariantAllele));
                if (variantGenotype.countAllele(variantAllele) != phasedVariantGenotype.countAllele(phasedVariantAllele)) {
                    return false;
                }

            }
        } else {
            for (final Allele phasedVariantAllele : phasedVariantGenotype.getAlleles()) {
                if (variantGenotype.countAllele(alleleMapping.remap(phasedVariantAllele)) != phasedVariantGenotype.countAllele(phasedVariantAllele)) {
                    return false;
                }
            }
        }
        return true;
    }

    @Override
    public void closeTool(){
        if(vcfWriter != null) {
            vcfWriter.close();
        }

        if (missingAllelesReportWriter != null) {
            missingAllelesReportWriter.close();
        }

        if (discordantGenotypesReportWriter != null) {
            discordantGenotypesReportWriter.close();
        }

        if (concordanceSummaryReportWriter != null) {
            concordanceSummary.write(concordanceSummaryReportWriter);
            concordanceSummaryReportWriter.close();
        }

        if (unmatchedVariantReportWriter != null) {
            unmatchedVariantReportWriter.close();
        }

    }

    static class ConcordanceSummary {
        private int variants;
        private int snpVariants;
        private int indelVariants;

        private int matchedVariants;
        private int matchedSnpVariants;
        private int matchedIndelVariants;

        private int genotypeAllelesMismatched;

        private int snpGenotypeAllelesMismatched;

        private int indelGenotypeAllelesMismatched;

        private int genotypesConcordant;
        private int genotypesDiscordant;

        private int snpGenotypesConcordant;
        private int snpGenotypesDiscordant;

        private int indelGenotypesConcordant;
        private int indelGenotypesDiscordant;

        public void sawVariant(final VariantContext variant) {
            variants++;
            if (variant.isSNP()) {
                snpVariants++;
            }
            if (variant.isIndel() || variant.isMixed()) {
                indelVariants++;
            }
        }

        public void variantSiteMatch(final VariantContext variant) {
            matchedVariants++;
            if (variant.isSNP()) {
                matchedSnpVariants++;
            }
            if (variant.isIndel() || variant.isMixed()) {
                matchedIndelVariants++;
            }

        }

        public void mismatchedAlleles(final VariantContext variant) {
            genotypeAllelesMismatched++;
            if (variant.isSNP()) {
                snpGenotypeAllelesMismatched++;
            }
            if (variant.isIndel()|| variant.isMixed()) {
                indelGenotypeAllelesMismatched++;
            }

        }

        public void discordantGenotype(final VariantContext variant) {
            genotypesDiscordant++;
            if (variant.isSNP()) {
                snpGenotypesDiscordant++;
            }
            if (variant.isIndel()|| variant.isMixed()) {
                indelGenotypesDiscordant++;
            }
        }

        public void concordantGenotype(final VariantContext variant) {
            genotypesConcordant++;
            if (variant.isSNP()) {
                snpGenotypesConcordant++;
            }
            if (variant.isIndel()|| variant.isMixed()) {
                indelGenotypesConcordant++;
            }
        }


        public void write(final PrintWriter concordanceSummaryReportWriter) {
            concordanceSummaryReportWriter.println("TYPE\tTOTAL\tSNP\tINDEL");
            concordanceSummaryReportWriter.println("VARIANTS\t" + variants + "\t" + snpVariants + "\t" + indelVariants);
            concordanceSummaryReportWriter.println("VARIANTS_MATCHED\t" + matchedVariants + "\t" + matchedSnpVariants + "\t" + matchedIndelVariants);
            concordanceSummaryReportWriter.println("GT_ALLELES_MISMATCHED\t" + genotypeAllelesMismatched + "\t" + snpGenotypeAllelesMismatched + "\t" + indelGenotypeAllelesMismatched);
            concordanceSummaryReportWriter.println("GT_CONCORDANT\t" + genotypesConcordant + "\t" + snpGenotypesConcordant + "\t" + indelGenotypesConcordant);
            concordanceSummaryReportWriter.println("GT_DISCORDANT\t" + genotypesDiscordant + "\t" + snpGenotypesDiscordant + "\t" + indelGenotypesDiscordant);
        }

    }

    // nicked from picard FilterVcf
    private static class Counts { int samples; int allele1; int allele2; }

    private Double getMinAB(final VariantContext variant) {
        if (variant.getHetCount() == 0) return null;
        final Map<List<Allele>, Counts> countsMap = new HashMap<List<Allele>, Counts>();

        for (final Genotype gt : variant.getGenotypesOrderedByName()) {
            if (gt.isNoCall() || !gt.isHet() || !gt.hasAD()) continue;

            final List<Allele> alleles = gt.getAlleles();
            Counts counts = countsMap.get(alleles);
            if (counts == null) {
                counts = new Counts();
                countsMap.put(alleles, counts);
            }

            counts.allele1 += gt.getAD()[variant.getAlleleIndex(alleles.get(0))];
            counts.allele2 += gt.getAD()[variant.getAlleleIndex(alleles.get(1))];
        }

        double minAB = Double.MAX_VALUE;
        for (final Counts counts : countsMap.values()) {
            final int total = counts.allele1 + counts.allele2;
            if (total > 0 && Math.min(counts.allele1, counts.allele2) / (double) total < minAB) {
                minAB = Math.min(counts.allele1, counts.allele2) / (double) total;
            }
        }
        return minAB;
    }

}
