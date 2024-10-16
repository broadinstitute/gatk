package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;

import java.util.Arrays;
import java.util.List;
import java.util.Set;
import java.util.Map;
import java.util.LinkedHashSet;
import java.util.HashSet;
import java.util.HashMap;
import java.util.stream.Collectors;

/**
 * Completes an initial series of cleaning steps for a VCF produced by the GATK-SV pipeline.
 *
 * <h3>Inputs</h3>
 * <ul>
 *     <li>
 *         VCF containing structural variant (SV) records from the GATK-SV pipeline.
 *     </li>
 *     <li>
 *         TODO
 *     </li>
 * </ul>
 *
 * <h3>Output</h3>
 * <ul>
 *     <li>
 *         Annotated VCF.
 *     </li>
 * </ul>
 *
 * <h3>Usage Example</h3>
 * <pre>
 *     gatk SVCleanPt1a \
 *       -V structural.vcf.gz \
 *       -O cleansed.vcf.gz
 *       --ped-file pedigree.ped
 *       --chrX chrX
 *       --chrY chrY
 *       --fail-list background_fail.txt
 *       --pass-list bothsides_pass.txt
 *       --sample-list sample_list.txt
 *       --revised-list revised_list.txt
 * </pre>
 *
 * <h3>Cleaning Steps</h3>
 * <ol>
 *     <li>
 *         Adds new FILTER and INFO tags to header.
 *     </li>
 *     <li>
 *         TODO
 *     </li>
 * </ol>
 */
@CommandLineProgramProperties(
        summary = "Clean and format structural variant VCFs per Step 1a",
        oneLineSummary = "Clean and format structural variant VCFs per Step 1a",
        programGroup = StructuralVariantDiscoveryProgramGroup.class
)
@BetaFeature
@DocumentedFeature
public final class SVCleanPt1a extends VariantWalker {
    public static final String PED_FILE_LONG_NAME = "ped-file";
    public static final String CHRX_LONG_NAME = "chrX";
    public static final String CHRY_LONG_NAME = "chrY";
    public static final String FAIL_LIST_LONG_NAME = "fail-list";
    public static final String PASS_LIST_LONG_NAME = "pass-list";
    public static final String OUTPUT_SAMPLES_LIST_LONG_NAME = "sample-list";
    public static final String OUTPUT_REVISED_EVENTS_LIST_LONG_NAME = "revised-list";

    @Argument(
            fullName = PED_FILE_LONG_NAME,
            doc = "Sample PED file"
    )
    private GATKPath pedFile;

    @Argument(
            fullName = CHRX_LONG_NAME,
            doc = "chrX column name"
    )
    private String chrX;

    @Argument(
            fullName = CHRY_LONG_NAME,
            doc = "chrY column name"
    )
    private String chrY;

    @Argument(
            fullName = FAIL_LIST_LONG_NAME,
            doc = "List of complex variants failing the background test"
    )
    private GATKPath failList;

    @Argument(
            fullName = PASS_LIST_LONG_NAME,
            doc = "List of complex variants passing both sides"
    )
    private GATKPath passList;

    @Argument(
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "Output VCF name"
    )
    private GATKPath outputVcf;

    @Argument(
            fullName = OUTPUT_SAMPLES_LIST_LONG_NAME,
            doc="Output list of samples"
    )
    private GATKPath outputSamplesList;

    @Argument(
            fullName = OUTPUT_REVISED_EVENTS_LIST_LONG_NAME,
            doc="Output list of revised genotyped events"
    )
    private GATKPath outputRevisedEventsList;

    private VariantContextWriter vcfWriter;
    private BufferedWriter samplesWriter;
    private BufferedWriter revisedEventsWriter;

    private Map<String, Integer> sampleSexMap;
    private Set<String> failSet;
    private Set<String> passSet;
    private final Set<String> writtenRevisedEvents = new HashSet<>();

    private static final int MIN_ALLOSOME_EVENT_SIZE = 5000;


    @Override
    public void onTraversalStart() {
        // Read supporting files
        sampleSexMap = readPedFile(pedFile);
        failSet = readLastColumn(failList);
        passSet = readLastColumn(passList);

        // Create header without the 'UNRESOLVED' INFO line
        final VCFHeader header = getHeaderForVariants();
        Set<VCFHeaderLine> newHeaderLines = new LinkedHashSet<>();
        for (final VCFHeaderLine line : header.getMetaDataInInputOrder()) {
            if (!(line instanceof VCFInfoHeaderLine) || !((VCFInfoHeaderLine) line).getID().equals(GATKSVVCFConstants.UNRESOLVED)) {
                newHeaderLines.add(line);
            }
        }

        // Add new header lines
        VCFHeader newHeader = new VCFHeader(newHeaderLines, header.getGenotypeSamples());
        newHeader.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.HIGH_SR_BACKGROUND, 0, VCFHeaderLineType.Flag, "High number of SR splits in background samples indicating messy region"));
        newHeader.addMetaDataLine(new VCFFilterHeaderLine(GATKSVVCFConstants.UNRESOLVED, "Variant is unresolved"));
        newHeader.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.BOTHSIDES_SUPPORT, 0, VCFHeaderLineType.Flag, "Variant has read-level support for both sides of breakpoint"));

        // Write header
        vcfWriter = createVCFWriter(outputVcf);
        vcfWriter.writeHeader(newHeader);

        // Create output writers
        try {
            revisedEventsWriter = new BufferedWriter(new FileWriter(outputRevisedEventsList.toPath().toFile()));
            samplesWriter = new BufferedWriter(new FileWriter(outputSamplesList.toPath().toFile()));

            for (final String sample : header.getGenotypeSamples()) {
                samplesWriter.write(sample);
                samplesWriter.newLine();
            }
            samplesWriter.flush();
        } catch (IOException e) {
            throw new RuntimeException("Error creating output file", e);
        }
    }

    @Override
    public void closeTool() {
        try {
            if (vcfWriter != null) {
                vcfWriter.close();
            }

            if (samplesWriter != null) {
                samplesWriter.close();
            }
            if (revisedEventsWriter != null) {
                revisedEventsWriter.close();
            }
        } catch (IOException e) {
            throw new RuntimeException("Error closing output file", e);
        }
    }

    @Override
    public void apply(final VariantContext variant, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext) {
        VariantContextBuilder variantBuilder = new VariantContextBuilder(variant);
        final List<Genotype> processedGenotypes = processGenotypes(variant);
        variantBuilder.genotypes(processedGenotypes);
        processVariant(variant, variantBuilder);
        vcfWriter.add(variantBuilder.make());
    }

    private List<Genotype> processGenotypes(VariantContext variant) {
        return variant.getGenotypes().stream()
                .map(genotype -> {
                    GenotypeBuilder genotypeBuilder = new GenotypeBuilder(genotype);
                    processEVGenotype(genotype, genotypeBuilder);
                    // processSVTypeGenotype(variant, genotype, genotypeBuilder);
                    processAllosomesGenotype(variant, genotype, genotypeBuilder);
                    return genotypeBuilder.make();
                })
                .collect(Collectors.toList());
    }

    private void processVariant(VariantContext variant, VariantContextBuilder builder) {
        // processSVType(variant, builder);
        processVarGQ(variant, builder);
        processMultiallelic(builder);
        processUnresolved(variant, builder);
        processNoisyEvents(variant, builder);
        processBothsidesSupportEvents(variant, builder);
    }

    private void processEVGenotype(Genotype genotype, GenotypeBuilder genotypeBuilder) {
        if (genotype.hasExtendedAttribute(GATKSVVCFConstants.EV)) {
            String evAttribute = (String) genotype.getExtendedAttribute(GATKSVVCFConstants.EV);
            final int evIndex = Integer.parseInt(evAttribute);
            if (evIndex >= 0 && evIndex < GATKSVVCFConstants.evValues.size()) {
                genotypeBuilder.attribute(GATKSVVCFConstants.EV, GATKSVVCFConstants.evValues.get(evIndex));
            }
        }
    }

    private void processSVTypeGenotype(VariantContext variant, Genotype genotype, GenotypeBuilder genotypeBuilder) {
        final String svType = variant.getAttributeAsString(GATKSVVCFConstants.SVTYPE, null);
        if (svType != null && variant.getAlleles().stream().noneMatch(allele -> allele.getDisplayString().contains(GATKSVVCFConstants.ME))) {
            List<Allele> newGenotypeAlleles = Arrays.asList(
                    variant.getReference(),
                    Allele.create("<" + svType + ">", false)
            );
            genotypeBuilder.alleles(newGenotypeAlleles);
        }
    }

    private void processAllosomesGenotype(VariantContext variant, Genotype genotype, GenotypeBuilder genotypeBuilder) {
        final String chromosome = variant.getContig();
        if (chromosome.equals(chrX) || chromosome.equals(chrY)) {
            final boolean isY = chromosome.equals(chrY);
            final String svType = variant.getAttributeAsString(GATKSVVCFConstants.SVTYPE, "");
            if ((svType.equals(GATKSVVCFConstants.SYMB_ALT_STRING_DEL) || svType.equals(GATKSVVCFConstants.SYMB_ALT_STRING_DUP)) &&
                    (variant.getEnd() - variant.getStart() >= MIN_ALLOSOME_EVENT_SIZE)) {
                final String sampleName = genotype.getSampleName();
                final int sex = sampleSexMap.get(sampleName);
                if (sex == 1 && isRevisableEvent(variant, isY)) { // Male
                    writeRevisedEvents(variant);
                    adjustMaleGenotype(genotype, genotypeBuilder, svType);
                } else if (sex == 2 && isY) { // Female
                    genotypeBuilder.alleles(Arrays.asList(Allele.NO_CALL, Allele.NO_CALL));
                } else if (sex == 0) { // Unknown
                    genotypeBuilder.alleles(Arrays.asList(Allele.NO_CALL, Allele.NO_CALL));
                }
            }
        }
    }

    private void adjustMaleGenotype(Genotype genotype, GenotypeBuilder genotypeBuilder, String svType) {
        if (genotype.hasExtendedAttribute(GATKSVVCFConstants.RD_CN)) {
            final int rdCN = Integer.parseInt(genotype.getExtendedAttribute(GATKSVVCFConstants.RD_CN).toString());
            genotypeBuilder.attribute(GATKSVVCFConstants.RD_CN, rdCN + 1);

            final Allele refAllele = genotype.getAllele(0);
            final Allele altAllele = genotype.getAllele(1);
            if (svType.equals(GATKSVVCFConstants.SYMB_ALT_STRING_DEL)) {
                if (rdCN >= 1) genotypeBuilder.alleles(Arrays.asList(refAllele, refAllele));
                else if (rdCN == 0) genotypeBuilder.alleles(Arrays.asList(refAllele, altAllele));
            } else if (svType.equals(GATKSVVCFConstants.SYMB_ALT_STRING_DUP)) {
                if (rdCN <= 1) genotypeBuilder.alleles(Arrays.asList(refAllele, refAllele));
                else if (rdCN == 2) genotypeBuilder.alleles(Arrays.asList(refAllele, altAllele));
                else genotypeBuilder.alleles(Arrays.asList(altAllele, altAllele));
            }
        }
    }

    private boolean isRevisableEvent(VariantContext variant, boolean isY) {
        final List<Genotype> genotypes = variant.getGenotypes();
        final int[] maleCounts = new int[4];
        final int[] femaleCounts = new int[4];
        for (final Genotype genotype : genotypes) {
            final String sampleName = genotype.getSampleName();
            final Integer sex = sampleSexMap.get(sampleName);
            if (sex == null) continue;

            final int rdCN = (int) genotype.getExtendedAttribute(GATKSVVCFConstants.RD_CN, -1);
            if (rdCN == -1) continue;

            final int rdCNVal = Math.min(rdCN, 3);
            if (sex == 1) {
                maleCounts[rdCNVal]++;
            } else if (sex == 2) {
                femaleCounts[rdCNVal]++;
            }
        }

        final double maleMedian = calcMedian(maleCounts);
        final double femaleMedian = calcMedian(femaleCounts);
        return maleMedian == 1.0 && (isY ? femaleMedian == 0.0 : femaleMedian == 2.0);
    }

    private double calcMedian(int[] counts) {
        final int total = Arrays.stream(counts).sum();
        if (total == 0) return Double.NaN;

        final double target = total / 2.0;
        int runningTotal = 0;
        for (int i = 0; i < 4; i++) {
            runningTotal += counts[i];
            if (runningTotal == target) {
                return i + 0.5;
            } else if (runningTotal > target) {
                return i;
            }
        }
        throw new RuntimeException("Error calculating median");
    }

    private void processSVType(VariantContext variant, VariantContextBuilder builder) {
        final String svType = variant.getAttributeAsString(GATKSVVCFConstants.SVTYPE, null);
        if (svType != null && variant.getAlleles().stream().noneMatch(allele -> allele.getDisplayString().contains(GATKSVVCFConstants.ME))) {
            final Allele refAllele = variant.getReference();
            final Allele altAllele = Allele.create("<" + svType + ">", false);
            List<Allele> newAlleles = Arrays.asList(refAllele, altAllele);
            builder.alleles(newAlleles);
        }
    }

    private void processVarGQ(VariantContext variant, VariantContextBuilder builder) {
        if (variant.hasAttribute(GATKSVVCFConstants.VAR_GQ)) {
            final double varGQ = variant.getAttributeAsDouble(GATKSVVCFConstants.VAR_GQ, 0);
            builder.rmAttribute(GATKSVVCFConstants.VAR_GQ);
            builder.log10PError(varGQ / -10.0);
        }
    }

    private void processMultiallelic(VariantContextBuilder builder) {
        builder.rmAttribute(GATKSVVCFConstants.MULTIALLELIC);
    }

    private void processUnresolved(VariantContext variant, VariantContextBuilder builder) {
        if (variant.hasAttribute(GATKSVVCFConstants.UNRESOLVED)) {
            builder.rmAttribute(GATKSVVCFConstants.UNRESOLVED);
            builder.filter(GATKSVVCFConstants.UNRESOLVED);
        }
    }

    private void processNoisyEvents(VariantContext variant, VariantContextBuilder builder) {
        if (failSet.contains(variant.getID())) {
            builder.attribute(GATKSVVCFConstants.HIGH_SR_BACKGROUND, true);
        }
    }

    private void processBothsidesSupportEvents(VariantContext variant, VariantContextBuilder builder) {
        if (passSet.contains(variant.getID())) {
            builder.attribute(GATKSVVCFConstants.BOTHSIDES_SUPPORT, true);
        }
    }

    private Set<String> readLastColumn(GATKPath filePath) {
        try {
            return Files.lines(Paths.get(filePath.toString()))
                    .filter(line -> !line.trim().isEmpty() && !line.startsWith("#"))
                    .map(line -> {
                        int lastTabIndex = line.lastIndexOf('\t');
                        return lastTabIndex != -1 ? line.substring(lastTabIndex + 1).trim() : line.trim();
                    })
                    .collect(Collectors.toSet());
        } catch (IOException e) {
            throw new RuntimeException("Error reading variant list file: " + filePath, e);
        }
    }

    private Map<String, Integer> readPedFile(GATKPath pedFile) {
        Map<String, Integer> sampleSexMap = new HashMap<>();
        try (BufferedReader reader = new BufferedReader(new FileReader(pedFile.toPath().toFile()))) {
            String line;
            while ((line = reader.readLine()) != null) {
                if (line.startsWith("#")) continue;
                final String[] fields = line.split("\t");
                if (fields.length >= 5) {
                    final String sampleName = fields[1];
                    final int sex = Integer.parseInt(fields[4]);
                    sampleSexMap.put(sampleName, sex);
                }
            }
        } catch (IOException e) {
            throw new RuntimeException("Error reading PED file", e);
        }
        return sampleSexMap;
    }

    private void writeRevisedEvents(VariantContext variant) {
        final String variantId = variant.getID();
        if (!writtenRevisedEvents.contains(variantId)) {
            try {
                revisedEventsWriter.write(variantId);
                revisedEventsWriter.newLine();
                writtenRevisedEvents.add(variantId);
            } catch (IOException e) {
                throw new RuntimeException("Error writing to revised events output file", e);
            }
        }
    }
}