package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFHeaderLines;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVCallRecordUtils;
import org.broadinstitute.hellbender.tools.sv.cluster.SVCollapser;
import org.broadinstitute.hellbender.tools.sv.cluster.SVDeduplicator;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Creates multi-sample structural variant (SV) VCF from a collection of SV VCFs. Supported types include biallelic DEL,
 * DUP, INS, INV, and BND. Input files may contain mutually exclusive samples and/or variants.
 *
 * Each record must have the following INFO fields:
 * <ul>
 *     <li>END Integer, end coordinate on CHROM, may not precede POS</li>
 *     <li>SVLEN Integer, the length of the event or -1 if undefined, e.g. for BND.</li>
 *     <li>SVTYPE String, type of structural variant</li>
 *     <li>STRANDS String, from the set {'++', '+-', '-+', '--'} (INV/BND only)</li>
 *     <li>ALGORITHMS String list, SV calling algorithms or "depth" if from a CNV caller</li>
 *     <li>CHR2 String and END2 Integer, mate coordinate (BND only)</li>
 * </ul>
 *
 * Additional INFO and FORMAT fields are cleared from each variant. Samples with non-reference genotypes are flagged
 * with a raw call (RC) field in the output VCF. Duplicate records of the same type, coordinates, and SVLEN are collapsed.
 *
 * <h3>Inputs</h3>
 *
 * <ul>
 *     <li>
 *         VCF(s)
 *     </li>
 * </ul>
 *
 * <h3>Output</h3>
 *
 * <ul>
 *     <li>
 *         Merged VCF
 *     </li>
 * </ul>
 *
 * <h3>Usage example</h3>
 *
 * <pre>
 *     gatk MergeSVRecords \
 *      -V sample_1.caller_A.vcf.gz \
 *      -V sample_1.caller_B.vcf.gz \
 *      -V sample_2.caller_A.vcf.gz
 *      -O merged.vcf.gz
 * </pre>
 *
 * @author Mark Walker &lt;markw@broadinstitute.org&gt;
 */

@CommandLineProgramProperties(
        summary = "Merges structural variant VCFs",
        oneLineSummary = "Merges structural variant VCFs",
        programGroup = StructuralVariantDiscoveryProgramGroup.class
)
@ExperimentalFeature
@DocumentedFeature
public final class SVPreprocessRecords extends MultiVariantWalker {
    @Argument(
            doc = "Output combined VCF",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private GATKPath outputFile;

    @Argument(
            doc = "Low memory mode. Drops hom-ref and no-call genotype fields and emits all as diploid no-calls.",
            fullName = "low-memory"
    )
    private boolean lowMemoryMode;

    final static List<Allele> DEFAULT_ALLELES = Arrays.asList(Allele.NO_CALL);

    private VariantContextWriter writer;
    private List<SVCallRecord> records;
    private Set<String> samples;
    private SAMSequenceDictionary dictionary;
    private SVDeduplicator<SVCallRecord> deduplicator;

    private String currentContig;
    private int currentPosition = 0;

    @Override
    public void onTraversalStart() {
        super.onTraversalStart();
        dictionary = getBestAvailableSequenceDictionary();
        if (dictionary == null) {
            throw new UserException("Reference sequence dictionary required");
        }
        records = new ArrayList<>();
        samples = getSamplesForVariants();
        writer = createVCFWriter(outputFile);
        writer.writeHeader(createVcfHeader());
        final SVCollapser collapser = new SVCollapser(null);
        deduplicator = new SVDeduplicator<>(collapser::collapse, dictionary);
    }

    @Override
    public Object onTraversalSuccess() {
        super.onTraversalSuccess();
        flushRecords();
        return null;
    }

    @Override
    public void apply(final VariantContext variant,
                      final ReadsContext readsContext,
                      final ReferenceContext referenceContext,
                      final FeatureContext featureContext) {
        Utils.validate(variant.getNAlleles() == 2, "Records must be biallelic");
        final SVCallRecord record = preprocessInputGenotypes(SVCallRecordUtils.create(variant));
        if (record.getPositionA() != currentPosition || !record.getContigA().equals(currentContig)) {
            flushRecords();
            currentPosition = record.getPositionA();
            currentContig = record.getContigA();
        }
        records.add(record);
        progressMeter.update(record.getPositionAInterval());
    }

    private void flushRecords() {
        deduplicator.deduplicateSortedItems(records).stream().map(this::createVariant).forEachOrdered(writer::add);
        records.clear();
    }

    @Override
    public void closeTool() {
        super.closeTool();
        if (writer != null) {
            writer.close();
        }
    }

    private SVCallRecord preprocessInputGenotypes(final SVCallRecord record) {
        // Saves memory to only retain non-ref genotypes
        if (lowMemoryMode) {
            final GenotypesContext genotypes = GenotypesContext.copy(record.getGenotypes().stream()
                    .filter(SVCallRecordUtils::isAltGenotype)
                    .map(GenotypeBuilder::new)
                    .map(GenotypeBuilder::noAttributes)
                    .map(GenotypeBuilder::make)
                    .collect(Collectors.toList()));
            return SVCallRecordUtils.copyCallWithNewGenotypes(record, genotypes);
        } else {
            return record;
        }
    }

    private VariantContext createVariant(final SVCallRecord call) {
        final VariantContextBuilder builder = SVCallRecordUtils.getVariantBuilder(call);
        builder.genotypes(SVCallRecordUtils.fillMissingSamplesWithGenotypes(builder.getGenotypes(), samples, DEFAULT_ALLELES, Collections.emptyMap()));
        return builder.make();
    }

    private VCFHeader createVcfHeader() {
        final VCFHeader header = new VCFHeader(getDefaultToolVCFHeaderLines(), samples);
        header.setVCFHeaderVersion(VCFHeaderVersion.VCF4_2);
        header.setSequenceDictionary(dictionary);

        // Copy from inputs
        getHeaderForVariants().getFormatHeaderLines().forEach(header::addMetaDataLine);
        getHeaderForVariants().getInfoHeaderLines().forEach(header::addMetaDataLine);

        // Info lines
        header.addMetaDataLine(VCFStandardHeaderLines.getInfoLine(VCFConstants.END_KEY));
        header.addMetaDataLine(GATKSVVCFHeaderLines.getInfoLine(GATKSVVCFConstants.SVLEN));
        header.addMetaDataLine(GATKSVVCFHeaderLines.getInfoLine(GATKSVVCFConstants.SVTYPE));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.END2_ATTRIBUTE, 1, VCFHeaderLineType.String, "Second position"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.CONTIG2_ATTRIBUTE, 1, VCFHeaderLineType.String, "Second contig"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.STRANDS_ATTRIBUTE, 1, VCFHeaderLineType.String, "First and second strands"));
        header.addMetaDataLine(new VCFInfoHeaderLine("##INFO=<ID=ALGORITHMS,Number=.,Type=String,Description=\"Source algorithms\">", header.getVCFHeaderVersion()));

        // Format lines
        header.addMetaDataLine(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_KEY));

        return header;
    }

}
