package org.broadinstitute.hellbender.tools.sv.cluster;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.LazyGenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFEncoder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderVersion;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Expands a sparse output record to an all-sample VCF record by text splicing instead of per-sample
 * {@link Genotype} construction and encoding.
 *
 * <p>Motivation: at biobank scale (hundreds of thousands of samples) the cost of writing an SV site is
 * dominated by building one default (non-carrier) genotype object per sample and running the VCF encoder over
 * all of them, roughly 90 ms per site at 240k samples, even though nearly every column is one of a handful of
 * identical default strings determined only by the sample's ploidy on the contig.</p>
 *
 * <p>Method. The caller builds a small {@link VariantContext} holding the site's carrier genotypes plus one
 * default genotype per ploidy class on the contig (see {@link #representativeSamples}), produced by the same
 * fill code the full record would use. That small record is encoded with htsjdk's own {@link VCFEncoder}
 * against a header restricted to its samples. Because the FORMAT key list is a union over genotypes and every
 * default genotype of a class is identical, the small record yields exactly the FORMAT column, exactly the
 * carrier columns, and exactly the default column of each class that the full record would produce. Default
 * columns are laid out once per (contig, FORMAT, default texts) into a cached template covering all samples in
 * header order; per site, the carrier columns are spliced into that template. The result is attached to the
 * record as a {@link LazyGenotypesContext} carrying the pre-rendered {@code FORMAT\tcol\t...\tcol} string,
 * which {@link VCFEncoder} appends verbatim when writing, so no genotype objects are ever created for
 * non-carriers. Should anything force a decode, the attached parser rebuilds the genotypes with
 * {@link VCFCodec}, so the record stays correct (only slower).</p>
 *
 * <p>Output is byte-identical to the full encode; see {@code TemplatedGenotypeColumnEncoderTest}.</p>
 */
final class TemplatedGenotypeColumnEncoder {
    private static final Logger logger = LogManager.getLogger(TemplatedGenotypeColumnEncoder.class);
    private static final char FIELD_SEPARATOR = '\t';
    /** CHROM POS ID REF ALT QUAL FILTER INFO FORMAT */
    private static final int NUM_FIXED_FIELDS = 9;

    private final VCFHeader header;
    private final List<String> headerSamples;
    private final Map<String, Integer> sampleIndex;
    private final PloidyTable ploidyTable;
    private final boolean allowMissingFieldsInHeader;
    private final VCFCodec fallbackCodec;
    private boolean warnedFallbackDecode = false;

    // ----- per-contig state -----
    private String contig = null;
    /** Ploidy class id of each header sample on {@link #contig}. */
    private int[] sampleClass;
    /** Header sample indices per ploidy class, in header order. */
    private int[][] classMembers;
    /**
     * Class id of samples absent from the ploidy table, or -1 if all samples have a ploidy. Such samples can never
     * receive a default genotype (the fill path rejects them), so they are never representatives and must appear
     * as carriers in every record; {@link #expand} enforces this to fail exactly where the full fill would.
     */
    private int unknownPloidyClass;
    /** Templates for the current contig, keyed by FORMAT column and per-class default column texts. */
    private final Map<String, Template> templates = new HashMap<>();

    /** All sample columns for one (contig, FORMAT, default texts), tab-separated in header order, with offsets. */
    private static final class Template {
        final String text;
        /** Column start offset per header sample index. */
        final int[] start;
        /** Column end offset (exclusive) per header sample index. */
        final int[] end;

        Template(final String text, final int[] start, final int[] end) {
            this.text = text;
            this.start = start;
            this.end = end;
        }
    }

    /**
     * @param header output header; its sample list defines column order
     * @param ploidyTable ploidy for every header sample on every contig that will be written
     * @param allowMissingFieldsInHeader must match the output writer's setting (GATK: {@code --lenient})
     */
    TemplatedGenotypeColumnEncoder(final VCFHeader header, final PloidyTable ploidyTable,
                                   final boolean allowMissingFieldsInHeader) {
        this.header = Utils.nonNull(header);
        Utils.validateArg(header.hasGenotypingData(), "Header has no samples");
        this.headerSamples = header.getGenotypeSamples();
        this.sampleIndex = new HashMap<>(SVUtils.hashMapCapacity(headerSamples.size()));
        for (int i = 0; i < headerSamples.size(); i++) {
            sampleIndex.put(headerSamples.get(i), i);
        }
        this.ploidyTable = Utils.nonNull(ploidyTable);
        this.allowMissingFieldsInHeader = allowMissingFieldsInHeader;
        this.fallbackCodec = new VCFCodec();
        this.fallbackCodec.setVCFHeader(header, VCFHeaderVersion.VCF4_2);
    }

    /**
     * Chooses, for each ploidy class present on {@code contig}, one sample that is not a carrier (the first such
     * sample in header order). Filling default genotypes for exactly these samples into a sparse record makes its
     * FORMAT keys and per-class default columns identical to those of the fully expanded record. A class whose
     * every sample is a carrier gets no representative, and needs none.
     */
    Set<String> representativeSamples(final String contig, final Set<String> carrierSamples) {
        ensureContig(contig);
        final Set<String> representatives = new LinkedHashSet<>();
        for (int c = 0; c < classMembers.length; c++) {
            if (c == unknownPloidyClass) {
                continue;
            }
            for (final int i : classMembers[c]) {
                final String sample = headerSamples.get(i);
                if (!carrierSamples.contains(sample)) {
                    representatives.add(sample);
                    break;
                }
            }
        }
        return representatives;
    }

    /**
     * Expands a sparse record to all header samples.
     *
     * @param sparse record whose genotypes are exactly its carriers plus default genotypes for
     *               {@code representatives}, built with the same fill code as the full record
     * @param representatives the result of {@link #representativeSamples} for this record's carriers
     * @return a record with identical site fields whose genotypes are a pre-rendered all-sample column string
     */
    VariantContext expand(final VariantContext sparse, final Set<String> representatives) {
        ensureContig(sparse.getContig());
        final int numSparse = sparse.getNSamples();

        // Order the sparse record's samples by header position.
        final int[] sparseIndices = new int[numSparse];
        int k = 0;
        for (final Genotype g : sparse.getGenotypes()) {
            final Integer index = sampleIndex.get(g.getSampleName());
            Utils.validate(index != null, () -> "Sample not in header: " + g.getSampleName());
            sparseIndices[k++] = index;
        }
        Arrays.sort(sparseIndices);
        requireUnknownPloidySamplesPresent(sparse.getContig(), sparseIndices);
        final List<String> sparseSamples = new ArrayList<>(numSparse);
        for (final int i : sparseIndices) {
            sparseSamples.add(headerSamples.get(i));
        }

        // Encode carriers + representatives with htsjdk against a header restricted to those samples. This is the
        // sole source of FORMAT, carrier, and default column text, so formatting cannot drift from the encoder.
        final VCFHeader sparseHeader = new VCFHeader(header.getMetaDataInInputOrder(), sparseSamples);
        final String line = new VCFEncoder(sparseHeader, allowMissingFieldsInHeader, false).encode(sparse);
        final String[] fields = splitFields(line, NUM_FIXED_FIELDS + numSparse);
        final String format = fields[NUM_FIXED_FIELDS - 1];

        // Default column text per ploidy class, read off the representatives' columns.
        final String[] classDefaults = new String[classMembers.length];
        for (k = 0; k < numSparse; k++) {
            if (representatives.contains(sparseSamples.get(k))) {
                classDefaults[sampleClass[sparseIndices[k]]] = fields[NUM_FIXED_FIELDS + k];
            }
        }
        final Template template = templateFor(format, classDefaults);

        // Splice every sparse column (carriers, and representatives whose text equals the template's) into the
        // template. Non-carrier columns are copied straight from the template text.
        final StringBuilder columns = new StringBuilder(format.length() + 1 + template.text.length() + 64);
        columns.append(format).append(FIELD_SEPARATOR);
        int cursor = 0;
        for (k = 0; k < numSparse; k++) {
            final int i = sparseIndices[k];
            columns.append(template.text, cursor, template.start[i]);
            columns.append(fields[NUM_FIXED_FIELDS + k]);
            cursor = template.end[i];
        }
        columns.append(template.text, cursor, template.text.length());

        final LazyGenotypesContext lazy = new LazyGenotypesContext(
                new FallbackParser(sparse.getAlleles(), sparse.getContig(), sparse.getStart()),
                columns.toString(), headerSamples.size());
        return new VariantContextBuilder(sparse).genotypesNoValidation(lazy).make();
    }

    /**
     * The full fill path looks up ploidy only for samples that lack a genotype, so a sample missing from the ploidy
     * table is legal as long as it is a carrier in every record. Mirror that: such samples are never given a
     * template default, and a record in which one of them is not a carrier fails with the same error the fill
     * path would raise.
     */
    private void requireUnknownPloidySamplesPresent(final String recordContig, final int[] sortedSparseIndices) {
        if (unknownPloidyClass < 0) {
            return;
        }
        for (final int member : classMembers[unknownPloidyClass]) {
            if (Arrays.binarySearch(sortedSparseIndices, member) < 0) {
                ploidyTable.get(headerSamples.get(member), recordContig); // throws: sample not in ploidy records
                throw new IllegalStateException("Sample " + headerSamples.get(member) + " has no ploidy on " + recordContig);
            }
        }
    }

    private void ensureContig(final String newContig) {
        if (newContig.equals(contig)) {
            return;
        }
        final int numSamples = headerSamples.size();
        final Integer[] ploidies = new Integer[numSamples];
        final Set<Integer> distinct = new LinkedHashSet<>();
        boolean anyUnknown = false;
        for (int i = 0; i < numSamples; i++) {
            final String sample = headerSamples.get(i);
            if (ploidyTable.contains(sample)) {
                ploidies[i] = ploidyTable.get(sample, newContig);
                distinct.add(ploidies[i]);
            } else {
                anyUnknown = true;
            }
        }
        final int[] classPloidies = distinct.stream().mapToInt(Integer::intValue).sorted().toArray();
        final int numClasses = classPloidies.length + (anyUnknown ? 1 : 0);
        unknownPloidyClass = anyUnknown ? classPloidies.length : -1;
        final Map<Integer, Integer> classOfPloidy = new HashMap<>();
        for (int c = 0; c < classPloidies.length; c++) {
            classOfPloidy.put(classPloidies[c], c);
        }
        final int[] classSizes = new int[numClasses];
        sampleClass = new int[numSamples];
        for (int i = 0; i < numSamples; i++) {
            sampleClass[i] = ploidies[i] == null ? unknownPloidyClass : classOfPloidy.get(ploidies[i]);
            classSizes[sampleClass[i]]++;
        }
        classMembers = new int[numClasses][];
        final int[] fill = new int[numClasses];
        for (int c = 0; c < numClasses; c++) {
            classMembers[c] = new int[classSizes[c]];
        }
        for (int i = 0; i < numSamples; i++) {
            final int c = sampleClass[i];
            classMembers[c][fill[c]++] = i;
        }
        templates.clear();
        contig = newContig;
    }

    private Template templateFor(final String format, final String[] classDefaults) {
        final StringBuilder keyBuilder = new StringBuilder(format);
        for (final String d : classDefaults) {
            keyBuilder.append('\n').append(d == null ? "" : d);
        }
        final String key = keyBuilder.toString();
        Template template = templates.get(key);
        if (template == null) {
            template = buildTemplate(classDefaults);
            templates.put(key, template);
        }
        return template;
    }

    private Template buildTemplate(final String[] classDefaults) {
        final int numSamples = headerSamples.size();
        int estimate = numSamples;
        for (int c = 0; c < classDefaults.length; c++) {
            estimate += classMembers[c].length * (classDefaults[c] == null ? 0 : classDefaults[c].length());
        }
        final StringBuilder text = new StringBuilder(estimate);
        final int[] start = new int[numSamples];
        final int[] end = new int[numSamples];
        for (int i = 0; i < numSamples; i++) {
            if (i > 0) {
                text.append(FIELD_SEPARATOR);
            }
            start[i] = text.length();
            final String d = classDefaults[sampleClass[i]];
            if (d != null) {
                text.append(d);
            }
            end[i] = text.length();
        }
        return new Template(text.toString(), start, end);
    }

    private static String[] splitFields(final String line, final int expectedFields) {
        final String[] fields = new String[expectedFields];
        int from = 0;
        for (int f = 0; f < expectedFields; f++) {
            int to = line.indexOf(FIELD_SEPARATOR, from);
            if (to < 0) {
                Utils.validate(f == expectedFields - 1, () -> "Encoded record has fewer than " + expectedFields + " fields");
                to = line.length();
                if (to > 0 && line.charAt(to - 1) == '\n') {
                    to--;
                }
            } else {
                Utils.validate(f < expectedFields - 1, () -> "Encoded record has more than " + expectedFields + " fields");
            }
            fields[f] = line.substring(from, to);
            from = to + 1;
        }
        return fields;
    }

    /**
     * Decodes the pre-rendered column string back into genotypes if a consumer asks for them. Uses the same codec
     * logic the VCF reader would apply, so the decoded genotypes match what re-reading the output would give.
     */
    private final class FallbackParser implements LazyGenotypesContext.LazyParser {
        private final List<Allele> alleles;
        private final String contig;
        private final int start;

        FallbackParser(final List<Allele> alleles, final String contig, final int start) {
            this.alleles = alleles;
            this.contig = contig;
            this.start = start;
        }

        @Override
        public LazyGenotypesContext.LazyData parse(final Object data) {
            if (!warnedFallbackDecode) {
                logger.warn("Pre-rendered genotype columns were decoded back into genotype objects; output is " +
                        "unaffected but this path is slow at large sample counts");
                warnedFallbackDecode = true;
            }
            return fallbackCodec.createGenotypeMap((String) data, alleles, contig, start);
        }
    }
}
