package org.broadinstitute.hellbender.tools.walkers.contamination;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SimpleSampleLocatableMetadata;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;

import java.io.File;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * Created by David Benjamin on 2/14/17.
 */
public class PileupSummary implements Locatable {
    private final String contig;
    private final int position;
    private final int refCount;
    private final int altCount;
    private final int otherAltsCount;
    private final int totalCount;
    private final byte refBase;
    private final byte altBase;

    private final double alleleFrequency;

    public PileupSummary(String contig, int position, int refCount, int altCount, int otherAltsCount, byte refBase, byte altBase, double alleleFrequency) {
        this.contig = contig;
        this.position = position;
        this.altCount = altCount;
        this.refCount = refCount;
        this.otherAltsCount = otherAltsCount;
        this.totalCount = refCount + altCount + otherAltsCount;
        this.refBase = refBase;
        this.altBase = altBase;
        this.alleleFrequency = alleleFrequency;
    }

    public PileupSummary(final VariantContext vc, final ReadPileup pileup) {
        contig = vc.getContig();
        position = vc.getStart();
        alleleFrequency = vc.getAttributeAsDouble(VCFConstants.ALLELE_FREQUENCY_KEY, 0);
        altBase = vc.getAlternateAllele(0).getBases()[0];
        refBase = vc.getReference().getBases()[0];
        final int[] baseCounts = pileup.getBaseCounts();
        altCount = baseCounts[BaseUtils.simpleBaseToBaseIndex(altBase)];
        refCount = baseCounts[BaseUtils.simpleBaseToBaseIndex(refBase)];
        totalCount = (int) MathUtils.sum(baseCounts);
        otherAltsCount = totalCount - altCount - refCount;
    }

    @Override
    public String getContig() { return contig; }

    @Override
    public int getStart() { return position; }

    @Override
    public int getEnd() { return position; }

    public int getAltCount() {
        return altCount;
    }
    public int getRefCount() {
        return refCount;
    }
    public int getOtherAltCount() { return otherAltsCount; }
    public int getTotalCount() {
        return totalCount;
    }
    public byte getRefBase() { return refBase; }
    public byte getAltBase() { return altBase; }
    public double getAlleleFrequency() {
        return alleleFrequency;
    }
    public double getRefFrequency() {
        return 1 - alleleFrequency;
    }
    public double getAltFraction() {
        return totalCount == 0 ? 0 : (double) altCount / totalCount;
    }
    public double getMinorAlleleFraction() {
        final double altFraction = getAltFraction();
        return FastMath.min(altFraction, 1 - altFraction);
    }

    // useful for tests when we want to write without worrying about the header stuff
    public static void writeWithAutomaticSequenceDictionary(final String sample, final List<PileupSummary> records, final File outputFile) {
        final SAMSequenceDictionary seqDict = new SAMSequenceDictionary();
        final Map<String, List<PileupSummary>> recordsByContig = records.stream().collect(Collectors.groupingBy(PileupSummary::getContig));
        for (final Map.Entry<String, List<PileupSummary>> entry : recordsByContig.entrySet()) {
            final String contig = entry.getKey();
            final int maxPosition = entry.getValue().stream().mapToInt(PileupSummary::getEnd).max().orElse(1);
            seqDict.addSequence(new SAMSequenceRecord(contig, maxPosition));
        }
        final SampleLocatableMetadata metadata = new SimpleSampleLocatableMetadata(sample, seqDict);
        new PileupSummaryCollection(metadata, records).write(outputFile);
    }

    // Takes a list of PileupSummaryTable files and write them all in one output file in order
    public static void writeToFile(final List<File> inputFiles, final File output) {
        Utils.nonEmpty(inputFiles);
        SampleLocatableMetadata metadata = null;
        final List<PileupSummary> records = new ArrayList<>();
        for (final File inputFile : inputFiles) {
            final PileupSummaryCollection collection = new PileupSummaryCollection(inputFile);
            records.addAll(collection.getRecords());
            final SampleLocatableMetadata fileMetadata = collection.getMetadata();
            metadata = metadata == null ? fileMetadata : metadata;

            Utils.validate(metadata.getSampleName() == fileMetadata.getSampleName(), "different sample names");
            metadata.getSequenceDictionary().assertSameDictionary(fileMetadata.getSequenceDictionary());
        }

        new PileupSummaryCollection(metadata, records).write(output);
    }

    public static ImmutablePair<String, List<PileupSummary>> readFromFile(final File tableFile) {
        final PileupSummaryCollection collection = new PileupSummaryCollection(tableFile);
        final List<PileupSummary> pileupSummaries = collection.getRecords();
        final String sample = collection.getMetadata().getSampleName();
        return ImmutablePair.of(sample, pileupSummaries);
    }

    public static class PileupSummaryComparator implements Comparator<PileupSummary> {
        final SAMSequenceDictionary sequenceDictionary;
        final List<String> contigsInOrder;

        public PileupSummaryComparator(final SAMSequenceDictionary sequenceDictionary){
            this.sequenceDictionary = sequenceDictionary;
            contigsInOrder = sequenceDictionary.getSequences().stream()
                    .map(ssr -> ssr.getSequenceName())
                    .collect(Collectors.toList());
        }

        @Override
        public int compare(PileupSummary ps1, PileupSummary ps2) {
            // Use Contig Index in case the contig name is e.g. chr2
            final int contigIndex1 = contigsInOrder.indexOf(ps1.getContig());
            final int contigIndex2 = contigsInOrder.indexOf(ps2.getContig());

            if (contigIndex1 != contigIndex2){
                return Integer.compare(contigIndex1, contigIndex2);
            } else {
                return Integer.compare(ps1.getStart(), ps2.getStart());
            }
        }
    }
}
