package org.broadinstitute.hellbender.tools.picard.analysis.directed;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.SequenceUtil;
import org.broadinstitute.hellbender.metrics.MetricAccumulationLevel;
import org.broadinstitute.hellbender.metrics.PerUnitMetricCollector;
import org.broadinstitute.hellbender.metrics.SAMRecordAndReference;
import org.broadinstitute.hellbender.metrics.SAMRecordAndReferenceMultiLevelCollector;

import java.util.*;

public final class RrbsMetricsCollector extends SAMRecordAndReferenceMultiLevelCollector<RrbsMetrics, Long> {
	private final int minReadLength;
	private final double maxMismatchRate;
	private final int cQualityThreshold;
	private final int nextBaseQualityThreshold;

	public RrbsMetricsCollector(final Set<MetricAccumulationLevel> accumulationLevels, final List<SAMReadGroupRecord> samRgRecords,
								final int cQualityThreshold, final int nextBaseQualityThreshold, final int minReadLength,
								final double maxMismatchRate) {
		this.cQualityThreshold = cQualityThreshold;
		this.nextBaseQualityThreshold = nextBaseQualityThreshold;
		this.minReadLength = minReadLength;
		this.maxMismatchRate = maxMismatchRate;
		setup(accumulationLevels, samRgRecords);
	}

	@Override
	protected PerUnitMetricCollector<RrbsMetrics, Long, SAMRecordAndReference> makeChildCollector(final String sample, final String library, final String readGroup) {
		return new PerUnitRrbsMetricsCollector(sample, library, readGroup);
	}

	private class PerUnitRrbsMetricsCollector implements PerUnitMetricCollector<RrbsMetrics, Long, SAMRecordAndReference> {
		final String sample;
		final String library;
		final String readGroup;

		// Counters for CpG & non-CpG seen/converted sites
		int nCytoConverted = 0;
		int nCytoTotal = 0;
		final Histogram<CpgLocation> cpgTotal = new Histogram<>();
		final Histogram<CpgLocation> cpgConverted = new Histogram<>();

		// Counters for QC filters used in the final metrics
		int mappedRecordCount = 0;
		int smallReadCount = 0;
		int mismatchCount = 0;
		int noCpgCount = 0;

		// Final metrics calculated once all the reads are done
		double cytoConversionRate;
		double cpgConversionRate;
		int nCpgSeen;
		int nCpgConverted;
		double coverageMean;
		int coverageMedian;

		public PerUnitRrbsMetricsCollector(final String sample, final String library, final String readGroup) {
			this.sample = sample;
			this.library = library;
			this.readGroup = readGroup;
		}

		public void acceptRecord(final SAMRecordAndReference args) {
			mappedRecordCount++;

			final SAMRecord samRecord = args.getSamRecord();
			final ReferenceSequence referenceSequence = args.getReferenceSequence();

			final byte[] readBases = samRecord.getReadBases();
			final byte[] readQualities = samRecord.getBaseQualities();
			final byte[] refBases = referenceSequence.getBases();

			if (samRecord.getReadLength() < minReadLength) {
				smallReadCount++;
				return;
			} else if (SequenceUtil.countMismatches(samRecord, refBases, true) > Math.round(samRecord.getReadLength() * maxMismatchRate)) {
				mismatchCount++;
				return;
			}

			// We only record non-CpG C sites if there was at least one CpG in the read, keep track of
			// the values for this record and then apply to the global totals if valid
			int recordCpgs = 0;

			for (final AlignmentBlock alignmentBlock : samRecord.getAlignmentBlocks()) {
				final int blockLength = alignmentBlock.getLength();
				final int refFragmentStart = alignmentBlock.getReferenceStart() - 1;
				final int readFragmentStart = alignmentBlock.getReadStart() - 1;

				final byte[] refFragment = getFragment(refBases, refFragmentStart, blockLength);
				final byte[] readFragment = getFragment(readBases, readFragmentStart, blockLength);
				final byte[] readQualityFragment = getFragment(readQualities, readFragmentStart, blockLength);

				if (samRecord.getReadNegativeStrandFlag()) {
					// In the case of a negative strand, reverse (and complement for base arrays) the reference,
					// reads & qualities so that it can be treated as a positive strand for the rest of the process
					SequenceUtil.reverseComplement(refFragment);
					SequenceUtil.reverseComplement(readFragment);
					SequenceUtil.reverseQualities(readQualityFragment);
				}

				for (int i=0; i < blockLength-1; i++) {
					final int curRefIndex = getCurRefIndex(refFragmentStart, blockLength, i, samRecord.getReadNegativeStrandFlag());

					// Look at a 2-base window to see if we're on a CpG site, and if so check for conversion
					// (CG -> TG). We do not consider ourselves to be on a CpG site if we're on the last base of a read
					if ((SequenceUtil.basesEqual(refFragment[i], SequenceUtil.C)) &&
							(SequenceUtil.basesEqual(refFragment[i + 1], SequenceUtil.G))) {
						// We want to catch the case where there's a CpG in the reference, even if it is not valid
						// to prevent the C showing up as a non-CpG C down below. Otherwise this could have been all
						// in one if statement
						if (isValidCpg(refFragment, readFragment, readQualityFragment, i)) {
							recordCpgs++;
							final CpgLocation curLocation = new CpgLocation(samRecord.getReferenceName(), curRefIndex);
							cpgTotal.increment(curLocation);
							if (SequenceUtil.isBisulfiteConverted(readFragment[i], refFragment[i])) {
								cpgConverted.increment(curLocation);
							}
						}
						i++;
					} else if (isC(refFragment[i], readFragment[i]) && isAboveCytoQcThreshold(readQualities, i) &&
							SequenceUtil.bisulfiteBasesEqual(false, readFragment[i + 1], refFragment[i + 1])) {
						// C base in the reference that's not associated with a CpG
						nCytoTotal++;
						if (SequenceUtil.isBisulfiteConverted(readFragment[i], refFragment[i])) {
							nCytoConverted++;
						}
					}
				}
			}

			if (recordCpgs == 0) {
				noCpgCount++;
			}
		}

		public void finish() {
			cytoConversionRate = nCytoTotal == 0 ? 0 : nCytoConverted / (double)nCytoTotal;
			nCpgSeen = (int)cpgTotal.getSumOfValues();
			nCpgConverted = (int)cpgConverted.getSumOfValues();
			cpgConversionRate = nCpgSeen == 0 ? 0 : nCpgConverted / (double)nCpgSeen;
			coverageMean = cpgTotal.getMeanBinSize();
			coverageMedian = (int)cpgTotal.getMedianBinSize();
		}

		@Override
		public void addMetricsToFile(final MetricsFile<RrbsMetrics, Long> metricsFile) {
			// Create both the summary and detail metrics & add them to the RrbsMetrics container class for
			// the downstream code to use as desired
			final RrbsSummaryMetrics summaryMetrics = buildSummaryMetrics();
			final List<RrbsCpgDetailMetrics> detailMetrics = buildDetailMetrics();
			final RrbsMetrics rrbsMetrics = new RrbsMetrics(summaryMetrics, detailMetrics);
			metricsFile.addMetric(rrbsMetrics);
		}

		private RrbsSummaryMetrics buildSummaryMetrics() {
			final RrbsSummaryMetrics summaryMetrics = new RrbsSummaryMetrics();
			summaryMetrics.SAMPLE = sample;
			summaryMetrics.READ_GROUP = readGroup;
			summaryMetrics.LIBRARY = library;
			summaryMetrics.READS_ALIGNED = mappedRecordCount;
			summaryMetrics.NON_CPG_BASES = nCytoTotal;
			summaryMetrics.NON_CPG_CONVERTED_BASES = nCytoConverted;
			summaryMetrics.PCT_NON_CPG_BASES_CONVERTED = cytoConversionRate;
			summaryMetrics.CPG_BASES_SEEN = nCpgSeen;
			summaryMetrics.CPG_BASES_CONVERTED = nCpgConverted;
			summaryMetrics.PCT_CPG_BASES_CONVERTED = cpgConversionRate;
			summaryMetrics.MEAN_CPG_COVERAGE = coverageMean;
			summaryMetrics.MEDIAN_CPG_COVERAGE = coverageMedian;
			summaryMetrics.READS_IGNORED_SHORT = smallReadCount;
			summaryMetrics.READS_WITH_NO_CPG = noCpgCount;
			summaryMetrics.READS_IGNORED_MISMATCHES = mismatchCount;
			return summaryMetrics;
		}

		private List<RrbsCpgDetailMetrics> buildDetailMetrics() {
			final List<RrbsCpgDetailMetrics> detailMetrics = new ArrayList<>();
			for (final CpgLocation key : cpgTotal.keySet()) {
				final RrbsCpgDetailMetrics cpgMetric = new RrbsCpgDetailMetrics();
				cpgMetric.SAMPLE = sample;
				cpgMetric.READ_GROUP = readGroup;
				cpgMetric.LIBRARY = library;
				cpgMetric.SEQUENCE_NAME = key.getSequence();
				cpgMetric.POSITION = key.getPosition();
				cpgMetric.TOTAL_SITES = (int)cpgTotal.get(key).getValue();
				cpgMetric.CONVERTED_SITES = cpgConverted.containsKey(key) ? (int)cpgConverted.get(key).getValue() : 0;
				cpgMetric.PCT_CONVERTED = cpgMetric.CONVERTED_SITES == 0 ? 0 : cpgMetric.CONVERTED_SITES / (double)cpgMetric.TOTAL_SITES;
				detailMetrics.add(cpgMetric);
			}
			return detailMetrics;
		}
	}

	private byte[] getFragment(final byte[] fullArray, final int fragmentStart, final int length) {
		return Arrays.copyOfRange(fullArray, fragmentStart, fragmentStart + length);
	}

	/**
	 * True if there's a C in the reference as well as read (possibly bisulfite converted)
	 */
	private boolean isC(final byte refBase, final byte readBase) {
		return (SequenceUtil.basesEqual(refBase, SequenceUtil.C) && SequenceUtil.bisulfiteBasesEqual(readBase, refBase));
	}

	/**
	 * Checks a pair of bases for CpG status & quality thresholds
	 */
	private boolean isValidCpg(final byte[] refBases, final byte[] readBases, final byte[] readQualities, final int index) {
		return isC(refBases[index], readBases[index]) && SequenceUtil.basesEqual(refBases[index + 1], readBases[index + 1]) &&
				isAboveCytoQcThreshold(readQualities, index);
	}

	/**
	 * Any cyto base (CpG context or not) needs to be above a specified quality threshold. Similarly, the neighboring
	 * base on the 3' side must also be above another threshold unless the cyto base is the final base in the 3'
	 * direction.
	 */
	private boolean isAboveCytoQcThreshold(final byte[] readQualities, final int index) {
		return ((index < readQualities.length - 1) && (readQualities[index] >= cQualityThreshold) &&
				(readQualities[index+1] >= nextBaseQualityThreshold));
	}

	/**
	 * Accounts for the fact that negative strand counts have been reversed
	 */
	private int getCurRefIndex(final int refStart, final int blockLength, final int idx, final boolean isNegative) {
		return isNegative ? refStart + (blockLength - 1) - idx - 1 : refStart + idx;
	}
}

/**
 * Used to keep track of the location of CpG sites
 */
final class CpgLocation implements Comparable<CpgLocation> {
	private final String sequence;
	private final Integer position;

	public CpgLocation(final String sequence, final int position) {
		this.sequence = sequence;
		this.position = position;
	}

	@Override
	public int compareTo(final CpgLocation other) {
		final int seqComp = sequence.compareTo(other.sequence);
		return seqComp == 0 ? position.compareTo(other.position) : seqComp;
	}

	@Override
	public boolean equals(final Object other) {
		if (this == other) {
			return true;
		}

		if (other == null || getClass() != other.getClass()) {
			return false;
		}

		final CpgLocation that = (CpgLocation)other;

		return (sequence.equals(that.sequence)) && (position.equals(that.position));
	}

	@Override
	public int hashCode() {
		int result = sequence.hashCode();
		result = 31 * result + position;
		return result;
	}

	public String getSequence() {
		return sequence;
	}

	public Integer getPosition() {
		return position;
	}
}
