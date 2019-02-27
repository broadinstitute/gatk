package org.broadinstitute.hellbender.utils.read.mergealignment;

import htsjdk.samtools.*;
import htsjdk.samtools.filter.FilteringSamIterator;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.reference.ReferenceSequenceFileWalker;
import htsjdk.samtools.util.*;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.LogManager;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.runtime.ProgressLogger;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

/**
 * Abstract class that coordinates the general task of taking in a set of alignment information,
 * possibly in SAM format, possibly in other formats, and merging that with the set of all reads
 * for which alignment was attempted, stored in an unmapped SAM file.
 * <p/>
 * The order of processing is as follows:
 * <p/>
 * 1.  Get records from the unmapped bam and the alignment data
 * 2.  Merge the alignment information and public tags ONLY from the aligned SAMRecords
 * 3.  Do additional modifications -- handle clipping, trimming, etc.
 * 4.  Fix up mate information on paired reads
 * 5.  Do a final calculation of the NM and UQ tags.
 * 6.  Write the records to the output file.
 * <p/>
 * Concrete subclasses which extend AbstractAlignmentMerger should implement getQueryNameSortedAlignedRecords.
 * If these records are not in queryname order, mergeAlignment will throw an IllegalStateException.
 * <p/>
 * Subclasses may optionally implement ignoreAlignment(), which can be used to skip over certain alignments.
 *
 * @author ktibbett@broadinstitute.org
 */
public abstract class AbstractAlignmentMerger {

    public static final int MAX_RECORDS_IN_RAM = 500000;

    private static final char[] RESERVED_ATTRIBUTE_STARTS = {'X', 'Y', 'Z'};

    private final Logger logger = LogManager.getLogger(AbstractAlignmentMerger.class);
    private final ProgressLogger progress = new ProgressLogger(this.logger, 1000000, "Written to sorting collection in queryname order", "records");

    private final File unmappedBamFile;
    private final File targetBamFile;
    private final SAMSequenceDictionary sequenceDictionary;
    private ReferenceSequenceFileWalker refSeq = null;
    private final boolean clipAdapters;
    private final boolean bisulfiteSequence;
    private SAMProgramRecord programRecord;
    private final boolean alignedReadsOnly;
    private final SAMFileHeader header;
    private final List<String> attributesToRetain = new ArrayList<>();
    private final List<String> attributesToRemove = new ArrayList<>();
    protected final File referenceFasta;
    private final Integer read1BasesTrimmed;
    private final Integer read2BasesTrimmed;
    private final List<SamPairUtil.PairOrientation> expectedOrientations;
    private final SAMFileHeader.SortOrder sortOrder;
    private MultiHitAlignedReadIterator alignedIterator = null;
    private boolean clipOverlappingReads = true;
    private final PrimaryAlignmentSelectionStrategy primaryAlignmentSelectionStrategy;
    private boolean keepAlignerProperPairFlags = false;
    private boolean addMateCigar = false;

    private final SamRecordFilter alignmentFilter = new SamRecordFilter() {
        @Override
        public boolean filterOut(final SAMRecord record) {
            return ignoreAlignment(record);
        }

        @Override
        public boolean filterOut(final SAMRecord first, final SAMRecord second) {
            throw new UnsupportedOperationException("Paired SamRecordFilter not implemented!");
        }
    };
    private boolean includeSecondaryAlignments = true;

    protected abstract CloseableIterator<SAMRecord> getQuerynameSortedAlignedRecords();

    protected boolean ignoreAlignment(final SAMRecord sam) { return false; } // default implementation

    /**
     * Constructor
     *
     * @param unmappedBamFile                   The BAM file that was used as the input to the aligner, which will
     *                                          include info on all the reads that did not map.  Required.
     * @param targetBamFile                     The file to which to write the merged SAM records. Required.
     * @param referenceFasta                    The reference sequence for the map files. Required.
     * @param clipAdapters                      Whether adapters marked in unmapped BAM file should be marked as
     *                                          soft clipped in the merged bam. Required.
     * @param bisulfiteSequence                 Whether the reads are bisulfite sequence (used when calculating the
     *                                          NM and UQ tags). Required.
     * @param alignedReadsOnly                  Whether to output only those reads that have alignment data
     * @param programRecord                     Program record for target file SAMRecords created.
     * @param attributesToRetain                private attributes from the alignment record that should be
     *                                          included when merging.  This overrides the exclusion of
     *                                          attributes whose tags start with the reserved characters
     *                                          of X, Y, and Z
     * @param attributesToRemove                attributes from the alignment record that should be
     *                                          removed when merging.  This overrides attributesToRetain if they share
     *                                          common tags.
     * @param read1BasesTrimmed                 The number of bases trimmed from start of read 1 prior to alignment.  Optional.
     * @param read2BasesTrimmed                 The number of bases trimmed from start of read 2 prior to alignment.  Optional.
     * @param expectedOrientations              A List of SamPairUtil.PairOrientations that are expected for
     *                                          aligned pairs.  Used to determine the properPair flag.
     * @param sortOrder                         The order in which the merged records should be output.  If null,
     *                                          output will be coordinate-sorted
     * @param primaryAlignmentSelectionStrategy What to do when there are multiple primary alignments, or multiple
     *                                          alignments but none primary, for a read or read pair.
     * @param addMateCigar                      True if we are to add or maintain the mate CIGAR (MC) tag, false if we are to remove or not include.
     */
    public AbstractAlignmentMerger(final File unmappedBamFile, final File targetBamFile,
                                   final File referenceFasta, final boolean clipAdapters,
                                   final boolean bisulfiteSequence, final boolean alignedReadsOnly,
                                   final SAMProgramRecord programRecord, final List<String> attributesToRetain,
                                   final List<String> attributesToRemove,
                                   final Integer read1BasesTrimmed, final Integer read2BasesTrimmed,
                                   final List<SamPairUtil.PairOrientation> expectedOrientations,
                                   final SAMFileHeader.SortOrder sortOrder,
                                   final PrimaryAlignmentSelectionStrategy primaryAlignmentSelectionStrategy,
                                   final boolean addMateCigar) {
        IOUtil.assertFileIsReadable(unmappedBamFile);
        IOUtil.assertFileIsWritable(targetBamFile);
        IOUtil.assertFileIsReadable(referenceFasta);

        this.unmappedBamFile = unmappedBamFile;
        this.targetBamFile = targetBamFile;
        this.referenceFasta = referenceFasta;

        this.refSeq = new ReferenceSequenceFileWalker(referenceFasta);
        this.sequenceDictionary = refSeq.getSequenceDictionary();
        if (this.sequenceDictionary == null) {
            throw new UserException("No sequence dictionary found for " + referenceFasta.getAbsolutePath() +
                    ".  Use CreateSequenceDictionary.jar to create a sequence dictionary.");
        }

        this.clipAdapters = clipAdapters;
        this.bisulfiteSequence = bisulfiteSequence;
        this.alignedReadsOnly = alignedReadsOnly;

        this.header = new SAMFileHeader();
        this.sortOrder = sortOrder != null ? sortOrder : SAMFileHeader.SortOrder.coordinate;
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        if (programRecord != null) {
            setProgramRecord(programRecord);
        }
        header.setSequenceDictionary(this.sequenceDictionary);
        if (attributesToRetain != null) {
            this.attributesToRetain.addAll(attributesToRetain);
        }
        if (attributesToRemove != null) {
            this.attributesToRemove.addAll(attributesToRemove);
            // attributesToRemove overrides attributesToRetain
            if (!this.attributesToRetain.isEmpty()) {
                for (String attribute : this.attributesToRemove) {
                    if (this.attributesToRetain.contains(attribute)) {
                        logger.info("Overriding retaining the " + attribute + " tag since remove overrides retain.");
                        this.attributesToRetain.remove(attribute);
                    }
                }
            }
        }
        this.read1BasesTrimmed = read1BasesTrimmed;
        this.read2BasesTrimmed = read2BasesTrimmed;
        this.expectedOrientations = expectedOrientations;

        this.primaryAlignmentSelectionStrategy = primaryAlignmentSelectionStrategy;

        this.addMateCigar = addMateCigar;
    }

    /**
     * Do this unconditionally, not just for aligned records, for two reasons:
     * - An unaligned read has been processed by the aligner, so it is more truthful.
     * - When chaining additional PG records, having all the records in the output file refer to the same PG
     * record means that only one chain will need to be created, rather than a chain for the mapped reads
     * and a separate chain for the unmapped reads.
     */
    private void maybeSetPgTag(final SAMRecord rec) {
        if (this.programRecord != null) {
            rec.setAttribute(ReservedTagConstants.PROGRAM_GROUP_ID, this.programRecord.getProgramGroupId());
        }
    }

    /**
     * /**
     * Merges the alignment data with the non-aligned records from the source BAM file.
     */
    public void mergeAlignment(final File referenceFasta) {
        // Open the file of unmapped records and write the read groups to the the header for the merged file
        final SamReader unmappedSam = SamReaderFactory.makeDefault().referenceSequence(referenceFasta).open(this.unmappedBamFile);

        final CloseableIterator<SAMRecord> unmappedIterator = unmappedSam.iterator();
        this.header.setReadGroups(unmappedSam.getFileHeader().getReadGroups());

        int aligned = 0;
        int unmapped = 0;

        // Get the aligned records and set up the first one
        alignedIterator = new MultiHitAlignedReadIterator(new FilteringSamIterator(getQuerynameSortedAlignedRecords(), alignmentFilter), primaryAlignmentSelectionStrategy);
        HitsForInsert nextAligned = nextAligned();

        // Check that the program record we are going to insert is not already used in the unmapped SAM
        // Must come after calling getQuerynameSortedAlignedRecords() in case opening the aligned records
        // sets the program group
        if (getProgramRecord() != null) {
            for (final SAMProgramRecord pg : unmappedSam.getFileHeader().getProgramRecords()) {
                if (pg.getId().equals(getProgramRecord().getId())) {
                    throw new GATKException("Program Record ID already in use in unmapped BAM file.");
                }
            }
        }

        // Create the sorting collection that will write the records in the coordinate order
        // to the final bam file
        final SortingCollection<SAMRecord> sorted = SortingCollection.newInstance(
                SAMRecord.class, new BAMRecordCodec(header), new SAMRecordCoordinateComparator(),
                MAX_RECORDS_IN_RAM);

        while (unmappedIterator.hasNext()) {
            // Load next unaligned read or read pair.
            final SAMRecord rec = unmappedIterator.next();

            rec.setHeaderStrict(this.header);
            maybeSetPgTag(rec);

            final SAMRecord secondOfPair;
            if (rec.getReadPairedFlag()) {
                secondOfPair = unmappedIterator.next();
                secondOfPair.setHeaderStrict(this.header);
                maybeSetPgTag(secondOfPair);

                // Validate that paired reads arrive as first of pair followed by second of pair
                if (!rec.getReadName().equals(secondOfPair.getReadName()))
                    throw new GATKException("Second read from pair not found in unmapped bam: " + rec.getReadName() + ", " + secondOfPair.getReadName());

                if (!rec.getFirstOfPairFlag())
                    throw new GATKException("First record in unmapped bam is not first of pair: " + rec.getReadName());
                if (!secondOfPair.getReadPairedFlag())
                    throw new GATKException("Second record in unmapped bam is not marked as paired: " + secondOfPair.getReadName());
                if (!secondOfPair.getSecondOfPairFlag())
                    throw new GATKException("Second record in unmapped bam is not second of pair: " + secondOfPair.getReadName());
            } else {
                secondOfPair = null;
            }

            // See if there are alignments for current unaligned read or read pair.
            if (nextAligned != null && rec.getReadName().equals(nextAligned.getReadName())) {
                // If there are multiple alignments for a read (pair), then the unaligned SAMRecord must be cloned
                // before copying info from the aligned record to the unaligned.
                final boolean clone = nextAligned.numHits() > 1 || nextAligned.hasSupplementalHits();
                SAMRecord r1Primary = null, r2Primary = null;

                if (rec.getReadPairedFlag()) {
                    for (int i = 0; i < nextAligned.numHits(); ++i) {
                        // firstAligned or secondAligned may be null, if there wasn't an alignment for the end,
                        // or if the alignment was rejected by ignoreAlignment.
                        final SAMRecord firstAligned = nextAligned.getFirstOfPair(i);
                        final SAMRecord secondAligned = nextAligned.getSecondOfPair(i);

                        final boolean isPrimaryAlignment = (firstAligned != null && !firstAligned.isSecondaryOrSupplementary()) ||
                                (secondAligned != null && !secondAligned.isSecondaryOrSupplementary());

                        final SAMRecord firstToWrite;
                        final SAMRecord secondToWrite;
                        if (clone) {
                            firstToWrite = ReadUtils.cloneSAMRecord(rec);
                            secondToWrite = ReadUtils.cloneSAMRecord(secondOfPair);
                        } else {
                            firstToWrite = rec;
                            secondToWrite = secondOfPair;
                        }

                        // If these are the primary alignments then stash them for use on any supplemental alignments
                        if (isPrimaryAlignment) {
                            r1Primary = firstToWrite;
                            r2Primary = secondToWrite;
                        }

                        transferAlignmentInfoToPairedRead(firstToWrite, secondToWrite, firstAligned, secondAligned);

                        // Only write unmapped read when it has the mate info from the primary alignment.
                        if (!firstToWrite.getReadUnmappedFlag() || isPrimaryAlignment) {
                            addIfNotFiltered(sorted, firstToWrite);
                            if (firstToWrite.getReadUnmappedFlag()) ++unmapped;
                            else ++aligned;
                        }
                        if (!secondToWrite.getReadUnmappedFlag() || isPrimaryAlignment) {
                            addIfNotFiltered(sorted, secondToWrite);
                            if (!secondToWrite.getReadUnmappedFlag()) ++aligned;
                            else ++unmapped;
                        }
                    }

                    // Take all of the supplemental reads which had been stashed and add them (as appropriate) to sorted
                    for (final boolean isRead1 : new boolean[]{true, false}) {
                        final List<SAMRecord> supplementals = isRead1 ? nextAligned.getSupplementalFirstOfPairOrFragment() : nextAligned.getSupplementalSecondOfPair();
                        final SAMRecord sourceRec = isRead1 ? rec : secondOfPair;
                        final SAMRecord matePrimary = isRead1 ? r2Primary : r1Primary;

                        for (final SAMRecord supp : supplementals) {
                            final SAMRecord out = ReadUtils.cloneSAMRecord(sourceRec);
                            transferAlignmentInfoToFragment(out, supp);
                            if (matePrimary != null) SamPairUtil.setMateInformationOnSupplementalAlignment(out, matePrimary, addMateCigar);
                            ++aligned;
                            addIfNotFiltered(sorted, out);
                        }
                    }
                } else {
                    for (int i = 0; i < nextAligned.numHits(); ++i) {
                        final SAMRecord recToWrite = clone ? ReadUtils.cloneSAMRecord(rec) : rec;
                        transferAlignmentInfoToFragment(recToWrite, nextAligned.getFragment(i));
                        addIfNotFiltered(sorted, recToWrite);
                        if (recToWrite.getReadUnmappedFlag()) ++unmapped;
                        else ++aligned;
                    }
                    // Take all of the supplemental reads which had been stashed and add them (as appropriate) to sorted
                    for (final SAMRecord supplementalRec : nextAligned.getSupplementalFirstOfPairOrFragment()) {
                        final SAMRecord recToWrite = ReadUtils.cloneSAMRecord(rec);
                        transferAlignmentInfoToFragment(recToWrite, supplementalRec);
                        addIfNotFiltered(sorted, recToWrite);
                        ++aligned;
                    }
                }
                nextAligned = nextAligned();
            } else {
                // There was no alignment for this read or read pair.
                if (nextAligned != null &&
                        SAMRecordQueryNameComparator.compareReadNames(rec.getReadName(), nextAligned.getReadName()) > 0) {
                    throw new IllegalStateException("Aligned record iterator (" + nextAligned.getReadName() +
                            ") is behind the unmapped reads (" + rec.getReadName() + ")");
                }
                // No matching read from alignedIterator -- just output reads as is.
                if (!alignedReadsOnly) {
                    sorted.add(rec);
                    ++unmapped;
                    if (secondOfPair != null) {
                        sorted.add(secondOfPair);
                        ++unmapped;
                    }
                }
            }
        }
        unmappedIterator.close();
        Utils.validate(!alignedIterator.hasNext(), () -> "Reads remaining on alignment iterator: " + alignedIterator.next().getReadName() + "!");
        alignedIterator.close();

        // Write the records to the output file in specified sorted order,
        header.setSortOrder(this.sortOrder);
        final boolean presorted = this.sortOrder == SAMFileHeader.SortOrder.coordinate;
        final SAMFileWriter writer = new SAMFileWriterFactory().makeWriter(header, presorted, this.targetBamFile, referenceFasta);
        writer.setProgressLogger(
                new ProgressLogger(logger, (int) 1e7, "Wrote", "records from a sorting collection"));
        final ProgressLogger finalProgress = new ProgressLogger(logger, 10000000, "Written in coordinate order to output", "records");

        for (final SAMRecord rec : sorted) {
            if (!rec.getReadUnmappedFlag()) {
                if (refSeq != null) {
                    final byte[] referenceBases = refSeq.get(sequenceDictionary.getSequenceIndex(rec.getReferenceName())).getBases();
                    rec.setAttribute(SAMTag.NM.name(), SequenceUtil.calculateSamNmTag(rec, referenceBases, 0, bisulfiteSequence));

                    if (rec.getBaseQualities() != SAMRecord.NULL_QUALS) {
                        rec.setAttribute(SAMTag.UQ.name(), SequenceUtil.sumQualitiesOfMismatches(rec, referenceBases, 0, bisulfiteSequence));
                    }
                }
            }
            writer.addAlignment(rec);
            finalProgress.record(rec);
        }
        writer.close();
        sorted.cleanup();

        CloserUtil.close(unmappedSam);
        logger.info("Wrote " + aligned + " alignment records and " + (alignedReadsOnly ? 0 : unmapped) + " unmapped reads.");
    }

    /**
     * Add record if it is primary or optionally secondary.
     */
    private void addIfNotFiltered(final SortingCollection<SAMRecord> sorted, final SAMRecord rec) {
        if (includeSecondaryAlignments || !rec.isSecondaryAlignment()) {
            sorted.add(rec);
            this.progress.record(rec);
        }
    }

    /**
     * @return Next read's alignment(s) from aligned input or null, if there are no more.
     * The alignments are run through ignoreAlignment() filter before being returned, which may result
     * in an entire read being skipped if all alignments for that read should be ignored.
     */
    private HitsForInsert nextAligned() {
        if (alignedIterator.hasNext()) return alignedIterator.next();
        return null;
    }

    /**
     * Copies alignment info from aligned to unaligned read, clips as appropriate, and sets PG ID.
     *
     * @param unaligned Original SAMRecord, and object into which values are copied.
     * @param aligned   Holds alignment info that will be copied into unaligned.
     */
    private void transferAlignmentInfoToFragment(final SAMRecord unaligned, final SAMRecord aligned) {
        setValuesFromAlignment(unaligned, aligned);
        updateCigarForTrimmedOrClippedBases(unaligned, aligned);
        if (SAMUtils.cigarMapsNoBasesToRef(unaligned.getCigar())) {
            SAMUtils.makeReadUnmapped(unaligned);
        } else if (SAMUtils.recordMapsEntirelyBeyondEndOfReference(aligned)) {
            logger.warn("Record mapped off end of reference; making unmapped: " + aligned);
            SAMUtils.makeReadUnmapped(unaligned);
        }
    }

    /**
     * Copies alignment info from aligned to unaligned read, if there is an alignment, and sets mate information.
     *
     * @param firstUnaligned  Original first of pair, into which alignment and pair info will be written.
     * @param secondUnaligned Original second of pair, into which alignment and pair info will be written.
     * @param firstAligned    Aligned first of pair, or null if no alignment.
     * @param secondAligned   Aligned second of pair, or null if no alignment.
     */
    private void transferAlignmentInfoToPairedRead(final SAMRecord firstUnaligned, final SAMRecord secondUnaligned, final SAMRecord firstAligned, final SAMRecord secondAligned) {
        if (firstAligned != null) transferAlignmentInfoToFragment(firstUnaligned, firstAligned);
        if (secondAligned != null) transferAlignmentInfoToFragment(secondUnaligned, secondAligned);
        if (isClipOverlappingReads()) clipForOverlappingReads(firstUnaligned, secondUnaligned);
        SamPairUtil.setMateInfo(secondUnaligned, firstUnaligned, addMateCigar);
        if (!keepAlignerProperPairFlags) {
            SamPairUtil.setProperPairFlags(secondUnaligned, firstUnaligned, expectedOrientations);
        }
    }


    /**
     * Checks to see whether the ends of the reads overlap and soft clips reads
     * them if necessary.
     */
    protected void clipForOverlappingReads(final SAMRecord read1, final SAMRecord read2) {
        // If both reads are mapped, see if we need to clip the ends due to small
        // insert size
        if (!(read1.getReadUnmappedFlag() || read2.getReadUnmappedFlag())) {

            if (read1.getReadNegativeStrandFlag() != read2.getReadNegativeStrandFlag()) {
                final SAMRecord pos = (read1.getReadNegativeStrandFlag()) ? read2 : read1;
                final SAMRecord neg = (read1.getReadNegativeStrandFlag()) ? read1 : read2;

                // Innies only -- do we need to do anything else about jumping libraries?
                if (pos.getAlignmentStart() < neg.getAlignmentEnd()) {
                    final int posDiff = pos.getAlignmentEnd() - neg.getAlignmentEnd();
                    final int negDiff = pos.getAlignmentStart() - neg.getAlignmentStart();

                    if (posDiff > 0) {
                        CigarUtil.softClip3PrimeEndOfRead(pos, Math.min(pos.getReadLength(),
                                pos.getReadLength() - posDiff + 1));
                    }

                    if (negDiff > 0) {
                        CigarUtil.softClip3PrimeEndOfRead(neg, Math.min(neg.getReadLength(),
                                neg.getReadLength() - negDiff + 1));
                    }

                }
            } else {
                // TODO: What about RR/FF pairs?
            }
        }

    }

    /**
     * Sets the values from the alignment record on the unaligned BAM record.  This
     * preserves all data from the unaligned record (ReadGroup, NoiseRead status, etc)
     * and adds all the alignment info
     *
     * @param rec       The unaligned read record
     * @param alignment The alignment record
     */
    protected void setValuesFromAlignment(final SAMRecord rec, final SAMRecord alignment) {
        for (final SAMRecord.SAMTagAndValue attr : alignment.getAttributes()) {
            // Copy over any non-reserved attributes.  attributesToRemove overrides attributesToRetain.
            if ((!isReservedTag(attr.tag) || this.attributesToRetain.contains(attr.tag)) && !this.attributesToRemove.contains(attr.tag)) {
                rec.setAttribute(attr.tag, attr.value);
            }
        }
        rec.setReadUnmappedFlag(alignment.getReadUnmappedFlag());

        // Note that it is important to get reference names rather than indices in case the sequence dictionaries
        // in the two files are in different orders.
        rec.setReferenceName(alignment.getReferenceName());

        rec.setAlignmentStart(alignment.getAlignmentStart());
        rec.setReadNegativeStrandFlag(alignment.getReadNegativeStrandFlag());
        rec.setSecondaryAlignment(alignment.isSecondaryAlignment());
        rec.setSupplementaryAlignmentFlag(alignment.getSupplementaryAlignmentFlag());
        if (!alignment.getReadUnmappedFlag()) {
            // only aligned reads should have cigar and mapping quality set
            rec.setCigar(alignment.getCigar());  // cigar may change when a
            // clipCigar called below
            rec.setMappingQuality(alignment.getMappingQuality());
        }
        if (rec.getReadPairedFlag()) {
            rec.setProperPairFlag(alignment.getProperPairFlag());
            // Mate info and alignment size will get set by the ClippedPairFixer.
        }

        // If it's on the negative strand, reverse complement the bases
        // and reverse the order of the qualities
        if (rec.getReadNegativeStrandFlag()) {
            rec.reverseComplement(true);
        }

    }

    private static Cigar createNewCigarIfMapsOffEndOfReference(SAMFileHeader header,
                                                               boolean isUnmapped,
                                                               int referenceIndex,
                                                               int alignmentEnd,
                                                               int readLength,
                                                               Cigar oldCigar) {
        Cigar newCigar = null;
        if (!isUnmapped) {
            final SAMSequenceRecord refseq = header.getSequence(referenceIndex);
            final int overhang = alignmentEnd - refseq.getSequenceLength();
            if (overhang > 0) {
                // 1-based index of first base in read to clip.
                int clipFrom = readLength - overhang + 1;
                // we have to check if the last element is soft-clipping, so we can subtract that from clipFrom
                final CigarElement cigarElement = oldCigar.getCigarElement(oldCigar.numCigarElements()-1);
                if (CigarOperator.SOFT_CLIP == cigarElement.getOperator()) clipFrom -= cigarElement.getLength();
                final List<CigarElement> newCigarElements = CigarUtil.softClipEndOfRead(clipFrom, oldCigar.getCigarElements());
                newCigar = new Cigar(newCigarElements);
            }
        }
        return newCigar;
    }

    /**
     * Soft-clip an alignment that hangs off the end of its reference sequence.  Checks both the read and its mate,
     * if available.
     *
     * @param rec
     */
    public static void createNewCigarsIfMapsOffEndOfReference(final SAMRecord rec) {
        // If the read maps off the end of the alignment, clip it
        if (!rec.getReadUnmappedFlag()) {
            final Cigar readCigar = createNewCigarIfMapsOffEndOfReference(rec.getHeader(),
                    rec.getReadUnmappedFlag(),
                    rec.getReferenceIndex(),
                    rec.getAlignmentEnd(),
                    rec.getReadLength(),
                    rec.getCigar());
            if (null != readCigar) rec.setCigar(readCigar);
        }

        // If the read's mate maps off the end of the alignment, clip it
        if (SAMUtils.hasMateCigar(rec)) {
            Cigar mateCigar = SAMUtils.getMateCigar(rec);
            mateCigar = createNewCigarIfMapsOffEndOfReference(rec.getHeader(),
                    rec.getMateUnmappedFlag(),
                    rec.getMateReferenceIndex(),
                    SAMUtils.getMateAlignmentEnd(rec), // NB: this could be computed without another call to getMateCigar
                    mateCigar.getReadLength(),
                    mateCigar);
            if (null != mateCigar) rec.setAttribute(SAMTag.MC.name(), mateCigar.toString());
        }
    }

    protected void updateCigarForTrimmedOrClippedBases(final SAMRecord rec, final SAMRecord alignment) {
        // If the read was trimmed or not all the bases were sent for alignment, clip it
        final int alignmentReadLength = alignment.getReadLength();
        final int originalReadLength = rec.getReadLength();
        final int trimmed = (!rec.getReadPairedFlag()) || rec.getFirstOfPairFlag()
                ? this.read1BasesTrimmed != null ? this.read1BasesTrimmed : 0
                : this.read2BasesTrimmed != null ? this.read2BasesTrimmed : 0;
        final int notWritten = originalReadLength - (alignmentReadLength + trimmed);

        // Update cigar if the mate maps off the reference
        createNewCigarsIfMapsOffEndOfReference(rec);

        rec.setCigar(CigarUtil.addSoftClippedBasesToEndsOfCigar(
                rec.getCigar(), rec.getReadNegativeStrandFlag(), notWritten, trimmed));

        // If the adapter sequence is marked and clipAdapter is true, clip it
        if (this.clipAdapters && rec.getAttribute(ReservedTagConstants.XT) != null) {
            CigarUtil.softClip3PrimeEndOfRead(rec, rec.getIntegerAttribute(ReservedTagConstants.XT));
        }
    }


    protected SAMProgramRecord getProgramRecord() { return programRecord; }

    protected void setProgramRecord(final SAMProgramRecord pg) {
        Utils.validate(programRecord == null, "Cannot set program record more than once on alignment merger.");
        programRecord = pg;
        header.addProgramRecord(pg);
        SAMUtils.chainSAMProgramRecord(header, pg);
    }

    protected boolean isReservedTag(final String tag) {
        final char firstCharOfTag = tag.charAt(0);

        // All tags that start with a lower-case letter are user defined and should not be overridden by aligner
        // unless explicitly specified in attributesToRetain.
        if (Character.isLowerCase(firstCharOfTag)) return true;

        for (final char c : RESERVED_ATTRIBUTE_STARTS) {
            if (firstCharOfTag == c) return true;
        }
        return false;
    }

    protected void resetRefSeqFileWalker() {
        this.refSeq = new ReferenceSequenceFileWalker(referenceFasta);
    }

    public boolean isClipOverlappingReads() {
        return clipOverlappingReads;
    }

    public void setClipOverlappingReads(final boolean clipOverlappingReads) {
        this.clipOverlappingReads = clipOverlappingReads;
    }

    /**
     * If true, keep the aligner's idea of proper pairs rather than letting alignment merger decide.
     */
    public void setKeepAlignerProperPairFlags(final boolean keepAlignerProperPairFlags) {
        this.keepAlignerProperPairFlags = keepAlignerProperPairFlags;
    }

    public void setIncludeSecondaryAlignments(final boolean includeSecondaryAlignments) {
        this.includeSecondaryAlignments = includeSecondaryAlignments;
    }

    public void close() {
        CloserUtil.close(this.refSeq);
    }
}
