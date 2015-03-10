package org.broadinstitute.hellbender.utils.iterators;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;
/**
 * <p/>
 * Interface GATKSAMIterator
 * <p/>
 * This is the standard interface for all iterators in the GATK package that iterate over SAMRecords
 */
public interface GATKSAMIterator extends CloseableIterator<SAMRecord>, Iterable<SAMRecord> {
}
