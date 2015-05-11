package org.broadinstitute.hellbender.utils.read.markduplicates;

import htsjdk.samtools.CoordinateSortedPairInfoMap;
import org.broadinstitute.hellbender.exceptions.GATKException;


import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.AbstractMap;
import java.util.Map;

/**
 * Disk-based implementation of ReadEndsForMarkDuplicatesMap.  A subdirectory of the system tmpdir is created to store
 * files, one for each reference sequence.  The reference sequence that is currently being queried (i.e. the
 * sequence for which remove() has been most recently called) is stored in RAM.  ReadEnds for all other sequences
 * are stored on disk.
 * <p/>
 * When put() is called for a sequence that is the current one in RAM, the ReadEnds object is merely put into the
 * in-memory map.  If put() is called for a sequence ID that is not the current RAM one, the ReadEnds object is
 * appended to the file for that sequence, creating the file if necessary.
 * <p/>
 * When remove() is called for a sequence that is the current one in RAM, remove() is called on the in-memory map.
 * If remove() is called for a sequence other than the current RAM sequence, then the current RAM sequence is written
 * to disk, the new sequence is read from disk into RAM map, and the file for the new sequence is deleted.
 * <p/>
 * If things work properly, and reads are processed in genomic order, records will be written for mates that are in
 * a later sequence.  When the mate is reached in the input SAM file, the file that was written will be deleted.
 * This should result in all temporary files being deleted by the time all the reads are processed.  The temp
 * directory is marked to be deleted on exit so everything should get cleaned up.
 *
 * @author alecw@broadinstitute.org
 */
public final class DiskBasedReadEndsForMarkDuplicatesMap implements ReadEndsForMarkDuplicatesMap {
    private final CoordinateSortedPairInfoMap<String, ReadEndsForMarkDuplicates> pairInfoMap;

    public DiskBasedReadEndsForMarkDuplicatesMap(int maxOpenFiles) {
        pairInfoMap = new CoordinateSortedPairInfoMap<>(maxOpenFiles, new Codec());
    }

    public ReadEndsForMarkDuplicates remove(int mateSequenceIndex, String key) {
        return pairInfoMap.remove(mateSequenceIndex, key);
    }

    public void put(int mateSequenceIndex, String key, ReadEndsForMarkDuplicates readEnds) {
        pairInfoMap.put(mateSequenceIndex, key, readEnds);
    }

    public int size() {
        return pairInfoMap.size();
    }

    public int sizeInRam() {
        return pairInfoMap.sizeInRam();
    }

    private static class Codec implements CoordinateSortedPairInfoMap.Codec<String, ReadEndsForMarkDuplicates> {
        private final ReadEndsForMarkDuplicatesCodec readEndsForMarkDuplicatesCodec = new ReadEndsForMarkDuplicatesCodec();

        public void setInputStream(final InputStream is) {
            readEndsForMarkDuplicatesCodec.setInputStream(is);
        }

        public void setOutputStream(final OutputStream os) {
            readEndsForMarkDuplicatesCodec.setOutputStream(os);
        }

        public Map.Entry<String, ReadEndsForMarkDuplicates> decode() {
            try {
                final String key = readEndsForMarkDuplicatesCodec.getInputStream().readUTF();
                final ReadEndsForMarkDuplicates record = readEndsForMarkDuplicatesCodec.decode();
                return new AbstractMap.SimpleEntry<>(key, record);
            } catch (IOException e) {
                throw new GATKException("Error loading ReadEndsForMarkDuplicatesMap from disk", e);
            }
        }

        public void encode(final String key, final ReadEndsForMarkDuplicates readEnds) {
            try {
                readEndsForMarkDuplicatesCodec.getOutputStream().writeUTF(key);
                readEndsForMarkDuplicatesCodec.encode(readEnds);
            } catch (IOException e) {
                throw new GATKException("Error spilling ReadEndsForMarkDuplicatesMap to disk.", e);
            }
        }
    }

}
