package org.broadinstitute.hellbender.utils.diffengine;

import htsjdk.samtools.*;
import htsjdk.samtools.util.BlockCompressedInputStream;
import org.broadinstitute.hellbender.exceptions.GATKException;

import java.io.*;
import java.util.Arrays;

/**
 * Class implementing diffnode reader for BAM/SAM/CRAM files.
 */
final class BAMDiffableReader implements DiffableReader {

    private static final byte[] BAM_MAGIC = "BAM\1".getBytes();

    @Override
    public String getName() { return "BAM"; }

    @Override
    public DiffElement readFromFile(final File file, final int maxElementsToRead) throws IOException{
        try(final SamReader reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(file)) {

            final DiffNode root = DiffNode.rooted(file.getName());
            final SAMRecordIterator iterator = reader.iterator();

            int count = 0;
            while (iterator.hasNext()) {
                final SAMRecord read = iterator.next();

                // name is the read name + first of pair
                String name = read.getReadName().replace('.', '_');
                if (read.getReadPairedFlag()) {
                    name += read.getFirstOfPairFlag() ? "_1" : "_2";
                }

                final DiffNode readRoot = DiffNode.empty(name, root);

                // add fields
                readRoot.add("NAME", read.getReadName());
                readRoot.add("FLAGS", read.getFlags());
                readRoot.add("RNAME", read.getReferenceName());
                readRoot.add("POS", read.getAlignmentStart());
                readRoot.add("MAPQ", read.getMappingQuality());
                readRoot.add("CIGAR", read.getCigarString());
                readRoot.add("RNEXT", read.getMateReferenceName());
                readRoot.add("PNEXT", read.getMateAlignmentStart());
                readRoot.add("TLEN", read.getInferredInsertSize());
                readRoot.add("SEQ", read.getReadString());
                readRoot.add("QUAL", read.getBaseQualityString());

                for (final SAMRecord.SAMTagAndValue xt : read.getAttributes()) {
                    readRoot.add(xt.tag, xt.value);
                }

                if (!root.hasElement(name)) { // protect ourselves from malformed files
                    root.add(readRoot);
                }
                count += readRoot.size();
                if (count > maxElementsToRead && maxElementsToRead != -1) {
                    break;
                }
            }
            return root.getBinding();
        }
    }

    @Override
    public boolean canRead(final File file) {
        final byte[] buffer = new byte[BAM_MAGIC.length];
        try (final InputStream fstream = new BufferedInputStream(new FileInputStream(file))){
            if ( !BlockCompressedInputStream.isValidFile(fstream) ) {
                return false;
            }
            final BlockCompressedInputStream bcis = new BlockCompressedInputStream(fstream);
            bcis.read(buffer, 0, BAM_MAGIC.length);
            return Arrays.equals(buffer, BAM_MAGIC);
        } catch ( IOException | FileTruncatedException e ) {
            return false;
        }
    }
}
