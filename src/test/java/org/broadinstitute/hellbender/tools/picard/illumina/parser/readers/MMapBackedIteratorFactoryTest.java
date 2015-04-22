package org.broadinstitute.hellbender.tools.picard.illumina.parser.readers;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.nio.Buffer;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.IntBuffer;
import java.util.*;

public final class MMapBackedIteratorFactoryTest {
    public static File TestDataDir = new File("src/test/resources/org/broadinstitute/hellbender/tools/picard/illumina/readerTests");
    public static File BinFile     = new File(TestDataDir, "binary_passing.bin");

    public static final int FileLength = 51;
    //remember that the MMappedBinaryFileReader assumes little endianness
    public byte[] fileAsBytes(final int start, final int end){
         final int [] bInts = {
                         0x31, 0x22, 0x41, 0x01,    0x45, 0x6E, 0x64, 0x4F,
                         0x66, 0x48, 0x65, 0x61,    0x64, 0x65, 0x72, 0x42,
                         0x6F, 0x64, 0x79, 0x50,    0x61, 0x72, 0x74, 0x6F,
                         0x66, 0x54, 0x68, 0x65,    0x46, 0x69, 0x6C, 0x65,
                         0x37, 0x37, 0x0A, 0x45,    0x6E, 0x64, 0x43, 0x6F,
                         0x6D, 0x6D, 0x75, 0x6E,    0x69, 0x63, 0x61, 0x74,
                         0x69, 0x6F, 0x6E
        };

        final int total = end - start + 1;
        final byte [] bytes = new byte[total];
        for(int i = 0; i < total; i++) {
            bytes[i] = (byte)bInts[start + i];
        }
        return bytes;
    }

    //Note these integers are different from the ones in fileAsBytes->bInts (because we're talking about 4 bytes at a time not one)
    @DataProvider(name="passing_bin_asTdBuffer")
    public ByteBuffer headerAsByteBuffer(final int headerSize) {
        final byte [] bytes = fileAsBytes(0, headerSize-1);
        final ByteBuffer bb = ByteBuffer.allocate(bytes.length);
        bb.put(bytes);
        bb.position(0);
        bb.order(ByteOrder.LITTLE_ENDIAN);
        return bb;
    }

    @DataProvider(name="passing_bin_asTdBuffer")
    public ByteBuffer fileAsByteBuffer(final int headerSize) {
        final byte [] bytes = fileAsBytes(headerSize, FileLength-1);
        final ByteBuffer bb = ByteBuffer.allocate(bytes.length);
        bb.put(bytes);
        bb.position(0);
        bb.order(ByteOrder.LITTLE_ENDIAN);
        return bb;
    }

    abstract class FileTestDef<T> {
        public final int headerSize;
        public final BinaryFileIterator bbIter;
        public final int numElements;

        public FileTestDef(final int headerSize, final BinaryFileIterator<T> bbIter) {
            this.headerSize  = headerSize;
            this.bbIter = bbIter;
            this.numElements = fileAsBytes(headerSize, FileLength-1).length / bbIter.getElementSize();
        }

        public void test() {
            final ByteBuffer testBuffer = fileAsByteBuffer(headerSize);


            if(headerSize > 0) {
                final ByteBuffer headerBuffer = headerAsByteBuffer(headerSize);
                testHeaderBytes(headerBuffer, bbIter.getHeaderBytes());
            }

            bbIter.assertTotalElementsEqual(numElements);
            final Iterator<T> testIter = getTestIterator(testBuffer);

            while(hasNext(testIter, bbIter)) {
                Assert.assertEquals(testIter.next(), bbIter.next());
            }
        }

        public abstract Iterator<T> getTestIterator(final ByteBuffer byteBuffer);

        public boolean hasNext(final Iterator<T> testIter, final Iterator<T> fileIter) {
            if(testIter.hasNext() && fileIter.hasNext())
                return true;

            if(testIter.hasNext()) {
                throw new RuntimeException("Test data (testIter) has more iterations while fileIter does not!");
            }

            if(fileIter.hasNext()) {
                throw new RuntimeException("File data (fileIter) has more iterations while testIter does not!");
            }

            return false;
        }
    }

    abstract class NoHeaderTestIter<T> implements Iterator<T> {
        public final Buffer buf;

        public NoHeaderTestIter(final Buffer buf) {
            this.buf = buf;
        }

        public boolean hasNext() {
            return buf.hasRemaining();
        }

        public void remove() {
            throw new UnsupportedOperationException();
        }
    }

    @DataProvider(name = "validTestDefs")
    public Object[][] validTestDefs() {
        return new Object[][] {
            {
                new FileTestDef<Integer>(15, MMapBackedIteratorFactory.getIntegerIterator(15, BinFile)) {
                    @Override
                    public Iterator<Integer> getTestIterator(final ByteBuffer byteBuffer) {
                        final IntBuffer ib = byteBuffer.asIntBuffer();
                        return new NoHeaderTestIter<Integer>(ib) {
                            @Override
                            public Integer next() {
                                return ib.get();
                            }
                        };
                    }
                }
            },
            {
                new FileTestDef<Byte>(2, MMapBackedIteratorFactory.getByteIterator(2, BinFile)) {
                    @Override
                    public Iterator<Byte> getTestIterator(final ByteBuffer byteBuffer) {
                        return new NoHeaderTestIter<Byte>(byteBuffer) {
                            @Override
                            public Byte next() {
                                return byteBuffer.get();
                            }
                        };
                    }
                }
            },
            {
                new FileTestDef<Float>(19, MMapBackedIteratorFactory.getFloatIterator(19, BinFile)) {
                    @Override
                    public Iterator<Float> getTestIterator(final ByteBuffer byteBuffer) {
                        return new NoHeaderTestIter<Float>(byteBuffer) {
                            @Override
                            public Float next() {
                                return byteBuffer.getFloat();
                            }
                        };
                    }
                }
            }
        };
    }

    @Test(dataProvider="validTestDefs")
    public void testValidConfigurations(final FileTestDef ftd) {
        ftd.test();
    }

    @Test
    public void onlyHeaderTest() {
        final BinaryFileIterator<Integer> bbIter = MMapBackedIteratorFactory.getIntegerIterator((int)BinFile.length(), BinFile);
        Assert.assertEquals(bbIter.getHeaderBytes(), headerAsByteBuffer((int) BinFile.length()));
        Assert.assertFalse(bbIter.hasNext());
    }

    @Test(expectedExceptions = IlluminaReaderException.class)
    public void tooLargeHeaderTest() {
        final BinaryFileIterator<Integer> bbIter = MMapBackedIteratorFactory.getIntegerIterator(FileLength + 10, BinFile);
        bbIter.getHeaderBytes();
    }

    @Test(expectedExceptions = IlluminaReaderException.class)
    public void negativeHeaderTest() {
        final BinaryFileIterator<Integer> bbIter = MMapBackedIteratorFactory.getIntegerIterator(-1, BinFile);
        bbIter.getHeaderBytes();
    }

    @DataProvider(name="invalidFileSizes")
    public Object[][] invalidFileSizes() {
        return new Object[][]{
            {1, 12}, //should result in two left over bytes
            {3, 11}, //should result in one left over elements
        };
    }

    @Test(expectedExceptions = IlluminaReaderException.class, dataProvider = "invalidFileSizes")
    public void invalidFileSizeTests(final int headerSize, final int expectedElements) {
        final BinaryFileIterator<Integer> bbIter = MMapBackedIteratorFactory.getIntegerIterator(headerSize, BinFile);
        bbIter.assertTotalElementsEqual(expectedElements);
    }

    public void testHeaderBytes(final ByteBuffer bb1, final ByteBuffer bb2) {
        Assert.assertTrue(bb1.equals(bb2), "Header bytes are not equal! " + bb1.toString() + "  !=  " + bb2.toString());
    }
}
