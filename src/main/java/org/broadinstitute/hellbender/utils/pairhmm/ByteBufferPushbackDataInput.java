package org.broadinstitute.hellbender.utils.pairhmm;

import org.broadinstitute.hellbender.utils.Utils;

import java.io.*;
import java.nio.BufferOverflowException;
import java.nio.ByteBuffer;

public class ByteBufferPushbackDataInput implements PushbackDataInput {

    private final byte[] bufferArray;
    private final ByteBuffer buffer;

    private final DataInput input;

    public ByteBufferPushbackDataInput(final DataInput in, final int bufferSize) {
        input = Utils.nonNull(in);
        bufferArray = new byte[bufferSize];
        buffer = ByteBuffer.wrap(bufferArray);
    }

    @Override
    public void unread(final byte b) throws IOException {
        try {
            buffer.put(b);
        } catch (final BufferOverflowException ex) {
            throw new IOException("buffer is full", ex);
        }
    }

    @Override
    public void unread(final byte[] bytes, final int offset, final int len) throws IOException {
        try {
            buffer.put(bytes, offset, len);
        } catch (final BufferOverflowException ex) {
            throw new IOException("buffer is full", ex);
        }
    }

    @Override
    public void unread(final short s) throws IOException {
        try {
            buffer.putShort(s);
        } catch (final BufferOverflowException ex) {
            throw new IOException("buffer is full", ex);
        }
    }

    @Override
    public void unread(final long l) throws IOException {
        try {
            buffer.putLong(l);
        } catch (final BufferOverflowException ex) {
            throw new IOException("buffer is full", ex);
        }
    }

    @Override
    public void unread(final int i) throws IOException {
        try {
            buffer.putInt(i);
        } catch (final BufferOverflowException ex) {
            throw new IOException("buffer is full", ex);
        }
    }

    @Override
    public void unread(final double d) throws IOException {
        try {
            buffer.putDouble(d);
        } catch (final BufferOverflowException ex) {
            throw new IOException("buffer is full", ex);
        }
    }

    @Override
    public void unread(final float f) throws IOException {
        try {
            buffer.putFloat(f);
        } catch (final BufferOverflowException ex) {
            throw new IOException("buffer is full", ex);
        }
    }

    @Override
    public void unread(final boolean b) throws IOException {
        try {
            buffer.put((byte) (b ? 0 : 1));
        } catch (final BufferOverflowException ex) {
            throw new IOException("buffer is full", ex);
        }
    }

    @Override
    public void unread(final char c) throws IOException {
        try {
            buffer.putChar(c);
        } catch (final BufferOverflowException ex) {
            throw new IOException("buffer is full", ex);
        }

    }

    @Override
    public void unreadLine(final String str) throws IOException {
        unread('\n');
        for (int i = str.length() - 1; i >= 0; i++) {
            unread(str.charAt(i));
        }
    }


    @Override
    public void readFully(final byte[] b) throws IOException {
        readFully(Utils.nonNull(b), 0, b.length);
    }

    @Override
    public void readFully(final byte[] b, final int off, final int len) throws IOException {
        Utils.nonNull(b);
        final int position = buffer.position();
        if (position == 0) {
            input.readFully(b, off, len);
        } else if (position >= len) {
            final int offset = position - len;
            Utils.flip(bufferArray, offset, position);
            System.arraycopy(bufferArray, offset, b, off, len);
            buffer.position(offset);
        } else {
            System.arraycopy(bufferArray, 0, b, off, position);
            input.readFully(b, off + position, len - position);
            buffer.clear();
        }
    }

    @Override
    public int skipBytes(final int n) throws IOException {
        if (n <= 0) {
            return 0;
        }
        final int position = buffer.position();
        if (position == 0) {
            return input.skipBytes(n);
        } else if (position >= n) {
            buffer.position(position - n);
            return n;
        } else {
            buffer.clear();
            return position + input.skipBytes(n - position);
        }
    }

    @Override
    public boolean readBoolean() throws IOException {
        final int position = buffer.position();
        if (position == 0) {
            return input.readBoolean();
        } else {
            final int newPosition = position - 1;
            final boolean result = buffer.get(newPosition) != 0;
            buffer.position(newPosition);
            return result;
        }
    }

    @Override
    public byte readByte() throws IOException {
        final int position = buffer.position();
        if (position == 0) {
            return input.readByte();
        } else {
            final int newPosition = position - 1;
            final byte result = buffer.get(newPosition);
            buffer.position(newPosition);
            return result;
        }
    }

    @Override
    public int readUnsignedByte() throws IOException {
        return 0xFF & readByte();
    }

    @Override
    public short readShort() throws IOException {
        final int position = buffer.position();
        if (position == 0) {
            return input.readShort();
        } else if (position == 1) {
            buffer.put(input.readByte());
        }
        final int newPosition = position - 1;
        final short result = buffer.getShort(newPosition);
        buffer.position(newPosition);
        return result;
    }

    @Override
    public int readUnsignedShort() throws IOException {
        return 0xFFFF & readShort();
    }

    @Override
    public char readChar() throws IOException {
        final int position = buffer.position();
        if (position == 0) {
            return input.readChar();
        } else if (position == 1) {
            buffer.put(input.readByte());
        }
        final int newPosition = position - 1;
        final char result = buffer.getChar(newPosition);
        buffer.position(newPosition);
        return result;
    }

    @Override
    public int readInt() throws IOException {
        final int position = buffer.position();
        if (position == 0) {
            return input.readInt();
        }
        final int newPosition = ensureBufferHasBytes(position,4);
        final int result = buffer.getInt(newPosition);
        buffer.position(newPosition);
        return result;
    }

    @Override
    public long readLong() throws IOException {
        final int position = buffer.position();
        if (position == 0) {
            return input.readLong();
        }
        final int newPosition = ensureBufferHasBytes(position,8);
        final long result = buffer.getLong(newPosition);
        buffer.position(newPosition);
        return result;
    }

    private int ensureBufferHasBytes(final int position, final int nb) throws IOException {
        if (position < nb) {
            try {
                input.readFully(bufferArray, position, nb - position);
            } catch (final IndexOutOfBoundsException ex) {
                throw new IOException("buffer is full", new BufferOverflowException());
            }
            buffer.position(nb);
            return 0;
        } else {
            return position - nb;
        }
    }

    @Override
    public float readFloat() throws IOException {
        final int position = buffer.position();
        if (position == 0) {
            return input.readFloat();
        }
        final int newPosition = ensureBufferHasBytes(position, 4);
        final float result = buffer.getFloat(newPosition);
        buffer.position(newPosition);
        return result;
    }

    @Override
    public double readDouble() throws IOException {
        final int position = buffer.position();
        if (position == 0) {
            return input.readFloat();
        }
        final int newPosition = ensureBufferHasBytes(position, 8);
        final double result = buffer.getDouble(newPosition);
        buffer.position(newPosition);
        return result;
    }

    @Override
    public String readLine() throws IOException {
        final int position = buffer.position();
        if (position == 0) {
            return input.readLine();
        } else {
            String result = null;
            int i;
            for (i = position - 1; i >= 0; i--) {
                byte b = bufferArray[i];
                if (b == '\n') {
                    Utils.flip(bufferArray, i + 1, position);
                    result = new String(bufferArray, i + 1, position - i - 1);
                    buffer.position(i);
                    break;
                } else if (b == '\r') {
                    Utils.flip(bufferArray, i + 1, position);
                    result = new String(bufferArray, i + 1, position - i - 1);
                    if (i > 0) {
                        buffer.position(bufferArray[i - 1] == '\n' ? i - 1 : i);
                    } else {
                        buffer.position(0);
                        try {
                            final byte next = input.readByte();
                            if (next != '\n') {
                                buffer.put(next);
                            }
                        } catch (final EOFException ex) {
                            // nothing to do.
                        }
                    }
                    break;
                }
            }
            if (result == null) {
                Utils.flip(bufferArray, 0, position);
                final String head = new String(bufferArray, 0, position);
                final String tail = input.readLine();
                result = head + tail;
                buffer.position(0);
            }
            return result;
        }
    }

    public void unreadUTF(final String str) throws IOException {
        final int len = str.length();
        final byte[] ba = new byte[len * 3];
        int i = 0;
        for (int j = 0; j < len ; j++) {
            final char ch = str.charAt(j);
            if ((ch & 0xFFFFFF80) == 0) {
                ba[i++] = ((byte) ch);
            } else if ((ch & 0xFFFFF800) == 0) {
                ba[i++] = ((byte) (192 | (ch & 0x07C0) >> 6));
                ba[i++] = ((byte) (128 | (ch & 0x3F)));
            } else {
                ba[i++] = ((byte) (224 | (ch & 0xF000) >> 12));
                ba[i++] = ((byte) (128 | (ch & 0xFC0) >> 6));
                ba[i++] = ((byte) (128 | (ch & 0x3F)));
            }
        }
        unread(ba, 0, i);
        unreadUnsignedShort(i);
    }

    @Override
    public String readUTF() throws IOException {
        final int length = readUnsignedShort();
        final byte[] body = new byte[length];
        readFully(body);
        final StringBuilder builder = new StringBuilder(length);
        for (int i = 0; i < body.length; i++) {
            final byte a = body[i];
            if ((a & 128) == 0) {
                builder.append((char) (a & 0xFF));
            } else if ((a & 192) == 192) {
                if (++i >= body.length) {
                    throw new UTFDataFormatException();
                }
                final byte b = body[i];
                if ((b & 192) == 128) {
                    builder.append((char)(((a & 0x1F) << 6) | (b & 0x3F)));
                } else {
                    throw new UTFDataFormatException();
                }
            } else if ((a & 0xF0) == 224) {
                if (i + 2 >= body.length) {
                    throw new UTFDataFormatException();
                }
                final byte b = body[++i];
                final byte c = body[++i];
                if ((b & 192) == 128 & (c & 192) == 128) {
                    builder.append((char)(((a & 0x0F) << 12) | ((b & 0x3F) << 6) | (c & 0x3F)));
                } else {
                    throw new UTFDataFormatException();
                }
            } else {
                throw new UTFDataFormatException();
            }
        }
        return builder.toString();
    }

    public boolean eof() throws IOException {
        if (buffer.position() > 0) {
            return false;
        }
        try {
            unread(readByte());
            return false;
        } catch (final EOFException ex) {
            return true;
        }
    }
}
