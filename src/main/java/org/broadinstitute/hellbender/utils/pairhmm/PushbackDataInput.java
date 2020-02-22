package org.broadinstitute.hellbender.utils.pairhmm;

import org.broadinstitute.hellbender.utils.Utils;

import java.io.DataInput;
import java.io.EOFException;
import java.io.IOException;

public interface PushbackDataInput extends DataInput {

    void unread(byte b) throws IOException;

    default void unread(byte[] bytes) throws IOException {
        unread(Utils.nonNull(bytes), 0, bytes.length);
    }

    void unread(byte[] bytes, int offset, int len) throws IOException;

    void unread(short s) throws IOException;

    void unread(long l) throws IOException;

    void unread(int i) throws IOException;

    void unread(double d) throws IOException;

    void unread(float f) throws IOException;

    void unread(boolean b) throws IOException;

    default void unreadUnsignedByte(int b) throws IOException {
        if (b != (b & 0xFF)) {
            throw new IllegalArgumentException();
        }
        unread((byte) b);
    }

    default void unreadUnsignedShort(int s) throws IOException {
        if (s != (s & 0xFFFF)) {
            throw new IllegalArgumentException("unsigned short overflow");
        }
        unread((short) s);
    }

    default boolean eof() throws IOException {
        try {
            unread(readByte());
            return false;
        } catch (final EOFException ex) {
            return true;
        }
    }

    void unread(char c) throws IOException;

    void unreadLine(String str) throws IOException;

    void unreadUTF(String utf) throws IOException;


}
