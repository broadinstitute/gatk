package org.broadinstitute.hellbender.engine.spark;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.Serializer;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import com.google.api.client.googleapis.util.Utils;
import com.google.api.client.json.JsonFactory;
import org.broadinstitute.hellbender.exceptions.GATKException;

import java.io.IOException;

/**
 * This is a simple JsonSerializer for use with Kryo. We need this because the Google
 * Genomics classes don't work "out of the box" with Kryo.
 * @param <T> The type of Json-like class.
 */
public class JsonSerializer<T> extends Serializer<T> {
    private static final JsonFactory JSON_FACTORY = Utils.getDefaultJsonFactory();

    @Override
    public void write(Kryo kryo, Output output, T object) {
        try {
            output.writeString(JSON_FACTORY.toString(object));
        } catch (IOException e) {
            throw new GATKException("Json writing error:" + e);
        }
    }

    @Override
    public T read(Kryo kryo, Input input, Class<T> type) {
        String s = input.readString();
        try {
            return JSON_FACTORY.fromString(s, type);
        } catch (IOException e) {
            throw new GATKException("Json reading error" + e);
        }
    }
}