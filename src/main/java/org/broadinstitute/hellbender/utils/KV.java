package org.broadinstitute.hellbender.utils;

import java.util.Objects;

/**
 * replacement for dataflow Key-Value class, don't use this anywhere new
 */
public class KV<K, V> {

    private final K key;
    private final V value;

    public static <K, V>  KV<K, V> of(K key, V value){
        return new KV<>(key, value);
    }

    private KV(K key, V value){
        this.key = key;
        this.value = value;
    }


    public K getKey() {
        return key;
    }

    public V getValue() {
        return value;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        KV<?, ?> kv = (KV<?, ?>) o;
        return Objects.equals(key, kv.key) &&
                Objects.equals(value, kv.value);
    }

    @Override
    public int hashCode() {
        return Objects.hash(key, value);
    }

    @Override
    public String toString() {
        return "(" + key + ", " + value + ")";
    }
}
