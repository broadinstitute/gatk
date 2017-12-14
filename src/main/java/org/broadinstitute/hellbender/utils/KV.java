package org.broadinstitute.hellbender.utils;

/**
 * replacement for dataflow Key-Value class
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
}
