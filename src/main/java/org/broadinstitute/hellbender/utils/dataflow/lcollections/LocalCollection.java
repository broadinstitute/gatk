package org.broadinstitute.hellbender.utils.dataflow.lcollections;

import com.google.common.collect.Lists;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.function.ToIntFunction;

public class LocalCollection<T> implements Serializable {
    private static final long serialVersionUID = 1l;
    ArrayList<T> data;

    /**
     * Creates a new LocalCollection that holds pointers to these values.
     * Do not mutate the values.
     */
    public static <T> LocalCollection<T> of(Iterable<T> values) {
        return new LocalCollection<>(Lists.newArrayList(values));
    }

    /**
     * Creates a new LocalCollection that holds the union of the contents
     * of the collections passed in argument.
     */
    public static <T> LocalCollection<T> union(LocalCollection<T>... lcs) {
        ArrayList<T> data = new ArrayList<T>();
        for (LocalCollection<T> c : lcs) {
            data.addAll(c.data);
        }
        return new LocalCollection<T>(data);
    }

    LocalCollection() {
        data = new ArrayList<>();
    }

    // don't change ts after passing it in.
    LocalCollection(ArrayList<T> ts) {
        data = ts;
    }

    public ArrayList<LocalCollection<T>> partitionBy(ToIntFunction<T> keyer) {
        return partitionBy(keyer, 0);
    }

    public ArrayList<LocalCollection<T>> partitionBy(ToIntFunction<T> keyer, int minSize) {
        ArrayList<LocalCollection<T>> ret = new ArrayList<LocalCollection<T>>(minSize);
        while (ret.size() < minSize) {
            ret.add(new LocalCollection<T>());
        }
        for (T el : data) {
            int key = keyer.applyAsInt(el);
            while (ret.size() <= key) {
                ret.add(new LocalCollection<T>());
            }
            // we're mutating the LocalCollection but that's OK because it's not closed yet.
            ret.get(key).data.add(el);
        }
        // no more mutation from this point on.
        return ret;
    }

    public LocalGroupedCollection<T> groupBy(Function<T,String> keyer) {
        return new LocalGroupedCollection<>(data,keyer);
    }

    // apply a function to every element of the collection
    public <U> LocalCollection<U> map(Function<T,U> func) {
        ArrayList<U> results = new ArrayList<U>(data.size());
        for (T el : data) {
            results.add(func.apply(el));
        }
        return new LocalCollection<>(results);
    }

    // apply a function to the whole collection
    public <U> LocalCollection<U> transform(Function<Iterable<T>, ArrayList<U>> func) {
        return new LocalCollection<>(func.apply(data));
    }

    // return a new collection with only the elements that satisfy the given predicate.
    public LocalCollection<T> filter(Predicate<T> filterPred) {
        ArrayList<T> results = new ArrayList<T>(data.size()/4);
        for (T el : data) {
            if (filterPred.test(el)) {
                results.add(el);
            }
        }
        return new LocalCollection<>(results);
    }

    public Iterable<T> iterable() {
        return data;
    }

    // number of elements in the collection. It's local, after all.
    public int size() {
        return data.size();
    }
}
