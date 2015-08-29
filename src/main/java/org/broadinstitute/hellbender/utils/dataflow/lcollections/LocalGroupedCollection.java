package org.broadinstitute.hellbender.utils.dataflow.lcollections;


import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.function.Function;

public class LocalGroupedCollection<T> implements Serializable {
    private static final long serialVersionUID = 1l;

    HashMap<String, LocalCollection<T>> groups;

    public LocalGroupedCollection(ArrayList<T> data, Function<T,String> keyer) {
        groups = new HashMap<>();
        for (T el : data) {
            String key = keyer.apply(el);
            LocalCollection<T> c;
            c = groups.getOrDefault(key, null);
            if (c==null) {
                c = new LocalCollection<>();
                groups.put(key, c);
            }
            c.data.add(el);
        }
    }

    private LocalGroupedCollection(HashMap<String, LocalCollection<T>> data) {
        groups = data;
    }

    // get a particular group
    public LocalCollection<T> get(String key) {
        return groups.get(key);
    }

    // apply a function to each group
    public <U> LocalGroupedCollection<U> map(Function<LocalCollection<T>, LocalCollection<U>> mapper) {
        HashMap<String, LocalCollection<U>> groups2 = new HashMap<String, LocalCollection<U>>();
        for (String k : groups.keySet()) {
            groups2.put(k, mapper.apply(groups.get(k)) );
        }
        return new LocalGroupedCollection<>(groups2);
    }

    // map and transform together
    public <U> LocalGroupedCollection<U> mapTransform(Function<Iterable<T>, ArrayList<U>> transformFn) {
        return map( (collection) -> collection.transform(transformFn) );
    }

    // return a single collection that contains the union of all the groups' content.
    public LocalCollection<T> flatten() {
        ArrayList<T> data = new ArrayList<>();
        for (LocalCollection<T> c : groups.values()) {
            data.addAll(c.data);
        }
        return new LocalCollection<>(data);
    }
}
