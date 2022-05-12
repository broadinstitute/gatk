package org.broadinstitute.hellbender.tools.walkers.groundtruth;

import htsjdk.samtools.util.Tuple;
import org.broadinstitute.hellbender.engine.GATKPath;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

public class SingleFileLocationTranslator {

    final private int[]           pos;
    final private int[]           offset;

    SingleFileLocationTranslator(final GATKPath path) throws IOException {

        // read the file in. we assume it is sorted and starts with pos=1
        final List<Tuple<Integer, Integer>>       data = new LinkedList<>();
        final BufferedReader      reader = new BufferedReader(new InputStreamReader(path.getInputStream()));
        String                    line = reader.readLine();       // ignore first line
        while ( (line = reader.readLine()) != null ) {
            final String[]        toks = line.split(",");
            data.add(new Tuple<>(Integer.parseInt(toks[0]), Integer.parseInt(toks[1])));
        }
        pos = data.stream().map(t->t.a).mapToInt(i->i).toArray();
        offset = data.stream().map(t->t.b).mapToInt(i->i).toArray();
    }

    int translate(final int from) {

        // search for starting point
        final int     index = Arrays.binarySearch(pos, from);
        if ( index >= 0 )
            return from + offset[index];
        else
            return from + offset[-index - 2];
    }
}
