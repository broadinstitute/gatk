package org.broadinstitute.hellbender.tools.walkers.featuremapping;

import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import org.json.JSONObject;

import java.util.Comparator;
import java.util.LinkedList;
import java.util.List;

public class Interpolator {

    final double[] x;
    final double[] y;
    final PolynomialSplineFunction psf;

    public Interpolator(JSONObject data) {

        // allocate data vectors
        x = new double[data.length()];
        y = new double[x.length];

        // create a sorted list of keys
        List<String> keys = new LinkedList<>(data.keySet());
        keys.sort(new Comparator<String>() {
            @Override
            public int compare(String o1, String o2) {
                return Double.compare(Double.parseDouble(o1), Double.parseDouble(o2));
            }
        });
        int i = 0;
        for ( final String key : keys ) {
            x[i] = Double.parseDouble(key);
            y[i] = data.getDouble(key);
            i++;
        }

        LinearInterpolator li = new LinearInterpolator();
        psf = li.interpolate(x, y);
    }

    public double value(double v) {
        if ( v < x[0] ) {
            return y[0];
        } else if ( v > x[x.length-1] ) {
            return y[y.length-1];
        } else {
            return psf.value(v);
        }
    }
}
