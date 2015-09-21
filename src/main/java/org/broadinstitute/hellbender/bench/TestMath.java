package org.broadinstitute.hellbender.bench;

/**
 * Demonstrates that adding small doubles in different orders doesn't always give the same answer.
 */
public class TestMath {

    public static void main(String[] args) {
        System.out.println("For the small numbers we care about, double accumulate numerical error.");
        System.out.println("So summing the same numbers in a different order may yield different answers. See:");
        one();
        two();
        System.out.println("Using a multiplier completely sidesteps this problem, as below demonstrates:");
        oneX();
        twoX();
    }

    public static void one() {
        double x = 0.0;
        x += 1.0;
        x += 0.025;
        x += 0.025;
        x += 0.025;
        x += 0.025;
        x += 0.125;
        System.out.println(String.format("x=%f. It truncates to: %.2f",x,x));
    }

    public static void two() {
        double x = 0.0;
        x += 0.025;
        x += 0.025;
        x += 0.025;
        x += 0.025;
        x += 0.125;
        x += 1.0;
        System.out.println(String.format("x=%f. It truncates to: %.2f", x, x));
    }

    public static void oneX() {
        double x = 0.0;
        x += 1000.0;
        x += 0025.0;
        x += 0025.0;
        x += 0025.0;
        x += 0025.0;
        x += 0125.0;
        x = x/1000.0;
        System.out.println(String.format("x=%f. It truncates to: %.2f",x,x));
    }

    public static void twoX() {
        double x = 0.0;
        x += 0025.0;
        x += 0025.0;
        x += 0025.0;
        x += 0025.0;
        x += 0125.0;
        x += 1000.0;
        x = x/1000.0;
        System.out.println(String.format("x=%f. It truncates to: %.2f",x,x));
    }

}
