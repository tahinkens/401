package COMMoDORE;

import static COMMoDORE.Main.*;

/**
 * MathUtil is a class of static methods used for basic mathematical functions
 * such as finding unit vectors and magnitudes.
 *
 * @author thinkens
 * on 4/19/2018
 * @since 0.1.0
 */
public class MathUtil {
    
    /**
     * Calculates the magnitude of a given vector.
     * 
     * @param vector the vector of interest
     * @return a double containing the length of the vector
     */
    public static double magnitude(double[] vector) {
        
        return Math.sqrt(vector[0]*vector[0] + vector[1]*vector[1] + vector[2]*vector[2]);
    }
    
    /**
     * Generates a normally distributed random variable with mean 0 and
     * variance 1. Used for assigning randomly generated velocities to 
     * particles. The Marsaglia method generates two paired values by default,
     * but since only one is needed only one is returned.
     * 
     * @return a random, normally distributed double between -1 and 1
     */
    public static double marsagliaRandomGenerator() {
        
        double x = 0, y, s = 2; //s must be initialized to something > 1
        double rv;
        
        while(s > 1) {
            x = RNG.nextDouble() * 2 - 1; //uniform double b/w -1 and 1
            y = RNG.nextDouble() * 2 - 1;
            s = x*x + y*y;
        }
        rv = x * Math.sqrt((-2 * Math.log(s)) / s);
        return rv;
    }
    
    /**
     * Takes a vector and divides it by its magnitude.
     * 
     * @param vector the vector to normalize
     * @return the corresponding unit vector to the vector argument
     */
    public static double[] unitVector(double[] vector) {
        
        double[] unit = new double[vector.length];
        double magnitude = Math.sqrt(vector[0]*vector[0] + vector[1]*vector[1] + vector[2]*vector[2]);
        
        for(int i = 0; i < vector.length; i++) {
            unit[i] = vector[i] / magnitude;
        }
        return unit;
    }
    
    /**
     * Creates a zero vector of a given length.
     * 
     * @param length how long the vector should be
     * @return a vector full of 0's
     */
    public static double[] zeroes(int length) {
        
        double[] zero = new double[length];
        for(int i = 0; i < length; i++) {
            zero[i] = 0;
        }
        return zero;
    }
}
