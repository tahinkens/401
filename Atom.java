package COMMoDORE;

import static COMMoDORE.Main.*;
import java.util.Random;

/**
 * Defines properties and functions of atoms.
 * 
 * @author thinkens
 * on 4/4/2018
 * @since 0.1.0
 */
public class Atom {
    
    public AtomInfo atomType;
    
    private double v, p, pe, ke, m;
    private double[] velocity;
    private double[] position = {0,0};
    
    private final Random rng = new Random();
    
    /**
     * Instantiates a new Atom of type
     * @param atomInfo 
     */
    public Atom(AtomInfo atomInfo) {
        
        this.atomType = atomInfo;
        
        this.v = 0;
        this.p = 0;
        this.ke = 0;
        this.pe = 0;
        this.m = atomInfo.m;
    }
    
    /**
     * Randomly generates an n-dimensional velocity using a scaled Marsaglia 
     * polar method. This method attempts to approximate a Maxwell-Boltzmann
     * distribution by creating randomly sampled normal values scaled by a 
     * factor of sqrt(kT/m) where k is Boltzmann's constant, T is thermodynamic 
     * temp, and m is the mass of the atom of a given species.
     * 
     * @param dimensions How many dimensions does this particle exist in?
     * @param T Thermodynamic temperature (temperature)
     * @param atomType The type of atom being simulated (AtomInfo)
     * @return a randomly generated particle velocity
     */
    public double[] initializeVelocityVector(int dimensions, double T, AtomInfo atomType) {
        
        velocity = new double[dimensions];
        final double a = Math.sqrt((K_B * T) / atomType.m);
        
        for(int i = 0; i < dimensions; i++) {
            velocity[i] = a * marsagliaRandomGenerator()[0];
        }
        return velocity;
    }
    
    /**
     * Generates a pair of normally distributed random variables with mean 0 and
     * variance 1. Used for assigning randomly generated velocities to 
     * particles.
     * 
     * @return an array containing two randomly generated numbers
     */
    private double[] marsagliaRandomGenerator() {
        
        double u = 0, v = 0, s = 2;
        double[] rv = new double[2];
        
        while(s > 1) {
            u = rng.nextDouble() * 2 - 1;
            v = rng.nextDouble() * 2 - 1;
            s = u*u + v*v;
        }
        rv[0] = u * Math.sqrt((-2 * Math.log(s)) / s);
        rv[1] = v * Math.sqrt((-2 * Math.log(s)) / s);
        return rv;
    }
    
    /**
     * 
     * @return an array containing each component of the velocity
     */
    public double[] getVelocity() {

        return velocity;
    }
    
    /**
     * 
     * @return the magnitude of the velocity vector
     */
    public double getSpeed() {
        
        return Math.sqrt(velocity[0]*velocity[0] + velocity[1]*velocity[1] + velocity[2]*velocity[2]);
    }
    
    public void setPosition(double[] pos) {
        
        this.position = pos;
    }
}
