package COMMoDORE;

import static COMMoDORE.Main.*;
import java.util.Arrays;
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
    
    private double v, p, pe, ke; //scalar properties
    private double[] velocity, momentum = {0,0,0}; //vector properties
    private double[] position = {0,0,0};
    
    private final double m;
    
    private final Random rng = new Random();
    
    /**
     * Instantiates a new Atom of species atomInfo. 
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
     * Reassigns values to the properties of an atom. Calculates new values for
     * velocity, momentum, their magnitudes, kinetic and potential energies,
     * and the position of the atom.
     * 
     * @param timestep the amount of time passed since last update (this should
     * generally be the universal timestep of the simulation)
     */
    public void update(double timestep) {
        
        v = Math.sqrt(velocity[0]*velocity[0] + velocity[1]*velocity[1] + velocity[2]*velocity[2]);
        p = m * v;
        ke = 0.5 * m * v*v;
        //pe = lennard jones?
        for(int i = 0; i < 3; i++) {
            momentum[i] = m * velocity[i];
        }
        for(int i = 0; i < 3; i++) { //x1 = x0 + vt
            position[i] = position[i] + velocity[i] * timestep;
        }
        //velocity vector needs to change due to particle interactions
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
        
        double x = 0, y = 0, s = 2; //s must be initialized to something > 1
        double[] rv = new double[2];
        
        while(s > 1) {
            x = rng.nextDouble() * 2 - 1; //uniform double b/w -1 and 1
            y = rng.nextDouble() * 2 - 1;
            s = x*x + y*y;
        }
        rv[0] = x * Math.sqrt((-2 * Math.log(s)) / s); //transform to normals
        rv[1] = y * Math.sqrt((-2 * Math.log(s)) / s);
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
        
        return v;
    }
    
    public void setPosition(double[] pos) {
        
        this.position = pos;
    }
    
    @Override
    public String toString() {
        
        return "\\mag(v)=" + v + " p=" + p + " ke=" + ke + " pe=" + pe + 
                " m=" + m + " v=" + Arrays.toString(velocity) + " r=" + 
                Arrays.toString(position);
    }
}
