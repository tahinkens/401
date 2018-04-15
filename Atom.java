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
     * @param neighbors an array of Atoms containing the atoms nearby to this
     */
    public void update(double timestep, Atom[] neighbors) {

        //update potential energy and differentiate to obtain force
        for(int i = 0; i < neighbors.length; i++) {
            if(this.getDistanceFromNeighbor(neighbors[i]) < 5e-10) //only calculate potential from atoms that are <5 Ang away
            pe = lennardJones(this.atomType.homopotentialWellDepth,this.getDistanceFromNeighbor(neighbors[i]),this.atomType.epsilon);
            //does this only work for one atom at a time? Might need to be a sum of all pe's
        }
        
        //velocity vector needs to change due to particle interactions
        v = Math.sqrt(velocity[0]*velocity[0] + velocity[1]*velocity[1] + velocity[2]*velocity[2]);
        ke = 0.5 * m * v*v;
        
        for(int i = 0; i < 3; i++) {
            momentum[i] = m * velocity[i];
        }
        p = m * v;
        for(int i = 0; i < 3; i++) { //x1 = x0 + vt
            position[i] = position[i] + velocity[i] * timestep;
        }
    }
    
    /**
     * Randomly generates an n-dimensional velocity using a scaled Marsaglia 
     * polar method. This method attempts to approximate a Maxwell-Boltzmann
     * distribution by creating randomly sampled normal values scaled by a 
     * factor of sqrt(kT/m) where k is Boltzmann's constant, T is thermodynamic 
     * temp, and m is the mass of the atom of a given species.
     * 
     * @param dimensions How many dimensions does this particle exist in? For
     * the extent of this project, this should be fixed at 3.
     * @param T Thermodynamic temperature (temperature)
     * @param atomType The type of atom being simulated (AtomInfo)
     * @return a randomly generated particle velocity
     */
    public double[] initializeVelocityVector(int dimensions, double T, AtomInfo atomType) {
        
        velocity = new double[dimensions];
        final double a = Math.sqrt((K_B * T) / atomType.m);
        
        for(int i = 0; i < dimensions; i++) {
            velocity[i] = a * marsagliaRandomGenerator();
        }
        return velocity;
    }
    
    /**
     * Calculates the potential between two atoms with a potential well depth e,
     * interatomic separation of s, and equilibrium separation of r.
     * 
     * @param e Potential well depth (energy)
     * @param s Interatomic separation (length)
     * @param r Equilibrium separation (length)
     * @return Interatomic potential (energy)
     */
    private double lennardJones(double e, double s, double r) {

        return 4 * e * (Math.pow((s/r),12) - Math.pow((s/r),6));
    }
    
    /**
     * Generates a pair of normally distributed random variables with mean 0 and
     * variance 1. Used for assigning randomly generated velocities to 
     * particles. The Marsaglia method generates two paired values by default,
     * but since only one is needed only one is returned.
     * 
     * @return a random, normally distributed double between -1 and 1
     */
    private double marsagliaRandomGenerator() {
        
        double x = 0, y, s = 2; //s must be initialized to something > 1
        //double[] rv = new double[2];
        double rv;
        
        while(s > 1) {
            x = rng.nextDouble() * 2 - 1; //uniform double b/w -1 and 1
            y = rng.nextDouble() * 2 - 1;
            s = x*x + y*y;
        }
        //rv[0] = x * Math.sqrt((-2 * Math.log(s)) / s); //transform to normals
        //rv[1] = y * Math.sqrt((-2 * Math.log(s)) / s);
        rv = x * Math.sqrt((-2 * Math.log(s)) / s);
        return rv;
    }
    
    /**
     * Generates a random particle speed based on the Maxwell-Boltzmann 
     * distribution based on thermodynamic temperature and particle mass. Speed
     * given in m/s.
     * 
     * @param T Thermodynamic temperature (temperature)
     * @param atomType The type of atom being simulated (AtomInfo)
     * @return a randomly generated particle speed in m/s
     * 
     * @deprecated This function isn't used and instead initial particle 
     * velocities are determined using a Marsaglia random generator and the 
     * initializeVelocityVector() function above. It's being left in for
     * possible future use and posterity.
     */
    @Deprecated
    private double maxwellBoltzmann(double T, AtomInfo atomType) {
       
        final double a = Math.sqrt((K_B * T) / atomType.m); //where K_B is Boltzmann's constant and m is particle mass
        final double x = rng.nextDouble();
        
        return Math.sqrt(2/PI) * ((x*x * Math.exp(-(x*x) / (2 * a*a))) / (a*a*a)); //FIXME this is the y-value of the pdf not speed
    }
    
    /**
     * 
     * @param neighbor the neighboring Atom
     * @return the distance between this Atom and the neighboring Atom
     */
    public double getDistanceFromNeighbor(Atom neighbor) {
        
        double d; //distance
        double[] distance = new double[3]; //difference vector between positions
        
        for(int i = 0; i < 3; i++) { //why 3?
            distance[i] = this.position[i] - neighbor.position[i];
        }
        d = Math.sqrt(distance[0]*distance[0] + distance[1]*distance[1] + distance[2]*distance[2]);
        return d;
    }
    
    public double[] getPosition() {
        
        return position;
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
