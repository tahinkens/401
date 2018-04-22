package COMMoDORE;

import static COMMoDORE.Main.*;
import java.util.Arrays;

/**
 * Defines properties and functions of atoms.
 * 
 * @author thinkens
 * on 4/4/2018
 * @since 0.1.0
 */
public class Atom {
    
    public final AtomInfo atomType;
    
    private double v, p, pe, ke; //scalar quantities
    private double[] velocity, momentum = {0,0,0}, position = {0,0,0}; //vector quantities
    
    private final double m;

    
    /**
     * Instantiates a new Atom of species atomInfo. 
     * 
     * @param atomInfo The species of atom to be instantiated
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
     * Calculates the acceleration on an atom given all of the other atoms
     * in the system.
     * 
     * @param atoms all of the atoms in the simulation
     * @return a vector representing the acceleration on this atom
     */
    protected synchronized double[] acceleration(Atom[] atoms) {
        
        double[] accel = new double[DIMENSIONS];
        double[] force = this.force(atoms);
        
        for(int i = 0; i < accel.length; i++) {
            accel[i] = (1/this.m) * force[i];
        }
        return accel;
    }
    
    /**
     * Calculates the first derivative of the potential between two atoms. That
     * is, the force of the atom this method is called on and the argument atom.
     * 
     * @param e Potential well depth (energy)
     * @param s Equilibrium separation (length)
     * @param r Interatomic separation (length)
     * @return a double representing the force between two atoms
     */
    private double diffLennardJones(double e, double s, double r) {
        
        return ((-48 * e * Math.pow(s,12)) * Math.pow(r,-13)) + ((24 * e * Math.pow(s,6)) * Math.pow(r,-7));
    }
    
    /**
     * Calculates the force on an atom using the derivative of the Lennard-Jones
     * potential on that atom. Negative forces are attractive.
     * 
     * @param atoms all of the atoms in the simulation
     * @return a double[] representing the force vector on this atom
     */
    protected synchronized double[] force(Atom[] atoms) {

        double[] f = {0,0,0}, lastf = {0,0,0}; //the force vector
        
        double dLJ; //the derivative of the lennard jones potential
        double[] r_hat;//position vector between two atoms
        double[] r1 = this.position, r2, diff = new double[DIMENSIONS]; //to find the relative pos vector between two atoms

        for(Atom atom : atoms) {
            if(!this.equals(atom)) { //perform calculations iff the atoms aren't the same
                r2 = atom.position;
                for(int i = 0; i < DIMENSIONS; i++) { //create distance vector between this atom and the neighbor atom in the array
                    diff[i] = r2[i] - r1[i];
                }
                r_hat = MathUtil.unitVector(diff);
                
                dLJ = diffLennardJones(this.atomType.homopotentialWellDepth,this.atomType.equilibriumSeparation,MathUtil.magnitude(diff));
                for(int i = 0; i < DIMENSIONS; i++) { //create force vector by multiplying force with unit vector
                    f[i] += dLJ * r_hat[i];
                }
            }
//            for(int i = 0; i < DIMENSIONS; i++) { //do we need this? might be necessary for >2 atom systems
//                f[i] += lastf[i];
//            }
        }
        return f;
    }
    
    /**
     * Returns the distance of one particle from another.
     * 
     * @param neighbor the neighboring Atom
     * @return the distance between this Atom and the neighboring Atom
     */
    public double getDistanceFromNeighbor(Atom neighbor) {
        
        //double d; //distance
        double[] distance = new double[DIMENSIONS]; //difference vector between positions
        
        for(int i = 0; i < DIMENSIONS; i++) {
            distance[i] = this.position[i] - neighbor.position[i];
        }
        //d = MathUtil.magnitude(distance);
        return MathUtil.magnitude(distance);
    }
    
    /**
     * Randomly generates an n-dimensional velocity using a scaled Marsaglia 
     * polar method. This method attempts to approximate a Maxwell-Boltzmann
 distribution by creating randomly sampled normal values scaled by a 
 factor of sqrt(kT/m) where k is Boltzmann'equilibriumSeparation constant, T is thermodynamic 
 temp, and m is the mass of the atom of the given species.
     * 
     * @param T Thermodynamic temperature (temperature)
     * @return a randomly generated particle velocity
     */
    public double[] initializeVelocityVector(double T) {
        
        velocity = new double[DIMENSIONS];
        final double a = Math.sqrt((K_B * T) / this.atomType.m);
        
        for(int i = 0; i < DIMENSIONS; i++) {
            velocity[i] = a * MathUtil.marsagliaRandomGenerator();
        }
        return velocity;
    }

    /**
     * Calculates the potential between two atoms with a potential well depth e,
 interatomic separation of equilibriumSeparation, and equilibrium separation of r. Note: it
     * may be more computationally efficient to calculate the 6 term first and
     * square it to yield the 12 term.
     * 
     * @param e Potential well depth (energy)
     * @param s Equilibrium separation (length)
     * @param r Interatomic separation (length)
     * @return Interatomic potential (energy)
     */
    private double lennardJones(double e, double s, double r) {

        return 4 * e * (Math.pow((s/r),12) - Math.pow((s/r),6));
    }
    
    /**
     * Calculates the potential energy of an atomic system as a function of the
     * positions of the particles. More precisely, the interatomic potential of
     * each atom on every other atom is calculated using a standard 
     * Lennard-Jones potential. This energy can then be differentiated with 
     * respect to a particle'equilibriumSeparation position to find the force on that particle.
     * 
     * @param atoms an array containing all of the atoms in the simulation
     * @return the potential energy 
     */
    protected static double potentialEnergy(Atom[] atoms) {
        
        double potentialEnergy = 0;
        AtomInfo species = atoms[0].atomType; //assume that every Atom in atoms is same species
        
        for(Atom atomi : atoms) {
            for(Atom atomj : atoms) {
                if(!atomj.equals(atomi)) //only sum contributions from different atoms
                    potentialEnergy += atomj.lennardJones(species.homopotentialWellDepth,species.equilibriumSeparation,atomj.getDistanceFromNeighbor(atomi));
            }
        }
        return 0.5 * potentialEnergy;
    }

    /**
     * 
     * @return the position vector of this Atom
     */
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
    
    /**
     * 
     * @param pos the desired position to move this atom to
     */
    public void setPosition(double[] pos) {
        
        this.position = pos;
    }
    
    /**
     * 
     * @return a string containing the state of all of this object'equilibriumSeparation fields
     */
    @Override
    public String toString() {
        
        return "\\mag(v)=" + v + " p=" + p + " ke=" + ke + " pe=" + pe + 
                " m=" + m + " v=" + Arrays.toString(velocity) + " r=" + 
                Arrays.toString(position);
    }
    
    /**
     * Generates a random particle speed based on the Maxwell-Boltzmann 
     * distribution based on thermodynamic temperature and particle mass. Speed
 given in m/equilibriumSeparation.
     * 
     * @param T Thermodynamic temperature (temperature)
     * @param atomType The type of atom being simulated (AtomInfo)
     * @return a randomly generated particle speed in m/equilibriumSeparation
 
 Deprecated: This function isn't used and instead initial particle 
 velocities are determined using a Marsaglia random generator and the 
 initializeVelocityVector() function above. It'equilibriumSeparation being left in for
 possible future use and posterity.
     */
    @Deprecated
    private double maxwellBoltzmann(double T, AtomInfo atomType) {
       
        final double a = Math.sqrt((K_B * T) / atomType.m); //where K_B is Boltzmann'equilibriumSeparation constant and m is particle mass
        final double x = RNG.nextDouble();
        
        return Math.sqrt(2/PI) * ((x*x * Math.exp(-(x*x) / (2 * a*a))) / (a*a*a));
    }
    
    /**
     * Reassigns values to the properties of an atom. Calculates new values for
     * velocity, momentum, their magnitudes, kinetic and potential energies,
     * and the position of the atom.
     * 
     * Deprecated: this method is going to be broken up into several smaller
     * functions to increase modularity and readability. Flagged for removal.
     * 
     * @param timestep the amount of time passed since last update (this should
     * generally be the universal timestep of the simulation)
     * @param neighbors an array of Atoms containing the atoms nearby to this
     */
    @Deprecated
    public void update(double timestep, Atom[] neighbors) {

        //update potential energy and differentiate to obtain force
        for(int i = 0; i < neighbors.length; i++) {
            if(this.getDistanceFromNeighbor(neighbors[i]) < 5e-10) //only calculate potential from atoms that are <5 Ang away
            pe = lennardJones(this.atomType.homopotentialWellDepth,this.getDistanceFromNeighbor(neighbors[i]),this.atomType.equilibriumSeparation);
            //does this only work for one atom at a time? Might need to be a sum of all pe'equilibriumSeparation
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
}
