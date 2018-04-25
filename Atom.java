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
    public double[][] equilibrationPositionList = new double[NUM_EQUIL_STEPS][DIMENSIONS];
    
    private double v, p, pe, ke; //scalar quantities
    private double[] velocity, momentum, position = new double[DIMENSIONS], prevPosition; //vector quantities
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
    protected double[] acceleration(Atom[] atoms) {
        
        double[] accel = new double[DIMENSIONS];
        double[] force = this.force(atoms);
        
        for(int i = 0; i < DIMENSIONS; i++) {
            accel[i] = force[i] / this.atomType.m;
        }
        //System.out.println("Acceleration: " + Arrays.toString(accel)); //debug
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

        double[] f = {0,0,0}; //the force vector
        double dLJ; //the derivative of the lennard jones potential
        double[] r_hat, r1 = this.position, r2, diff = new double[DIMENSIONS]; //to find the relative pos vector between two atoms

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
        }
        //System.out.println("\nForce on atom: " + Arrays.toString(f)); //debug
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
     * distribution by creating randomly sampled normal values scaled by a 
     * factor of sqrt(kT/m) where k is Boltzmann's constant, T is thermodynamic 
     * temp, and m is the mass of the atom of the given species.
     * 
     * @param T Thermodynamic temperature (temperature)
     * @return a randomly generated particle vel
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
     * interatomic separation of r, and equilibrium separation of s. Note: it
     * may be more computationally efficient to calculate the 6 term first and
     * square it to yield the 12 term.
     * 
     * @param e Potential well depth (energy)
     * @param s Equilibrium separation (length)
     * @param r Interatomic separation (length)
     * @return Interatomic potential (energy)
     */
    protected double lennardJones(double e, double s, double r) {

        return 4 * e * (Math.pow((s/r),12) - Math.pow((s/r),6));
    }
    
    /**
     * Calculates the momentum of this atom using its mass and velocity.
     * 
     * @return a vector containing this particle's momentum
     */
    protected double[] momentum() {
        
        momentum = new double[DIMENSIONS];
        
        for(int i = 0; i < DIMENSIONS; i++) {
            momentum[i] = this.atomType.m * velocity[i]; //p = mv
        }
        return momentum;
    }

    /**
     * Implements a basic Verlet algorithm to evolve the position of an atom
     * over time. Since 
     * 
     * @param r_initial the current position of the atom
     * @param r_previous the position of the atom one timestep previously
     * @param accel the acceleration vector acting on the atom
     * @return a new position vector for the atom representing the forces 
     * acted on it
     */
    protected synchronized double[] verlet(double[] r_initial, double[] r_previous, double[] accel) {
        
        double[] r_final = new double[DIMENSIONS];
        
        for(int i = 0; i < DIMENSIONS; i++) {
            r_final[i] = (2d * r_initial[i]) - r_previous[i] + (accel[i] * (TIMESTEP*TIMESTEP));
        }
        return r_final;
    }
    
    protected synchronized double[] verlet2(Atom[] atoms) {
        
        //double dt = TIMESTEP;
        //double[] r_t = positionList[step], r_tminusdt = positionList[step - 1], r_tplusdt = new double[DIMENSIONS];
        double[] r_t = this.position, r_tminusdt = this.prevPosition, r_tplusdt = new double[DIMENSIONS];
        double[] accel = this.acceleration(atoms);
        
        for(int i = 0; i < DIMENSIONS; i++) {
            r_tplusdt[i] = (2d * r_t[i]) - r_tminusdt[i] + (accel[i] * (TIMESTEP*TIMESTEP));
        }
        this.prevPosition = r_t;
        this.position = r_tplusdt;
        
        return r_tplusdt;
    }
    
    /**
     * Approximates the velocity of an atom after Verlet evolution. Note that
     * velocity is not a function of current position, rather it is a function
     * of the position at the previous and next times.
     * 
     * @param r_final the position of the atom at time t+dt
     * @param r_previous the position of the atom at time t-dt
     * @return the approximate velocity after evolution
     */
    protected synchronized double[] verletVelocity(double[] r_final, double[] r_previous) {
        
        double[] vel = new double[DIMENSIONS];
        
        for(int i = 0; i < DIMENSIONS; i++) {
            vel[i] = (r_final[i] - r_previous[i]) / (2f * TIMESTEP);
        }
        this.velocity = vel;
        //System.out.println("verletVel: " + Arrays.toString(vel));
        return vel;
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
        
        return MathUtil.magnitude(velocity);
    }
    
    /**
     * 
     * @return the position of this atom one timestep previously
     */
    public double[] getPrevPosition() {
        
        return prevPosition;
    }

    /**
     * 
     * @param prevPosition 
     */
    public void setPrevPosition(double[] prevPosition) {
        
        this.prevPosition = prevPosition;
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
     * @param velocity 
     */
    public void setVelocity(double[] velocity) {
        
        this.velocity = velocity;
    }
    /**
     * 
     * @return a string containing the state of all of this object's fields
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
     * given in m/s.
     * 
     * @deprecated This function isn't used and instead initial particle 
     * velocities are determined using a Marsaglia random generator and the 
     * initializeVelocityVector() function above. It's being left in for
     * possible future use and posterity.
     * 
     * @param T Thermodynamic temperature (temperature)
     * @param atomType The type of atom being simulated (AtomInfo)
     * @return a randomly generated particle speed in m/s
     * 
     */
    @Deprecated
    private double maxwellBoltzmann(double T, AtomInfo atomType) {
       
        final double a = Math.sqrt((K_B * T) / atomType.m); //where K_B is Boltzmann'equilibriumSeparation constant and m is particle mass
        final double x = RNG.nextDouble();
        
        return Math.sqrt(2/PI) * ((x*x * Math.exp(-(x*x) / (2 * a*a))) / (a*a*a));
    }
    
    /**
     * Reassigns values to the properties of an atom. Calculates new values for
 vel, momentum, their magnitudes, kinetic and potential energies,
 and the position of the atom.
     * 
     * @deprecated This method has been broken into smaller pieces and is
     * currently non-functional.
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
