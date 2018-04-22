package COMMoDORE;

import java.util.Random;
import java.util.Arrays;

/**
 * Discussion notes:
 * 
 * 
 */
/**
 * Contains main logic. Currently just testing, but will eventually contain main
 * equilibration loop and possibly rescaling code.
 * 
 * @author thinkens
 * on 3/21/2018
 * @since 0.1.0
 */
public class Main {
    
    /**
     * Number of dimensions this simulation exists in, unitless.
     * Acceptable values 2, 3.
     */
    static final int DIMENSIONS = 3;
    /**
     * Number of atoms to be simulated, unitless.
     */
    static final int N = 2;
    /**
     * Timestep used during Verlet evolution, equilibriumSeparation. Default 5 fs.
     */
    static final double TIMESTEP = 5e-15;
    /**
     * The number of timesteps to evolve through.
     */
    static final int NUM_TIMESTEPS = 10000;
    /**
     * Avogadro'equilibriumSeparation number, unitless.
     */
    static final double N_A = 6.022e23;
    /**
     * Boltzmann'equilibriumSeparation constant, J-K^-1.
     */
    static final double K_B = 1.38064852e-23;
    /**
     * The mathematical constant pi, unitless.
     */
    static final double PI = 3.1415926;
    /**
     * This factor can be used to multiply an atom's equilibrium
     * separation to obtain a force of 0.
     */
    static final double R_MIN_FACTOR = 1.12246204830937298;
    /**
     * Cartesian unit vectors, m.
     */
    static final double[] E_X = {1,0,0}, E_Y = {0,1,0} ,E_Z = {0,0,1};
    /**
     * The Bravais lattice of the material to be simulated.
     */
    static final CrystalLattice LATTICE = CrystalLattice.FCC;
    /**
     * Uniform random number generator.
     */
    static final Random RNG = new Random();
    
    /**
     * 
     * 
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        
        boolean verbose = true;
        
        double systemPotential;
        
        Atom[] atoms = new Atom[2];             //instantiate new atoms
        for(int i = 0; i < atoms.length; i++) {
            atoms[i] = new Atom(AtomInfo.ALUMINUM);
        }
        System.out.println("Atoms instantiated: " + atoms.length);
        
        double[][] atomPositions = {{0,0,0},{4.04e-10,0,0}};
        atoms[0].setPosition(atomPositions[0]); //initialize atom positions
        atoms[1].setPosition(atomPositions[1]);
        System.out.println("\nAtom positions set according to given parameters");
        if(verbose) System.out.println(Arrays.toString(atomPositions));
        
        for(Atom atom : atoms) {                //initialize atom velocities and momenta
            atom.initializeVelocityVector(298);
            atom.momentum();
            if(verbose) System.out.println("\nVelocity (m/s): " + Arrays.toString(atom.getVelocity()) + "\nMomentum (J): " + Arrays.toString(atom.momentum()));
        }
        System.out.println("Atom velocities generated and momenta initialized");
        
        systemPotential = Atom.potentialEnergy(atoms);
        if(verbose) System.out.println("\nPotential energy of system (J): " + systemPotential);
        
        //on step b of LeSar
        
        for(int step = 0; step < NUM_TIMESTEPS; step++) {
            
        }
        
        //test(); //for testing
    }
    
    private static void test() {
        
        final double testE = 0.997; // kJ/mol (energy/mol)
        final double testS = 3.4e-10, testR = 4e-10; //m (length)
        
        //double potential = lennardJones(testE,testS,testR);
        //System.out.println("Test potential = " + potential + " J");
        
        
        Lattice testLattice = new Lattice(LATTICE,2,AtomInfo.ALUMINUM,true);
        Atom[] atoms = new Atom[testLattice.getInhabitants().length];
        for(int i = 0; i < testLattice.getInhabitants().length; i++) {
            //testLattice.setInhabitants(atoms);
            atoms[i] = new Atom(AtomInfo.ALUMINUM);
            atoms[i].initializeVelocityVector(298);
            //atoms[i].update(TIMESTEP,atoms); //atoms is not the right argument here
            System.out.println(Arrays.toString(atoms[i].getVelocity()) + ", " + atoms[i].getSpeed());
        }
        testLattice.setInhabitants(atoms);
        
        //System.out.println(Arrays.toString(testLattice.getInhabitants()));
        
        //double distance = atoms[0].getDistanceFromNeighbor(atoms[1]);
        //System.out.println("Interatomic sep/m : " + distance);
        
        double[] pos = {0,0,0};
        double[] pos2 = {2e-10,0,0};
        atoms[0].setPosition(pos);
        atoms[1].setPosition(pos2);
        double distance = atoms[0].getDistanceFromNeighbor(atoms[1]);
        System.out.println("Interatomic sep : " + distance*1e10 + " \u212B"); //UNICODE ANGSTROMS
        //System.out.println(atoms[0].toString());
        //System.out.println(atoms[0].getDistanceFromNeighbor(atoms[1]));
        System.out.println("Force vector between atoms/N : " + Arrays.toString(atoms[0].force(atoms)));
        System.out.println("Acceleration on atom i/m-s^-2 : " + Arrays.toString(atoms[0].acceleration(atoms)));
        System.out.println("Potential energy of system/J : " + Atom.potentialEnergy(atoms));
        //unit vector testing
        
        //double[] vector = {2,0,-3};
        //System.out.println(Arrays.toString(unitVector(vector)));
        
        //testing of velocity initialization for atoms
//        Atom al = new Atom(AtomInfo.ALUMINUM);
//        al.initializeVelocityVector(3,298,al.atomType);
//        
//        System.out.println(Arrays.toString(al.getVelocity()));
//        System.out.println(al.getSpeed());

        //Old velocity distribution testing
//        double[] v;
//        double magV = 0;
//        for(int i = 0; i < 20; i++) {
//            v = initializeVelocityVector(3,300,XENON);
//            System.out.print(Arrays.toString(v));
//            magV = Math.sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
//            System.out.println("   V = " + magV + " m/s");
//        }
        
        //Testing random particle speed generation using Maxwell Boltzmann
//        double maxSpeed = 0, minSpeed = 1e25;
//        for(int i = 0; i < 30; i++) {
//            double speed = maxwellBoltzmann(298,COPPER);
//            if(speed > maxSpeed) maxSpeed = speed;
//            if(speed < minSpeed) minSpeed = speed;
//            System.out.println("Randomly sampled particle speed = " + speed + " m/s");
//        }
//        System.out.println("Maximum particle speed = " + maxSpeed + " m/s");
//        System.out.println("Minimum particle speed = " + minSpeed + " m/s");
    }
}
