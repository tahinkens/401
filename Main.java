package COMMoDORE;

import java.util.Random;
import java.util.Arrays;

/**
 * Discussion notes:
 * 
 * 
 */
/**
 * Contains main logic.
 * 
 * @author thinkens
 * on 3/21/2018
 * @since 0.1.0
 */
public class Main {
    
    /**
     * Number of dimensions this simulation exists in, unitless.
     * Acceptable values 1, 2, 3.
     */
    static final int DIMENSIONS = 1;
    /**
     * Number of atoms to be simulated, unitless.
     */
    static final int N = 2;
    /**
     * Timestep used during Verlet evolution, s. Default 1 fs.
     */
    static final double TIMESTEP = 1e-15;
    /**
     * The number of timesteps to equilibrate for, unitless.
     */
    static final int NUM_EQUIL_STEPS = 5000; //1500
    /**
     * The number of timesteps to evolve through, unitless.
     */
    static final int NUM_TIMESTEPS = 3000;
    /**
     * Avogadro's number, unitless.
     */
    static final double N_A = 6.022e23;
    /**
     * Boltzmann's constant, J-K^-1.
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
    static final double[] E_X = {1,0,0}, E_Y = {0,1,0}, E_Z = {0,0,1};
    /**
     * The Bravais lattice of the material to be simulated.
     */
    static final CrystalLattice LATTICE = CrystalLattice.FCC;
    /**
     * Uniform random number generator.
     */
    static final Random RNG = new Random(666+1);
    
    /**
     * 
     * 
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        
        //Data
        /**
         * Contains the potential energy of the system at every timestep during
         * equilibration.
         */
        double[] eqPotentialData = new double[NUM_EQUIL_STEPS];
        /**
         * Contains the kinetic energy of the system at every timestep during
         * equilibration.
         */
        double[] eqKineticData = new double[NUM_EQUIL_STEPS];
        /**
         * Contains the total (potential + kinetic) energy of the system at every 
         * timestep during equilibration.
         */
        double[] eqTotalEnergyData = new double[NUM_EQUIL_STEPS];
        /**
         * Contains temperature at every timestep during equilibration.
         */
        double[] eqTemperatureData = new double[NUM_EQUIL_STEPS];
        /**
         * Contains system potential energy at every timestep during simulation.
         */
        double[] potentialData = new double[NUM_TIMESTEPS];
        /**
         * Contains system kinetic energy at every timestep during simulation.
         */
        double[] kineticData = new double[NUM_TIMESTEPS];
        /**
         * Contains system total (U + K) energy at every timestep during simulation.
         */
        double[] totalEnergyData = new double[NUM_TIMESTEPS];
        /**
         * Contains temperature at every timestep.
         */
        double[] temperatureData = new double[NUM_TIMESTEPS];
        
        
        //Program variables
        /**
         * Print excessive amounts of simulation information?.
         */
        boolean verbose = false;
        /**
         * Potential energy of the system, J.
         */
        double systemPotential;
        /**
         * Kinetic energy of the system, J.
         */
        double systemKinetic;
        /**
         * The thermodynamic temperature of the system, K.
         */
        double temperature = 298;
        /**
         * The sum of the kinetic and potential energy of the system, J.
         */
        double totalEnergy;
        /**
         * Position of an atom at time t. Used for Verlet integration.
         */
        double[] currPosition;
        /**
         * Position of an atom at time t-dt. Used for Verlet integration.
         */
        double[] prevPosition;
        /**
         * Position of an atom at time t+dt. Used for Verlet integration.
         */
        double[] newPosition;
        /**
         * Velocity of an atom at time t+dt. Used for obtaining temperature.
         */
        double[] newVelocity;
        
        
        //Begin logic
        Atom[] atoms = new Atom[N];             //instantiate new atoms
        for(int i = 0; i < atoms.length; i++) {
            atoms[i] = new Atom(AtomInfo.ALUMINUM);
        }
        System.out.println("Atoms instantiated: " + atoms.length);
        
        double[][] atomPositions = {{0},{1.23e-10*R_MIN_FACTOR}};
        atoms[0].setPosition(atomPositions[0]); //initialize atom positions
        atoms[1].setPosition(atomPositions[1]);
        System.out.println("\nAtom positions set according to given parameters");
        if(verbose) System.out.println(Arrays.toString(atomPositions));
        
        for(int j = 0; j < N; j++) {                //initialize atom velocities and momenta
            atoms[j].initializeVelocityVector(temperature);
            atoms[j].momentum();
            if(verbose) System.out.println("\nVelocity (m/s): " + Arrays.toString(atoms[j].getVelocity()) + 
                    "\nMomentum (J): " + Arrays.toString(atoms[j].momentum()));
        }
        System.out.println("Atom velocities generated and momenta initialized");
        
        systemPotential = Simulation.potentialEnergy(atoms); //calculate initial system energies
        if(verbose) System.out.println("\nInitial potential energy of system (J): " + systemPotential);
        systemKinetic = Simulation.kineticEnergy(atoms);
        if(verbose) System.out.println("Initial kinetic energy of system (J): " + systemKinetic);
        totalEnergy = systemPotential + systemKinetic;
        System.out.println("Total energy of system (J): " + totalEnergy);

        System.out.println("\nBeginning equilibration run");

        for(int step = 0; step < NUM_EQUIL_STEPS - 1; step++) {
            if(step == 0) { //change positions of atoms for step 0 of verlet
                for(int j = 0; j < N; j++) { 
                    Atom atom = atoms[j];
                    currPosition = atom.getPosition();
                    double[] velocity = atom.getVelocity();
                    double[] accel = atom.acceleration(atoms);
                    newPosition = new double[DIMENSIONS];
                    
                    atom.equilibrationPositionList[step] = currPosition;
                    
                    atom.setPrevPosition(currPosition); //set previous position to their starting position on the lattice
                    for(int i = 0; i < DIMENSIONS; i++) { //use Taylor polynomial to estimate first evolution
                        newPosition[i] = currPosition[i] - (velocity[i] * TIMESTEP) + (0.5d * accel[i] * (TIMESTEP*TIMESTEP));
                    }
                    atom.equilibrationPositionList[step+1] = newPosition;
                    atom.setPosition(newPosition); //assign evolved position to atom
                    atom.setVelocity(atom.verletVelocity(atom.getPosition(),atom.getPrevPosition()));

                    atoms[j] = atom;
                    if(verbose) System.out.println("\nAtoms pre-evolved for initial Verlet integration\nBeginning Verlet evolution");
                }
                
            }
            else {
                for(int j = 0; j < N; j++) { //evolve atoms using Verlet algorithm
                    currPosition = atoms[j].equilibrationPositionList[step];    //r(t)
                    prevPosition = atoms[j].equilibrationPositionList[step - 1];//r(t - dt)

                    newPosition = atoms[j].verlet(currPosition,prevPosition,atoms[j].acceleration(atoms)); //r(t + dt)
                    atoms[j].equilibrationPositionList[step + 1] = newPosition; //r(t + dt)

                    newVelocity = atoms[j].verletVelocity(newPosition,prevPosition);
                    atoms[j].setPosition(newPosition);
                    
                    atoms[j].setVelocity(newVelocity);
                    
                    atoms[j].setPrevPosition(atoms[j].equilibrationPositionList[step]);
                }
            }
            //System.out.println(Arrays.toString(atoms));

            systemKinetic = Simulation.kineticEnergy(atoms);
            systemPotential = Simulation.potentialEnergy(atoms);
            totalEnergy = systemKinetic + systemPotential;
            temperature = Simulation.temperature(systemKinetic,N);
            //System.out.println(systemKinetic / (AtomInfo.ALUMINUM.homopotentialWellDepth * N));
            //System.out.println(systemPotential / (AtomInfo.ALUMINUM.homopotentialWellDepth * N));
            //System.out.println(totalEnergy / (AtomInfo.ALUMINUM.homopotentialWellDepth * N));
            //System.out.println(temperature);
            
            //add energy data/thermo averages to arrays
            eqKineticData[step] = systemKinetic;
            eqPotentialData[step] = systemPotential;
            eqTotalEnergyData[step] = systemKinetic + systemPotential;
        }
        //System.out.println("Equilibration finished, beginning simulation steps");
        
        for(int step = 0; step < NUM_TIMESTEPS; step++) {
            
        }
        
        for(int i = 0; i < atoms[0].equilibrationPositionList.length; i++) {
            //System.out.println(Arrays.toString(atoms[0].equilibrationPositionList[i]));
            //System.out.println(Arrays.toString(atoms[1].equilibrationPositionList[i]));
            System.out.println(atoms[1].equilibrationPositionList[i][0]);
            //System.out.println(atoms[0].equilibrationPositionList[i][0] - atoms[1].equilibrationPositionList[i][0]); //interatomic distance
        }
        //System.out.println(atoms[0].getSpeed());
        //test(); //for testing
    }
    
    private static void test() {
        
        final double testE = 0.997; // kJ/mol (energy/mol)
        final double testS = 3.4e-10, testR = 4e-10; //m (length)
        
        //double potential = lennardJones(testE,testS,testR);
        //System.out.println("Test potential = " + potential + " J");
        
        
        Lattice testLattice = new Lattice(LATTICE,2,AtomInfo.ALUMINUM,true);
        //Atom[] atoms = new Atom[testLattice.getInhabitants().length];
        for(int i = 0; i < testLattice.getInhabitants().length; i++) {
            //testLattice.setInhabitants(atoms);
            //atoms[i] = new Atom(AtomInfo.ALUMINUM);
            //atoms[i].initializeVelocityVector(298);
            //atoms[i].update(TIMESTEP,atoms); //atoms is not the right argument here
            //System.out.println(Arrays.toString(atoms[i].getVelocity()) + ", " + atoms[i].getSpeed());
        }
        //testLattice.setInhabitants(atoms);
        
        //System.out.println(Arrays.toString(testLattice.getInhabitants()));
        
        //double distance = atoms[0].getDistanceFromNeighbor(atoms[1]);
        //System.out.println("Interatomic sep/m : " + distance);
        
        double[] pos = {0,0,0};
        double[] pos2 = {2e-10,0,0};
        
        Atom[] atoms = new Atom[N];
        for(int i = 0; i < atoms.length; i++) {
            atoms[i] = new Atom(AtomInfo.ALUMINUM);
            atoms[i].initializeVelocityVector(298);
        }
        atoms[0].setPosition(pos);
        atoms[1].setPosition(pos2);
        
        for(int i = 0; i < NUM_EQUIL_STEPS; i++) {
            for(int j = 0; j < N; j++) {
                if(i != 0) atoms[j].verlet2(atoms);
                else {
                    double[] newPos = new double[DIMENSIONS];
                    for(int k = 0; k < DIMENSIONS; k++) {
                        newPos[k] = atoms[j].getPosition()[k] - (atoms[j].getVelocity()[k] * TIMESTEP) + (0.5d * atoms[j].acceleration(atoms)[k] * (TIMESTEP*TIMESTEP));
                    }
                    atoms[j].setPrevPosition(pos2);
                    atoms[j].setPosition(newPos);
                }
                System.out.println(Arrays.toString(atoms[j].getPosition()));
            }
        }
        
//        double distance = atoms[0].getDistanceFromNeighbor(atoms[1]);
//        System.out.println("Interatomic sep : " + distance*1e10 + " \u212B"); //UNICODE ANGSTROMS
//        //System.out.println(atoms[0].toString());
//        //System.out.println(atoms[0].getDistanceFromNeighbor(atoms[1]));
//        System.out.println("Force vector between atoms/N : " + Arrays.toString(atoms[0].force(atoms)));
//        System.out.println("Acceleration on atom i/m-s^-2 : " + Arrays.toString(atoms[0].acceleration(atoms)));
//        //System.out.println("Potential energy of system/J : " + Atom.potentialEnergy(atoms));
//        //System.out.println("Total kinetic energy of system (J): " + Atom.kineticEnergy(atoms));
//        System.out.println("Total energy of system (J): ");
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
