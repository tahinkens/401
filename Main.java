package COMMoDORE;

import java.util.Random;
import java.util.Arrays;

/**
 * Contains main logic
 * @author thinkens
 * on 3/21/2018
 * @since 0.1.0
 */
public class Main {
    
    static final double TIMESTEP = 5e-15; //universal simulation timestep/s
    static final double K_B = 1.38064852e-23; // Boltzmann's constant/J-K^-1
    static final double PI = 3.1415926; //mathematical constant pi
    static final double[] E_X = {1,0,0}, E_Y = {0,1,0} ,E_Z = {0,0,1}; //unit vectors/m,m,m
    static final CrystalLattice FCC = CrystalLattice.FCC;
    
    static Random rng = new Random();
    
    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        
        final double testE = 0.997; // kJ/mol (energy/mol)
        final double testS = 3.4e-10, testR = 4e-10; //m (length)
        
        final double alE = 0.368; //kJ/mol
        final double alInteratomicR = 1.23e-10; //m (eq. sep from Totten and Mackenzie)
        
        //final Atom[] = new
        
        //double potential = lennardJones(testE,testS,testR);
        //System.out.println("Test potential = " + potential + " J");
        
        
        Lattice testLattice = new Lattice(FCC,2,AtomInfo.ALUMINUM,true);
        Atom[] atoms = new Atom[testLattice.getInhabitants().length];
        for(int i = 0; i < testLattice.getInhabitants().length; i++) {
            //testLattice.setInhabitants(atoms);
            atoms[i] = new Atom(AtomInfo.ALUMINUM);
            atoms[i].initializeVelocityVector(3,298);
            //atoms[i].update(TIMESTEP,atoms); //atoms is not the right argument here
            System.out.println(Arrays.toString(atoms[i].getVelocity()) + ", " + atoms[i].getSpeed());
        }
        testLattice.setInhabitants(atoms);
        
        System.out.println(Arrays.toString(testLattice.getInhabitants()));
        
        double distance = testLattice.getInhabitants()[0].getDistanceFromNeighbor(testLattice.getInhabitants()[1]);
        System.out.println(distance);
        
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
