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

    static final double K_B = 1.38064852e-23; // J/K
    static final double PI = 3.1415926;
    static final CrystalLattice FCC = CrystalLattice.FCC;
    static final AtomInfo ALUMINUM = AtomInfo.ALUMINUM;
    static final AtomInfo HELIUM = AtomInfo.HELIUM;
    static final AtomInfo XENON = AtomInfo.XENON;
    
    static Random rng = new Random();
    
    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        
        double e = 0.997; // kJ/mol (energy)
        double s = 3.4e-10, r = 4e-10; //m (length)
        
        double potential = lennardJones(e,s,r);
        System.out.println("Interatomic potential = " + potential + " J");
        
        Atom[] atoms = {};
        Lattice testLattice = new Lattice(FCC,2,AtomInfo.ALUMINUM);
        testLattice.setInhabitants(atoms);
        
        System.out.println(Arrays.toString(testLattice.getInhabitants()));
        
        Atom al = new Atom(ALUMINUM);
        al.initializeVelocityVector(3, 298, al.atomType);
        
        System.out.println(Arrays.toString(al.getVelocity()));
        System.out.println(al.getSpeed());

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
    
    /**
     * Calculates the potential between two atoms with a potential well depth e,
     * interatomic separation of s, and equilibrium separation of r. Interatomic
     * potential given in Joules.
     * 
     * @param e Potential well depth (energy)
     * @param s Interatomic separation (length)
     * @param r Equilibrium separation (length)
     * @return Interatomic potential (J)
     */
    private static double lennardJones(double e, double s, double r) {

        return 4 * e * (Math.pow((s/r),12) - Math.pow((s/r),6));
    }
    
    /**
     * Generates a random particle speed based on the Maxwell-Boltzmann 
     * distribution based on thermodynamic temperature and particle mass. Speed
     * given in m/s.
     * 
     * @param T Thermodynamic temperature (temperature)
     * @param atomType The type of atom being simulated (AtomInfo)
     * @return a randomly generated particle speed in m/s
     */
    private static double maxwellBoltzmann(double T, AtomInfo atomType) {
       
        final double a = Math.sqrt((K_B * T) / atomType.m); //where K_B is Boltzmann's constant and m is particle mass
        final double x = rng.nextDouble();
        
        return Math.sqrt(2/PI) * ((x*x * Math.exp(-(x*x) / (2 * a*a))) / (a*a*a)); //FIXME this is the y-value of the pdf not speed
    }
    
    
    
}
