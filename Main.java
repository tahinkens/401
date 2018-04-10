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
    
    static Random rng = new Random();
    
    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        
        double e = 0.997; // kJ/mol (energy)
        double s = 3.4e-10, r = 4e-10; //m (length)
        
        double potential = lennardJones(e,s,r);
        System.out.println("Interatomic potential = " + potential + " J");
        
        Lattice testLattice = new Lattice(FCC,2);
        
        System.out.println(Arrays.toString(testLattice.getInhabitants()));
        double[] v;
        double magV = 0;
        
        for(int i = 0; i < 20; i++) {
            v = initializeVelocityVector(3,300,HELIUM);
            System.out.print(Arrays.toString(v));
            magV = Math.sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
            System.out.println("   V = " + magV + " m/s");
        }
        
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
    
    /**
     * Randomly generates an n-dimensional velocity using a scaled Marsaglia 
     * polar method. This method attempts to approximate a Maxwell-Boltzmann
     * distribution by creating randomly sampled normal values scaled by a factor
     * of sqrt(kT/m) where k is Boltzmann's constant, T is thermodynamic temp,
     * and m is the mass of the atom of a given species.
     * 
     * @param dimensions How many dimensions does this particle exist in?
     * @param T Thermodynamic temperature (temperature)
     * @param atomType The type of atom being simulated (AtomInfo)
     * @return a randomly generated particle velocity
     */
    private static double[] initializeVelocityVector(int dimensions, double T, AtomInfo atomType) {
        
        double[] velocity = new double[dimensions];
        final double a = Math.sqrt((K_B * T) / atomType.m);
        
        for(int i = 0; i < dimensions; i++) {
            velocity[i] = a * marsagliaRandomGenerator()[0];
        }        
        return velocity;
    }
    
    /**
     * Generates a pair of normally distributed random variables with mean 0 and
     * variance 1.
     * 
     * @return an array containing two randomly generated numbers
     */
    private static double[] marsagliaRandomGenerator() {
        
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
    
}
