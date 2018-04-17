package COMMoDORE;

/**
 *
 * @author thinkens
 * on 4/4/2018
 * @since 0.1.0
 */
public class Lattice {
    
    private Atom[] inhabitants;
    
    private double a = 0, b = 0, c = 0; //lattice parameters
    
    
    /**
     * 
     * @param latticeType desired Bravais lattice (SimpleCubic, FCC, BCC, HCP)
     * @param N desired number of atoms to simulate
     * @param atomType the species to be simulated
     * @param autoAssignLatticeParameters obtain lattice parameters from AtomInfo?
     */
    public Lattice(CrystalLattice latticeType, int N, AtomInfo atomType, boolean autoAssignLatticeParameters) {
        
        inhabitants = new Atom[N];
        inhabitants = initializePositions(latticeType);
        if(autoAssignLatticeParameters) {
            a = atomType.a;
            b = atomType.b;
            c = atomType.c;
        }
        
    }
    
    /**
     * 
     * @return an array of atoms with all of their position set to the correct
     * spots on the desired crystal lattice
     */
    private Atom[] initializePositions(CrystalLattice latticeType) {
        
        //double latticeWidth;
        //double latticeHeight;
        
        for(Atom atom : inhabitants) {
            switch(latticeType) {
                case SimpleCubic:
                    
                    break;
                case FCC:
                    
                    break;
                case BCC:
                    
                    break;
                case HCP:
                    
                    break;
            }
        }
        return inhabitants;
    }
    
    /**
     * 
     * @return an array containing the lattice parameters a, b, and c in indices
     * 0, 1, and 2.
     */
    public double[] getDimensions() {
        
        double[] dimensions = new double[3];
        dimensions[0] = a;
        dimensions[1] = b;
        dimensions[2] = c;
        return dimensions;
    }

    public Atom[] getInhabitants() {
        
        return inhabitants;
    }
    
    public void setInhabitants(Atom[] inhabitants) {
        
        this.inhabitants = inhabitants;
    }
}
