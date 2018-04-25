package COMMoDORE;

/**
 * Enumeration AtomInfo contains basic properties of atoms such as atomic number,
 * particle mass, equilibrium separation, etc. 
 * 
 * @author thinkens
 * on 4/4/2018
 * @since 0.1.0
 */
public enum AtomInfo {
    
    HELIUM(2,6.646632348e-27,-1,-1,-1,-1,-1),
    ALUMINUM(13,4.48039e-26,1.23e-10,6.111e-22,4.04e-10,4.04e-10,4.04e-10), //eq sep (eq. sep from Totten and Mackenzie)
    COPPER(29,1.0552308203254732e-25,-1,-1,-1,-1,-1),
    XENON(54,2.1802e-25,-1,-1,-1,-1,-1);
    
    /**
     * Atomic number, unitless.
     */
    public final int Z;
    /**
     * Mass of atom, kg.
     */
    public final double m;
    /**
     * Atomic separation at equilibrium, m.
     */
    public final double equilibriumSeparation;
    /**
     * Depth of the potential well for two atoms of the same species, J. Multiply
     * by N_A to find literature value in J/mol.
     */
    public final double homopotentialWellDepth;
    /**
     * Lattice constants, m.
     */
    public final double a, b, c;
    
    private AtomInfo(int Z, double m, double s, double homopotentialWellDepth, double a, double b, double c) {
        
        this.Z = Z;
        this.m = m;
        this.equilibriumSeparation = s;
        this.homopotentialWellDepth = homopotentialWellDepth;
        this.a = a;
        this.b = b;
        this.c = c;
    }
    
    
}
