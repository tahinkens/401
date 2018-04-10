/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package COMMoDORE;

/**
 *
 * @author thinkens
 */
public enum AtomInfo {
    
    HELIUM(2,6.646632348e-27),
    ALUMINUM(13,4.48039e-27),
    COPPER(29,1.0552308203254732e-25),
    XENON(54,2.1802e-25);
    
    public final int Z; //atomic number
    public final double m; //mass per atom
    
    private AtomInfo(int Z, double m) {
        
        this.Z = Z;
        this.m = m;
    }
    
    
}
