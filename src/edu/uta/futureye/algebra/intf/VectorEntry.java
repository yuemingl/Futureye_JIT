package edu.uta.futureye.algebra.intf;

/**
 * An entry of a sparse vector. Returned by the iterators over a sparse vector
 * 
 * @author liuyueming
 * 
 */
public interface VectorEntry {
    /**
     * Return the current vector index
     */
    int getIndex();

    /**
     * Return the current value
     */
    double getValue();

    /**
     * Set the current value
     */
    void setValue(double value);
}