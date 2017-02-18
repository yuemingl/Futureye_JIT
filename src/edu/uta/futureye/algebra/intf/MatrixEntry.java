/**
 * Copyright (c) 2010, nkliuyueming@gmail.com. All rights reserved.
 * 
 * 
 */
package edu.uta.futureye.algebra.intf;

/**
 * An entry of a sparse matrix. Returned by the iterators over a sparse matrix
 * 
 * @author liuyueming
 * 
 */
public interface MatrixEntry {
    /**
     * Return the current number of row (row index)
     */
    int getRow();

    /**
     * Return the current number of column (column index)
     */
    int getCol();

    /**
     * Return the current value
     */
    double getValue();

    /**
     * Set the current value
     */
    void setValue(double value);
}