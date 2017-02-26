/**
 * Copyright (c) 2010, nkliuyueming@gmail.com. All rights reserved.
 * 
 * 
 */
package edu.uta.futureye.lib.assembler;

import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.core.intf.FiniteElement;
import edu.uta.futureye.core.intf.VecFiniteElement;

/**
 * Parameters passed into functions in the process of assembly
 * 
 */
public class AssembleParam {
	public VecFiniteElement fe;
	public Element element; //current element which is being assemble
	public int trialDOFIdx; //local trial function index
	public int testDOFIdx;  //local test function index
	
	public Node node; //In some case, node is used instead of array of coordinates
	
	public AssembleParam(Element e, int i, int j) {
		this.element = e;
		this.trialDOFIdx = i;
		this.testDOFIdx = j;
	}
	
	public AssembleParam(Element e, int i, int j, Node n) {
		this.element = e;
		this.trialDOFIdx = i;
		this.testDOFIdx = j;
		this.node = n;
	}

}
