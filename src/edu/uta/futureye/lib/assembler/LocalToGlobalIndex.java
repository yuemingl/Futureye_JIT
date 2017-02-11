package edu.uta.futureye.lib.assembler;

import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.core.intf.FiniteElement;
import edu.uta.futureye.core.intf.VecFiniteElement;

/**
 * A class for the calculation of local index to global index
 * of DOF
 * @author yueming.liu
 *
 */
public class LocalToGlobalIndex {
	Mesh mesh;
	FiniteElement fe;
	public LocalToGlobalIndex(FiniteElement fe, Mesh mesh) {
		this.mesh = mesh;
		this.fe = fe;
	}
	
	public LocalToGlobalIndex(VecFiniteElement fe, Mesh mesh) {
		this.mesh = mesh;
		this.fe = fe; //Do we need a interface for FiniteElement and VecFiniteElement?
	}
	
	public int getGlobalIndex(int idxElement, int localIndex) {
		Element e = mesh.at(idxElement);
		Node node = e.nodes.at(localIndex);
		
		fe.getNumberOfDOFs()
	}
}
