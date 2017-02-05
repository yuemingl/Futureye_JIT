package edu.uta.futureye.lib.assembler;

import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Node;

public class AssembleParam {
	public Element element;
	public int trialDOFIdx;
	public int testDOFIdx;
	public Node node;
	
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
