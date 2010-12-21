package edu.uta.futureye.util;

import edu.uta.futureye.algebra.Matrix;
import edu.uta.futureye.core.Element;

public class PairElementMatrix {
	public Element element;
	public Matrix localMatrix;
	public PairElementMatrix(Element element,Matrix localMatrix) {
		this.element = element;
		this.localMatrix = localMatrix;
	}
}
