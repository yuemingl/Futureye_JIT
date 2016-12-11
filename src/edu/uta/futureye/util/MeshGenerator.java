package edu.uta.futureye.util;

import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.util.container.NodeList;

public class MeshGenerator {
	/*
	 * Generate a rectangle mesh on [x0,x1]*[y0,y1] with nx points in X
	 * direction and ny points in Y direction
	 */
	public static Mesh rectangle(double x0, double x1, double y0, double y1,
			int nx, int ny) {
		Mesh mesh = new Mesh();
		double stepx = (x1 - x0) / nx;
		double stepy = (y1 - y0) / ny;
		// generate nodes
		/**
		 * nx=3, ny=3 
		 * node index: 1,2,...,9 
		 * element index: (1),(2),...,(4)
		 * 7-----8-----9 
		 * | (3) | (4) | 
		 * 4-----5-----6 
		 * | (1) | (2) | 
		 * 1-----2-----3
		 */
		for (int i = 0; i < ny; i++) {
			double y = y0 + i * stepy;
			for (int j = 0; j < nx; j++) {
				double x = x0 + j * stepx;
				int index = i * ny + j + 1;
				Node node = new Node(index, x, y);
				mesh.addNode(node);
			}
		}

		// generate elements
		NodeList nodes = mesh.getNodeList();
		
		for (int i = 0; i < nx - 1; i++) {
			for (int j = 0; j < ny - 1; j++) {
				int n1 = i * nx + j;
				int n2 = n1 + 1;
				int n3 = (i + 1) * nx + j;
				NodeList list = new NodeList();
				list.add(nodes.at(n1+1));
				list.add(nodes.at(n2+1));
				list.add(nodes.at(n3+1));
				Element e = new Element(list);
				mesh.addElement(e);

				e = new Element();
				n1 = i * nx + j + 1;
				n2 = (i + 1) * nx + j + 1;
				n3 = n2 - 1;
				list = new NodeList();
				list.add(nodes.at(n1+1));
				list.add(nodes.at(n2+1));
				list.add(nodes.at(n3+1));
				e = new Element(list);
				mesh.addElement(e);
			}
		}
		//mesh.printMeshInfo();
		return mesh;
	}

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		rectangle(0,1,0,1,3,3);
	}
}
