package edu.uta.futureye.io;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;

import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.util.container.ElementList;
import edu.uta.futureye.util.container.NodeList;

/**
 * Read .node and .ele files from the output of 
 * Triangle (see https://www.cs.cmu.edu/~quake/triangle.html)
 * 
 * 
 * @author yuemingl
 *
 */
public class MeshReaderForTriangle {
	
	Mesh mesh = new Mesh();
	String nodeFile = null;
	String eleFile = null;
	
	public boolean debug = false;
	
	public MeshReaderForTriangle(String nodeFile, String eleFile) {
		this.nodeFile = nodeFile;
		this.eleFile = eleFile;
	}
	
	public Mesh read2DMesh() {
		FileInputStream in;
		try {
			in = new FileInputStream(nodeFile);

			InputStreamReader reader = new InputStreamReader(in,"UTF-8");
			BufferedReader br = new BufferedReader(reader);
	
			String str = null;
			int nNode = 0;
			int nElement = 0;
			mesh.clearAll();
			while((str = br.readLine()) != null){
				if(debug)
					System.out.println(str);
				if(str.startsWith("#")) continue;
				String[] word = str.trim().split("(\\s)+");
				if(nNode == 0) {
					nNode = Integer.valueOf(word[0]);
				} else if(mesh.getNodeList().size() < nNode) {
						int index = Integer.valueOf(word[0]);
						double x = Double.valueOf(word[1]);
						double y = Double.valueOf(word[2]);
						Node node = new Node(index, x, y);
						mesh.addNode(node);
				} else
					break;
			}
			br.close();
			in.close();
			
			in = new FileInputStream(eleFile);
			reader = new InputStreamReader(in,"UTF-8");
			br = new BufferedReader(reader);
			while((str = br.readLine()) != null){
				if(debug)
					System.out.println(str);
				if(str.startsWith("#")) continue;
				String[] word = str.trim().split("(\\s)+");
				if(nElement == 0) {
					nElement = Integer.valueOf(word[0]);
				} else if(mesh.getElementList().size() < nElement) {
					NodeList list = new NodeList();
					//index: word[0]
					list.add(mesh.getNodeList().at(Integer.valueOf(word[1])));
					list.add(mesh.getNodeList().at(Integer.valueOf(word[2])));
					list.add(mesh.getNodeList().at(Integer.valueOf(word[3])));
					Element ele = new Element(list);
					mesh.addElement(ele);
				} else
					break;
			}
			br.close();
			in.close();

			ElementList nEList = mesh.getElementList();
			int nE = nEList.size();
			for(int i=1;i<=nE;i++) {
				Element e = nEList.at(i);
				e.adjustVerticeToCounterClockwise();
			}
			return mesh;
		
		} catch (Exception e) {
			e.printStackTrace();
		}
		return null;
	}
	
	
	public NodeList getNodeList() {
		return mesh.getNodeList();
	}
	
	public ElementList getElementList() {
		return mesh.getElementList();
	}
	
	public static void main(String[] args) {
		MeshReaderForTriangle r1 = new MeshReaderForTriangle("./iphone/hand1.1.node","./iphone/hand1.1.ele");
		Mesh m = r1.read2DMesh();
		m.computeNodeBelongsToElements();
		System.out.println("nodes read: "+m.getNodeList().size());
		System.out.println("elements read: "+m.getElementList().size());
	}

}
