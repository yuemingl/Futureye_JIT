package edu.uta.futureye.io;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;


import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.util.NodeList;

public class MeshReader {
	Mesh mesh = new Mesh();
	String file = null;
	
	public MeshReader(String fileName) {
		this.file = fileName;
	}
	
	public Mesh read2DMesh() {
		FileInputStream in;
		try {
			in = new FileInputStream(file);

			InputStreamReader reader = new InputStreamReader(in,"UTF-8");
			BufferedReader br = new BufferedReader(reader);
	
			String str = null;
			int nNode = 0;
			int nElement = 0;
			mesh.clearAll();
			while((str = br.readLine()) != null){
				System.out.println(str);				
				if(str.startsWith("#")) continue;
				String[] line = str.split("(\\s)+");
				if(nNode == 0) {
					nNode = Integer.valueOf(line[0]);
					nElement = Integer.valueOf(line[1]);
				} else {
					if(mesh.getNodeList().size() < nNode) {
						int index = Integer.valueOf(line[0]);
						double x = Double.valueOf(line[1]);
						double y = Double.valueOf(line[2]);
						Node node = new Node(2);
						node.set(index, x, y);
						mesh.addNode(node);
					} else if(mesh.getElementList().size() < nElement) {
						String type = line[2];
						if(type.equalsIgnoreCase("tri")) {
							Element ele = new Element();
							Node node;
							int idxa = Integer.valueOf(line[3]);
							node = mesh.getNodeList().at(idxa);
							ele.addNode(node,true);
							int idxb = Integer.valueOf(line[4]);
							node = mesh.getNodeList().at(idxb);
							ele.addNode(node,true);						
							int idxc = Integer.valueOf(line[5]);
							node = mesh.getNodeList().at(idxc);
							ele.addNode(node,true);
							mesh.addElement(ele);
						} else if(type.equalsIgnoreCase("quad")) {
							Element ele = new Element();
							Node node;
							int idxa = Integer.valueOf(line[3]);
							node = mesh.getNodeList().at(idxa);
							ele.addNode(node,true);
							int idxb = Integer.valueOf(line[4]);
							node = mesh.getNodeList().at(idxb);
							ele.addNode(node,true);						
							int idxc = Integer.valueOf(line[5]);
							node = mesh.getNodeList().at(idxc);
							ele.addNode(node,true);
							int idxd = Integer.valueOf(line[6]);
							node = mesh.getNodeList().at(idxd);
							ele.addNode(node,true);
							mesh.addElement(ele);						
						}
					}
				}
			}
			
			return mesh;
		
		} catch (Exception e) {
			e.printStackTrace();
		}
		return null;
	}
	
	public NodeList getNodeList() {
		return null;
	}
	
	public static void main(String[] args) {
		MeshReader r1 = new MeshReader("mixed.grd");
		Mesh m = r1.read2DMesh();
		System.out.println("nodes read: "+m.getNodeList().size());
		System.out.println("elements read: "+m.getElementList().size());
		
	}
}
