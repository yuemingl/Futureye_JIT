package edu.uta.futureye.io;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;

import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.util.container.ElementList;
import edu.uta.futureye.util.container.NodeList;

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
						Node node = new Node(index, x, y);
						mesh.addNode(node);
					} else if(mesh.getElementList().size() < nElement) {
						String type = line[2];
						if(type.equalsIgnoreCase("tri")) {
							NodeList list = new NodeList();
							list.add(mesh.getNodeList().at(Integer.valueOf(line[3])));
							list.add(mesh.getNodeList().at(Integer.valueOf(line[4])));
							list.add(mesh.getNodeList().at(Integer.valueOf(line[5])));
							Element ele = new Element(list);
							mesh.addElement(ele);
						} else if(type.equalsIgnoreCase("quad")) {
							NodeList list = new NodeList();
							list.add(mesh.getNodeList().at(Integer.valueOf(line[3])));
							list.add(mesh.getNodeList().at(Integer.valueOf(line[4])));
							list.add(mesh.getNodeList().at(Integer.valueOf(line[5])));
							list.add(mesh.getNodeList().at(Integer.valueOf(line[6])));
							Element ele = new Element(list);
							mesh.addElement(ele);						
						}
					}
				}
			}
			br.close();
			in.close();

			return mesh;
		
		} catch (Exception e) {
			e.printStackTrace();
		}
		return null;
	}
	
	public Mesh read3DMesh() {
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
				if(str.startsWith("#")) continue;
				String[] line = str.split("(\\s)+");
				if(nNode == 0) {
					//Read in node number and element number
					nNode = Integer.valueOf(line[0]);
					nElement = Integer.valueOf(line[1]);
				} else {
					if(mesh.getNodeList().size() < nNode) {
						int index = Integer.valueOf(line[0]);
						double x = Double.valueOf(line[1]);
						double y = Double.valueOf(line[2]);
						double z = Double.valueOf(line[3]);
						Node node = new Node(index, x, y, z);
						mesh.addNode(node);
					} else if(mesh.getElementList().size() < nElement) {
						String type = line[2];
						if(type.equalsIgnoreCase("tet")) {
							NodeList list = new NodeList();
							list.add(mesh.getNodeList().at(Integer.valueOf(line[3])));
							list.add(mesh.getNodeList().at(Integer.valueOf(line[4])));
							list.add(mesh.getNodeList().at(Integer.valueOf(line[5])));
							list.add(mesh.getNodeList().at(Integer.valueOf(line[6])));
							Element ele = new Element(list);
							mesh.addElement(ele);						
						}
					}
				}
			}
			br.close();
			in.close();
			
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
		//MeshReader r1 = new MeshReader("mixed.grd");
		MeshReader r1 = new MeshReader("block1.grd");
		//Mesh m = r1.read2DMesh();
		Mesh m = r1.read3DMesh();
		System.out.println("nodes read: "+m.getNodeList().size());
		System.out.println("elements read: "+m.getElementList().size());
		
	}
}
