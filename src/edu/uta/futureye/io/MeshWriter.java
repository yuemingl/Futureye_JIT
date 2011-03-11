package edu.uta.futureye.io;

import java.io.File;
import java.io.FileOutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;

import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.util.list.ElementList;
import edu.uta.futureye.util.list.NodeList;

public class MeshWriter {
	Mesh mesh = null;
	
	public MeshWriter(Mesh mesh) {
		this.mesh = mesh;
	}
	
	public void writeTechplot(String fileName, Vector u) {
		FileOutputStream out;
		try {
			File file = new File(fileName);
			out = new FileOutputStream(file);
			OutputStreamWriter writer = new OutputStreamWriter(out, "UTF-8");
			PrintWriter br = new PrintWriter(writer);
			
			NodeList nodes = mesh.getNodeList();
			ElementList elements = mesh.getElementList();
			int nNode = nodes.size();
			int nElement = elements.size();
			//找到包含结点最多的单元的结点数目，因为有些网格是混合单元网格
			int nMaxNodes = 0;
			for(int i=1;i<=elements.size();i++) {
				if(elements.at(i).nodes.size() > nMaxNodes)
					nMaxNodes = elements.at(i).nodes.size();
			}
			int dim = nodes.at(1).dim();
			
			if(dim == 2) { //二维单元
				if(nMaxNodes % 3 == 0) {
					br.println("VARIABLES=\"X\",\"Y\",\"U\"");
					
					if(nMaxNodes == 3)
						br.println(String.format("ZONE F=FEPOINT ET=TRIANGLE N=%d E=%d",nNode,nElement));
					else if(nMaxNodes == 6)
						br.println(String.format("ZONE F=FEPOINT ET=TRIANGLE N=%d E=%d",nNode,4*nElement));
						
					for(int i=1;i<=nNode;i++) {
						Node node = nodes.at(i);
						br.println(String.format("%f    %f    %f", 
								node.coord(1),
								node.coord(2),
								u.get(i)));			
					}
					for(int i=1;i<=nElement;i++) {
						Element e = elements.at(i);
						if(e.nodes.size() == 3) {
							br.println(String.format("%d    %d    %d", 
									e.nodes.at(1).globalIndex,
									e.nodes.at(2).globalIndex,
									e.nodes.at(3).globalIndex
									));
						} else if(e.nodes.size() == 6) {
							br.println(String.format("%d    %d    %d", 
									e.nodes.at(1).globalIndex,
									e.nodes.at(4).globalIndex,
									e.nodes.at(6).globalIndex
									));						
							br.println(String.format("%d    %d    %d", 
									e.nodes.at(2).globalIndex,
									e.nodes.at(5).globalIndex,
									e.nodes.at(4).globalIndex
									));	
							br.println(String.format("%d    %d    %d", 
									e.nodes.at(3).globalIndex,
									e.nodes.at(6).globalIndex,
									e.nodes.at(5).globalIndex
									));	
							br.println(String.format("%d    %d    %d", 
									e.nodes.at(4).globalIndex,
									e.nodes.at(5).globalIndex,
									e.nodes.at(6).globalIndex
									));	
						} else {
							System.out.println("Error: TRIANGLE nodes number="+e.nodes.size());
						}
						
					}
				} else if(nMaxNodes % 4 == 0) {
					br.println("VARIABLES=\"X\",\"Y\",\"U\"");
					
					if(nMaxNodes == 4)
						br.println(String.format("ZONE F=FEPOINT ET=QUADRILATERAL N=%d E=%d",nNode,nElement));
					else if(nMaxNodes == 8)
						br.println(String.format("ZONE F=FEPOINT ET=QUADRILATERAL N=%d E=%d",nNode,5*nElement));
	
					for(int i=1;i<=nNode;i++) {
						Node node = nodes.at(i);
						br.println(String.format("%f    %f    %f", 
								node.coord(1),
								node.coord(2),
								u.get(i)));			
					}
					for(int i=1;i<=nElement;i++) {
						Element e = elements.at(i);
						if(e.nodes.size() == 4) {
							br.println(String.format("%d    %d    %d    %d", 
									e.nodes.at(1).globalIndex,
									e.nodes.at(2).globalIndex,
									e.nodes.at(3).globalIndex,
									e.nodes.at(4).globalIndex
									));
						} else if(e.nodes.size() == 8) {
							br.println(String.format("%d    %d    %d    %d", 
									e.nodes.at(1).globalIndex,
									e.nodes.at(5).globalIndex,
									e.nodes.at(8).globalIndex,
									e.nodes.at(1).globalIndex
									));						
							br.println(String.format("%d    %d    %d    %d", 
									e.nodes.at(2).globalIndex,
									e.nodes.at(6).globalIndex,
									e.nodes.at(5).globalIndex,
									e.nodes.at(2).globalIndex
									));						
							br.println(String.format("%d    %d    %d    %d", 
									e.nodes.at(3).globalIndex,
									e.nodes.at(7).globalIndex,
									e.nodes.at(6).globalIndex,
									e.nodes.at(3).globalIndex
									));
							br.println(String.format("%d    %d    %d    %d", 
									e.nodes.at(4).globalIndex,
									e.nodes.at(8).globalIndex,
									e.nodes.at(7).globalIndex,
									e.nodes.at(4).globalIndex
									));
							br.println(String.format("%d    %d    %d    %d", 
									e.nodes.at(5).globalIndex,
									e.nodes.at(6).globalIndex,
									e.nodes.at(7).globalIndex,
									e.nodes.at(8).globalIndex
									));
						} else if(e.nodes.size() == 3) {
							br.println(String.format("%d    %d    %d    %d", 
									e.nodes.at(1).globalIndex,
									e.nodes.at(2).globalIndex,
									e.nodes.at(3).globalIndex,
									e.nodes.at(1).globalIndex
									));
						} else {
							System.out.println("Error: QUADRILATERAL nodes number="+e.nodes.size());
						}
						
					}				
				}
			} else if(dim == 3) { //三维单元
				br.println("VARIABLES=\"X\",\"Y\",\"Z\",\"U\"");
				//四面体单元
				if(nMaxNodes == 4)
					br.println(String.format("ZONE F=FEPOINT ET=TETRAHEDRON N=%d E=%d",nNode,nElement));
				for(int i=1;i<=nNode;i++) {
					Node node = nodes.at(i);
					br.println(String.format("%f    %f    %f    %f", 
							node.coord(1),
							node.coord(2),
							node.coord(3),
							u.get(i)));			
				}
				for(int i=1;i<=nElement;i++) {
					Element e = elements.at(i);
					if(e.nodes.size() == 4) {
						br.println(String.format("%d    %d    %d    %d", 
								e.nodes.at(1).globalIndex,
								e.nodes.at(2).globalIndex,
								e.nodes.at(3).globalIndex,
								e.nodes.at(4).globalIndex
								));
					}
				}
				
			}
			br.close();
			out.close();
			
		
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	
	public void writeTechplotLine(String fileName, Vector x, Vector y) {
		FileOutputStream out;
		try {
			File file = new File(fileName);
			out = new FileOutputStream(file);
			OutputStreamWriter writer = new OutputStreamWriter(out, "UTF-8");
			PrintWriter br = new PrintWriter(writer);
			
			for(int i=1;i<=x.getDim();i++)
				br.println(x.get(i)+"\t"+y.get(i));
				
			br.close();
			out.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
