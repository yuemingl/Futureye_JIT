package edu.uta.futureye.core;

import java.util.HashMap;
import java.util.Map;

import edu.uta.futureye.util.FutureyeException;
import edu.uta.futureye.util.container.ElementList;
import edu.uta.futureye.util.container.NodeList;
import edu.uta.futureye.util.container.ObjList;
import edu.uta.futureye.util.container.VertexList;

public class Refiner {

	/**
	 * 检查是否需要链式加密周围网格，当进行第二次加密单元时，为了保持一条边界上只有一个hanging node，
	 * 可能需要加密当前单元周围的单元
	 * @param mesh
	 * @param eToRefine
	 */
	public static ElementList checkNeighborRefinement(ElementList eToRefine) {
		Map<Integer, Element> map = new HashMap<Integer,Element>();
		ElementList findNeighbors = new ElementList();
		
		for(int i=1;i<=eToRefine.size();i++) {
			Element e = eToRefine.at(i);
			ElementList eNeighbors = e.neighbors;
			for(int j=1;j<=eNeighbors.size();j++) {
				Element findNeighbor = eNeighbors.at(j);
				if(e.getLevel() - findNeighbor.getLevel()>=1) {
					map.put(findNeighbor.globalIndex, findNeighbor);
				}
			}
//有问题，还没标记level，无法判断
//			//如果一个单元相邻单元都细化了，该单元自动细化
//			for(int j=1;j<=eNeighbors.size();j++) {
//				Element findNeighbor = eNeighbors.at(j);
//				ElementList findNbrNbr = findNeighbor.neighbors;
//				for(int k=1;k<=findNbrNbr.size();k++) {
//					boolean all = true;
//					if(findNeighbor.getLevel() >= findNbrNbr.at(k).getLevel()) {
//						all = false;
//						break;
//					}
//					if(all)
//						map.put(findNeighbor.globalIndex, findNeighbor);
//				}
//			}
		}
		for(Element e : map.values())
			findNeighbors.add(e);
		return findNeighbors;
	}
	
	public static boolean directRefine(Mesh mesh, ElementList eToRefine) {
		ElementList oldEleList = mesh.getElementList();
		int refinedNum = 0;
		
		for(int i=1;i<=eToRefine.size();i++) {
			Element e = eToRefine.at(i);
			//防止重复refine
			if(e.isRefined()) {
				refinedNum++;
				continue;
			}
			
			ElementList newEleList = new ElementList();
			VertexList vList = e.vertices();
			
			if(vList.size() == 3) {
				/*
				 * 3
				 * | \
				 * |  \
				 * 6---5
				 * |\  |\
				 * | \ | \
				 * 1--4---2
				 */
				Vertex v1 = vList.at(1);
				Vertex v2 = vList.at(2);
				Vertex v3 = vList.at(3);
				
				Node node4 = new NodeRefined(2);
				node4.setCoord(1, (v1.coord(1)+v2.coord(1))/2.0);
				node4.setCoord(2, (v1.coord(2)+v2.coord(2))/2.0);
				node4.setRefineLevel(e.getLevel()+1);
				Node tmpNode = mesh.findNode(node4);
				if(tmpNode != null) node4 = tmpNode;
				
				Node node5 = new NodeRefined(2);
				node5.setCoord(1, (v2.coord(1)+v3.coord(1))/2.0);
				node5.setCoord(2, (v2.coord(2)+v3.coord(2))/2.0);
				node5.setRefineLevel(e.getLevel()+1);
				tmpNode = mesh.findNode(node5);
				if(tmpNode != null) node5 = tmpNode;
				
				Node node6 = new NodeRefined(2);
				node6.setCoord(1, (v3.coord(1)+v1.coord(1))/2.0);
				node6.setCoord(2, (v3.coord(2)+v1.coord(2))/2.0);
				node6.setRefineLevel(e.getLevel()+1);
				tmpNode = mesh.findNode(node6);
				if(tmpNode != null) node6 = tmpNode;
				
				Node node1 = e.nodes.at(v1.localIndex);
				Node node2 = e.nodes.at(v2.localIndex);
				Node node3 = e.nodes.at(v3.localIndex);
				
				NodeList nList = new NodeList();
				
				nList.clear();
				nList.add(node1);
				nList.add(node4);
				nList.add(node6);
				Element e1 = new Element(nList);
				e1.parent = e;
				e1.setLevel(e.getLevel()+1);
				
				nList.clear();
				nList.add(node4);
				nList.add(node2);
				nList.add(node5);
				Element e2 = new Element(nList);
				e2.parent = e;
				e2.setLevel(e.getLevel()+1);

				nList.clear();
				nList.add(node5);
				nList.add(node3);
				nList.add(node6);
				Element e3 = new Element(nList);
				e3.parent = e;
				e3.setLevel(e.getLevel()+1);

				nList.clear();
				nList.add(node4);
				nList.add(node5);
				nList.add(node6);
				Element e4 = new Element(nList);
				e4.parent = e;
				e4.setLevel(e.getLevel()+1);

				newEleList.add(e1);
				newEleList.add(e2);
				newEleList.add(e3);
				newEleList.add(e4);
			} else if(vList.size() == 4) {
				/*
				 * 	4--7--3
				 *  |  |  |
				 *  8--9--6
				 *  |  |  |
				 *  1--5--2
				 */
				Vertex v1 = vList.at(1);
				Vertex v2 = vList.at(2);
				Vertex v3 = vList.at(3);
				Vertex v4 = vList.at(4);
				
				Node node5 = new NodeRefined(2);
				node5.setCoord(1, (v1.coord(1)+v2.coord(1))/2.0);
				node5.setCoord(2, (v1.coord(2)+v2.coord(2))/2.0);
				node5.setRefineLevel(e.getLevel()+1);
				Node tmpNode = mesh.findNode(node5);
				if(tmpNode != null) node5 = tmpNode;
				
				Node node6 = new NodeRefined(2);
				node6.setCoord(1, (v2.coord(1)+v3.coord(1))/2.0);
				node6.setCoord(2, (v2.coord(2)+v3.coord(2))/2.0);
				node6.setRefineLevel(e.getLevel()+1);
				tmpNode = mesh.findNode(node6);
				if(tmpNode != null) node6 = tmpNode;
				
				Node node7 = new NodeRefined(2);
				node7.setCoord(1, (v3.coord(1)+v4.coord(1))/2.0);
				node7.setCoord(2, (v3.coord(2)+v4.coord(2))/2.0);
				node7.setRefineLevel(e.getLevel()+1);
				tmpNode = mesh.findNode(node7);
				if(tmpNode != null) node7 = tmpNode;
				
				Node node8 = new NodeRefined(2);
				node8.setCoord(1, (v4.coord(1)+v1.coord(1))/2.0);
				node8.setCoord(2, (v4.coord(2)+v1.coord(2))/2.0);
				node8.setRefineLevel(e.getLevel()+1);
				tmpNode = mesh.findNode(node8);
				if(tmpNode != null) node8 = tmpNode;
				
				Node node9 = new NodeRefined(2);
				node9.setCoord(1, (node5.coord(1)+node7.coord(1))/2.0);
				node9.setCoord(2, (node5.coord(2)+node7.coord(2))/2.0);
				node9.setRefineLevel(e.getLevel()+1);
				tmpNode = mesh.findNode(node9);
				if(tmpNode != null) node9 = tmpNode;
				
				Node node1 = e.nodes.at(v1.localIndex);
				Node node2 = e.nodes.at(v2.localIndex);
				Node node3 = e.nodes.at(v3.localIndex);
				Node node4 = e.nodes.at(v4.localIndex);
				
				NodeList nList = new NodeList();
				
				nList.clear();
				nList.add(node1);
				nList.add(node5);
				nList.add(node9);
				nList.add(node8);
				Element e1 = new Element(nList);
				e1.parent = e;
				e1.setLevel(e.getLevel()+1);
				
				nList.clear();
				nList.add(node5);
				nList.add(node2);
				nList.add(node6);
				nList.add(node9);
				Element e2 = new Element(nList);
				e2.parent = e;
				e2.setLevel(e.getLevel()+1);

				nList.clear();
				nList.add(node9);
				nList.add(node6);
				nList.add(node3);
				nList.add(node7);
				Element e3 = new Element(nList);
				e3.parent = e;
				e3.setLevel(e.getLevel()+1);

				nList.clear();
				nList.add(node8);
				nList.add(node9);
				nList.add(node7);
				nList.add(node4);
				Element e4 = new Element(nList);
				e4.parent = e;
				e4.setLevel(e.getLevel()+1);

				newEleList.add(e1);
				newEleList.add(e2);
				newEleList.add(e3);
				newEleList.add(e4);
			} else {
				FutureyeException ex = new FutureyeException("Unsupported element type for refinement!");
				ex.printStackTrace();
			}
			
			NodeList oldNodeList = mesh.getNodeList();
			
			//单元编号有问题，如果两次加密，会产生重复的单元编号
			//解决方案：每细化一个单元，整个单元列表都重新编号
			//结点列表由于只有增加结点，因此保持原方案
//			for(int j=1;j<=eList.size();j++) {
//				Element eNew = eList.at(j);
//				if(eNew.globalIndex == 0) {
//					eNew.globalIndex = oldEleList.size()+1;
//					oldEleList.add(eNew);
//					for(int k=1;k<=eNew.nodes.size();k++) {
//						Node nNew = eNew.nodes.at(k);
//						if(nNew.globalIndex == 0) {
//							nNew.globalIndex = oldNodeList.size()+1;
//							oldNodeList.add(nNew);
//						}
//					}
//				}
//			}
			
			//细化新增的单元与老单元放入同一个list，按照在列表中存储索引，统一重新编号
			//注意：此时被细化的父单元还没有删除
			//此时将新结点加入老结点列表，细化后，不需要删除多余结点
			oldEleList.addAll(newEleList);
			for(int j=1;j<=oldEleList.size();j++) {
				Element ele = oldEleList.at(j);
				if(ele.globalIndex == 0) {
					ele.globalIndex = j;
					for(int k=1;k<=ele.nodes.size();k++) {
						Node newNode = ele.nodes.at(k);
						if(newNode.globalIndex == 0) {
							newNode.globalIndex = oldNodeList.size()+1;
							oldNodeList.add(newNode);
						}
					}
				}  else {
					ele.globalIndex = j;
				}
			}
			e.childs = newEleList;
		}

		return refinedNum==eToRefine.size();
	}

	
	public static void reNumberElements(Mesh mesh) {
		ElementList eList = mesh.getElementList();
		for(int j=1;j<=eList.size();j++) {
			Element ele = eList.at(j);
			ele.globalIndex = j;
		}
	}
	
	public static void refineOnce(Mesh mesh, ElementList eToRefine) {
		ElementList eList = mesh.getElementList();
		while(true) {
			ElementList eNeighbors = checkNeighborRefinement(eToRefine);
			if(eNeighbors.size()>0) {
				boolean stop = directRefine(mesh,eNeighbors);
				computeHangingNode(eNeighbors);
				for(int iToRe=1;iToRe<=eNeighbors.size();iToRe++) {
					eList.remove(eNeighbors.at(iToRe));
				}
				reNumberElements(mesh);
				mesh.computeNodeBelongsToElements();
				mesh.computeNeighborNodes();
				mesh.computeGlobalEdge();
				mesh.computeNeighborElements();
				if(stop) break;
			} else {
				directRefine(mesh,eToRefine);
				computeHangingNode(eToRefine);
				for(int iToRe=1;iToRe<=eToRefine.size();iToRe++) {
					eList.remove(eToRefine.at(iToRe));
				}
				reNumberElements(mesh);
				mesh.computeNodeBelongsToElements();
				mesh.computeNeighborNodes();
				mesh.computeGlobalEdge();
				mesh.computeNeighborElements();
				break;
			}
		}
	}
	
	//计算hanging node
	public static void computeHangingNode(ElementList eToRefine) {
		for(int iToRe=1;iToRe<=eToRefine.size();iToRe++) {
			ElementList eChilds = eToRefine.at(iToRe).childs;
			for(int i=1;i<=eChilds.size();i++) {
				Element eChild = eChilds.at(i);
				for(int j=1;j<=eChild.nodes.size();j++) {
					Node nNew = eChild.nodes.at(j);
					//!!! Node的level与单元level相同时，清空hanging node限制值条件，准备重新计算
					if(nNew.getRefineLevel() == eChild.getLevel()) {
						NodeRefined nRefined = (NodeRefined)nNew;
						nRefined.clearConstrainNodes();
					}
				}
			}
		}
		
		for(int iToRe=1;iToRe<=eToRefine.size();iToRe++) {
			ElementList eChilds = eToRefine.at(iToRe).childs;
			ElementList eNeighbors = eToRefine.at(iToRe).neighbors;
			//循环大单元中的每个子单元
			for(int i=1;i<=eChilds.size();i++) {
				Element eChild = eChilds.at(i);
				//循环每个子单元的结点
				for(int j=1;j<=eChild.nodes.size();j++) {
					Node nNew = eChild.nodes.at(j);
					if(nNew.getRefineLevel() == eChild.getLevel()) {
						NodeRefined nRefined = (NodeRefined)nNew;
						//循环大单元的相邻单元，判断该节点是否Hanging node，
						//如果是边界单元的边界上加密，nRefined将没有ConstrainNode，默认为非Hanging node
						for(int k=1;k<=eNeighbors.size();k++) {
							Element eNeighbor = eNeighbors.at(k);
							//!!! 相邻单元没有标记为加密  并且 相邻单元的层次要低于该单元
							if(!eNeighbor.isRefined() && eNeighbor.getLevel()<eChild.getLevel()) {
								ObjList<EdgeLocal> edges = eNeighbors.at(k).edges();
								for(int kk=1;kk<=edges.size();kk++) {
									//判断结点是否在大单元边上
									if(edges.at(kk).globalEdge.isCoordOnEdge(nRefined.coords())) {
										NodeList endNodes = edges.at(kk).globalEdge.getEndNodes();
										//System.out.println("Hanging node:"+nRefined.globalIndex+
										//		" on "+endNodes.at(1).globalIndex+" "+endNodes.at(2).globalIndex);
										nRefined.addConstrainNode(endNodes.at(1));
										nRefined.addConstrainNode(endNodes.at(2));
									}
								}
							}
						}
					}
				}
			}
		}
	}

}
