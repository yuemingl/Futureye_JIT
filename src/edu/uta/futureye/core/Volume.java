package edu.uta.futureye.core;

import edu.uta.futureye.core.geometry.GeoEntity3D;

/**
 * Global block of an element
 * 三维体
 * 
 * TODO 如果一个单元包含多个体该怎么办？例如：TriBilinearElement
 * 
 * @author liuyueming
 *
 */
public class Volume extends GeoEntity3D<FaceLocal,EdgeLocal,NodeLocal> {
	public Element owner;
	
	
	
}
